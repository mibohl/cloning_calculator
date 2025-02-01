import streamlit as st
import numpy as np

def calculate_dna_volumes(lengths, concentrations, ratio, n_total):
    """
    Calculate the required DNA volumes for a molecular cloning reaction.
    
    Parameters:
    lengths (list): List of fragment lengths (kbp) including the backbone.
    concentrations (list): List of DNA concentrations (ng/µL) for each fragment.
    ratio (float): Insert-to-backbone molar ratio.
    n_total (float): Total amount of DNA in the reaction (pmol).
    
    Returns:
    tuple: Arrays of volumes (µL), weights (ng), and molar amounts (pmol) for each fragment.
    """
    # Calculate the number of pmol for the backbone
    n_backbone = n_total / (1 + ratio * (len(lengths) - 1))
    
    vol = np.zeros(len(lengths))
    weight = np.zeros(len(lengths))
    n = np.zeros(len(lengths))

    for i in range(len(lengths)):
        # Calculate pmol for inserts based on ratio
        n[i] = ratio * n_backbone if i > 0 else n_backbone
        # Convert pmol to weight in ng (10^-12 = conversion factor from pmol to mol,
        # 660 = average weight of a bp,
        # 10^9 = conversion factor from g to ng,
        # 10^3 = conversion factor from bp to kbp)
        weight[i] = n[i] * 10**-12 * 660 * 10**9 * lengths[i] * 10**3
        # Calculate volume (µL) based on concentration (ng/µL)
        vol[i] = weight[i] / concentrations[i]
    
    return vol, weight, n
       
# Streamlit UI setup
st.set_page_config(initial_sidebar_state='expanded')
st.title("Molecular Cloning Calculator")
st.write("This app calculates DNA fragment volumes for molecular cloning reactions such as Gibson or Golden Gate assemblies.")

# General reaction inputs
st.sidebar.header("Input Parameters")
num_inserts = st.sidebar.number_input("Number of inserts", min_value=1, value=1)
ratio = st.sidebar.number_input("Insert : Backbone ratio", min_value=0.5, value=2.0, step=0.5)
n_total = st.sidebar.number_input("Total amount of DNA in reaction (pmol)", min_value=0.05, value=0.5, step=0.01)

# Backbone input
st.header("DNA Fragment Information")
col1, col2 = st.columns(2)
with col1:
    bb_length = st.number_input("Length of the backbone (kbp)", min_value=0.1, value=5.0, step=0.1)
with col2:
    bb_concentration = st.number_input("Concentration of the backbone (ng/µL)", min_value=0.1, value=100.0, step=0.1)

# Store lengths and concentrations
lengths = [bb_length]
concentrations = [bb_concentration]

# Spacer for UI readability
st.write("")

# Insert inputs
for i in range(num_inserts):
    col1, col2 = st.columns(2)
    with col1:
        length = st.number_input(f"Length of insert {i+1} (kbp)", min_value=0.1, value=1.0, step=0.1, key=f"length_{i}")
        lengths.append(length)
    with col2:
        concentration = st.number_input(f"Concentration of insert {i+1} (ng/µL)", min_value=0.1, value=150.0, step=0.1, key=f"conc_{i}")
        concentrations.append(concentration)

# Calculation button
if st.button("Calculate"):
    vol, weight, n = calculate_dna_volumes(lengths, concentrations, ratio, n_total)
    
    # Backbone results
    st.subheader("Backbone")
    col1, col2, col3 = st.columns(3)
    col1.metric("Amount (pmol)", f"{n[0]:.2f}")
    col2.metric("Weight (ng)", f"{weight[0]:.2f}")
    col3.metric("Volume (µL)", f"{vol[0]:.2f}")

    # Insert results
    for i in range(1, len(lengths)):
        st.subheader(f"DNA Fragment {i}")
        col1, col2, col3 = st.columns(3)
        col1.metric("Amount (pmol)", f"{n[i]:.2f}")
        col2.metric("Weight (ng)", f"{weight[i]:.2f}")
        col3.metric("Volume (µL)", f"{vol[i]:.2f}")
    
    # Total volume
    total_volume = np.sum(vol)
    st.subheader("Total Volume")
    st.metric("Total DNA solution volume in reaction (µL)", f"{total_volume:.2f}")