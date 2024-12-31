#!/usr/bin/env python
# coding: utf-8

# In[6]:

pip install scanpy
import streamlit as st
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import leidenalg
from sklearn.decomposition import PCA
import os

# Streamlit settings
mpl.rcParams['figure.dpi'] = 70
sc.settings.verbosity = 0
sc.settings.set_figure_params(dpi=80, facecolor="white", frameon=False)

# Sidebar
st.sidebar.title("Scanpy Clustering App")
uploaded_file = st.sidebar.file_uploader("Upload .h5ad file", type="h5ad")
resolution = st.sidebar.slider("Leiden Resolution", 0.1, 2.0, step=0.1, value=0.5)
run_analysis = st.sidebar.button("Run Clustering")

# Main content
st.title("Scanpy Clustering Visualization")
if uploaded_file:
    st.write("## Input File Details")
    st.write(f"Filename: {uploaded_file.name}")
    file_path = f"./{uploaded_file.name}"
    with open(file_path, "wb") as f:
        f.write(uploaded_file.getbuffer())
    
    # Run analysis
    if run_analysis:
        st.write("### Running Clustering...")
        adata = sc.read(file_path)

        # Preprocessing
        sc.pp.neighbors(adata, n_pcs=30)
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=resolution, key_added=f"leiden_res{resolution}")
        
        # Plotting
        fig, ax = plt.subplots()
        sc.pl.umap(
            adata,
            color=[f"leiden_res{resolution}"],
            legend_loc="on data",
            show=False,
            ax=ax
        )
        st.pyplot(fig)

        # Save the output
        output_path = "Scanpy_analysis.h5ad"
        adata.write_h5ad(output_path)
        st.success(f"Analysis complete. File saved as {output_path}.")
        with open(output_path, "rb") as f:
            st.download_button(
                label="Download Processed File",
                data=f,
                file_name=output_path,
                mime="application/octet-stream"
            )
else:
    st.write("Please upload a .h5ad file to start the analysis.")


# In[ ]:




