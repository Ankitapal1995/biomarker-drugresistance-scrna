# biomarker-drugresistance-scrna
Single-Cell RNA-Seq Analysis of Drug Resistance Biomarkers in Breast Cancer
Single-Cell RNA-Seq Analysis of Drug Resistance Biomarkers in Breast Cancer
Project Overview

This repository contains the scripts, documentation, and key results for a project focused on the identification of cell type-specific biomarkers for drug resistance in breast cancer using single-cell RNA-sequencing (scRNA-seq) data.
A total of 40 human breast cancer samples (20 drug-resistant and 20 susceptible) were analyzed to uncover gene expression signatures and enriched pathways associated with drug response and resistance.
Repository Structure


Motivation

Drug resistance in breast cancer is a major challenge in precision medicine. This project aims to dissect the cellular and molecular features distinguishing drug-resistant from drug-sensitive tumor cells by leveraging high-resolution single-cell transcriptomics and robust gene set enrichment analysis.
Methods

    scRNA-seq Pipeline:
    Reads underwent QC, cell type annotation using CellTypist, and marker gene identification for each cell type.

    Gene Set Enrichment Analysis (GSEA):
    Marker gene lists from resistant and susceptible samples were compared across hallmark gene sets to identify enriched pathways and gene signatures.

    Visualization:
    Resulting figures, including pathway and expression plots, support interpretation and presentation.

Key Results
Drug-Sensitive/Susceptible Samples

    Enriched Pathways: MYC and E2F target activation, oxidative phosphorylation, cell cycle progression, mTORC1 signaling.

    Marker Genes: MYC, NPM1, CCNB1, MCM5, ATP5B, CDC45, BRCA1.

Drug-Resistant Samples

    Enriched Pathways: TNF-alpha/NFkB, inflammatory response, JAK/STAT signaling, EMT, hypoxia, apoptosis resistance.

    Marker Genes: NFKBIA, IL6, STAT3, SNAI1, RELA, BCL2, CXCL8, HIF1A.

Refer to the enrichment analysis HTML reports in the results/ directory and related figures in figures/ for complete details.
How to Run the Analysis

    Requirements:

        Python >= 3.8

        Required packages listed in the script (Scanpy, CellTypist, Pandas, Numpy, etc.)

    Run the Pipeline:

    bash
    python scrna_40sample_pipeline_v3.py

    Configure paths and parameters as needed within the script.

    Examine Results:

        Check processed CSVs and enrichment analysis reports in results/

        Open widgets and figures in figures/ for pathway visualization

Figures

Example figures to include (place in figures/):

    UMAP/tSNE cluster plots

    Marker gene heatmaps

    Pathway enrichment barplots (from GSEA)

    Differential gene expression violin/boxplots

Reproducibility

All relevant scripts and processed data summaries are included. For full data replicability, provide necessary data files (if permitted) or link to public sources.

Contact

Principal Investigator: Ankita Pal
