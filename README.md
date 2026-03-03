Sodium Valproate Induces Pancreatic Injury
This repository contains the code and data analysis workflows for the study investigating the mechanisms of pancreatic injury induced by Sodium Valproate (VPA).

📖 Project Overview
This project aims to elucidate the molecular mechanisms underlying Sodium Valproate-induced pancreatic damage through integrated multi-omics analysis. The study combines transcriptomic and metabolomic approaches to provide a comprehensive view of the pathological process.

Key analytical approaches include:

RNA-Seq Analysis: To identify differentially expressed genes and key biological pathways.

Pseudo-targeted Metabolomic Analysis: To profile metabolic changes and identify perturbed metabolites.

📁 Repository Structure
This repository is organized as follows:

RNA-SEQ analysis 1st.R

Primary R script for processing and analyzing RNA sequencing data. This script likely includes steps for quality control, alignment, quantification, and initial differential expression analysis.

RNA-SEQ analysis 2nd.R

Secondary R script for downstream RNA-seq analysis. This may involve more detailed pathway analysis (e.g., GO, KEGG), visualization (e.g., heatmaps, volcano plots), and integration with other data types.

pseudo-targeted metabolomic analysis.R

R script for analyzing pseudo-targeted metabolomics data. This script handles data preprocessing, identification of differential metabolites, and metabolic pathway enrichment analysis.

README.md

This file, providing an overview of the project.

🚀 Getting Started
Prerequisites
The analysis scripts are written in R. To run them, you will need to have R installed on your system. It is highly recommended to use RStudio for an integrated development environment.

Usage
Clone the repository:

bash
git clone https://github.com/Pancreatologist/Sodium-valproate-induces-pancreatic-injury.git
Prepare your data: Ensure your input data files (e.g., count matrices, metabolite tables) are placed in the appropriate directories as expected by the scripts. You may need to modify the file paths within the .R scripts.

Run the scripts: Open and execute the R scripts in RStudio or via the R command line. It is advisable to run them in the logical order of analysis.

📊 Data
The raw and processed data supporting this study (e.g., RNA-seq fastq files, metabolomic spectral data) are likely deposited in a public repository (e.g., GEO, MetaboLights). Please check the paper or后续的补充材料 for accession numbers.

Note: For now, input data formats are as expected by the provided R scripts.

📝 License
This project is licensed under the terms of the MIT license. See the LICENSE file for details.

📧 Contact
Project Author: Dr. Wenhao Cai and Dr.Di Wu
GitHub: @Pancreatologist

For any questions or issues regarding the code, please feel free to open an issue on this repository.
