# Lung Cancer Immunotherapy Analysis

# Data anlaysis for manuscript 

# Prerequisites
1. MacOS/Windows/Linux
2. R
3. Python

# Installing R [typical time = 5-10 minutes]
1. Open an internet browser and go to www.r-project.org.
2. Click the "download R" link in the middle of the page under "Getting Started."
3. Select a CRAN location (a mirror site) and click the corresponding link.
4. Click on the "Download R for [your operating system]" link at the top of the page.
5. Click on the file containing the latest version of R under "Files."
6. Save the .pkg file, double-click it to open, and follow the installation instructions.
7. Running the Analysis
8. Install required R packages from within R. install.packages(c("tidyverse", "data.table", "survival", "ggplot2")). Add other packages as specified in each analysis script

# Python (≥3.8)
1. Download and install from: https://www.python.org/downloads/
2. Ensure pip is installed for managing dependencies.
3. Installation
4. Install required Python packages: pip install -r requirements.txt
(Each folder may contain its own requirements.txt file depending on the analysis.)

# This repository contains code and data used for metagenomic (MG), metatranscriptomic (MT), host transcriptomic (RNA), mouse model, and multi-omics analyses related to lung cancer and immunotherapy outcomes.
# Repository Structure
The repository is organized into five main folders:
- MG/ – Metagenomics analyses (R scripts + datasets)
- MT/ – Metatranscriptomics analyses (R scripts + datasets)
- RNA/ – Host transcriptomics analyses (R scripts + datasets)
- MOUSE/ – Mouse model analyses (R + datasets)
- MULTIOMICS/ – Integrated multi-omics analyses (Python, requires data from MG, MT, and RNA folders)

Each folder includes:
R or Python code for the analysis.
The dataset(s) needed to reproduce results.
Instructions for file/folder adjustment if you download datasets separately.

# Usage
Clone this repository:
git clone https://github.com/YourUsername/YourRepoName.git
cd YourRepoName
Navigate to the folder of interest (e.g., MG/, MT/, RNA/, MOUSE/, or MULTIOMICS/).
Adjust file paths in the scripts according to where your data files are stored.
By default, datasets should be placed in the same folder as the corresponding script.
If you store them elsewhere, update the script paths accordingly.
Run the analysis. 

# Multi-Omics Analysis
The MULTIOMICS/ folder requires input from MG, MT, and RNA analyses. Ensure you have run each of these beforehand and that the output files are correctly placed in the paths specified in the multi-omics scripts.
Citation
If you use this repository, please cite our work
