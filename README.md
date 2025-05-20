# Introduction
Patients with von Hippel–Lindau (VHL) disease are highly susceptible to renal cell carcinoma (RCC) due to loss of VHL function and consequent constitutive activation of hypoxia-inducible factor 2α (HIF-2α) (Eric Jonasch et al, 2021, NEJM). In 2023, the [FDA](https://www.fda.gov/drugs/resources-information-approved-drugs/fda-approves-belzutifan-advanced-renal-cell-carcinoma#:~:text=The%20most%20common%20adverse%20reactions,(formerly%20Twitter)%20@FDAOncology%20) approved belzutifan (Welireg, Merck), a selective HIF-2α inhibitor, for the treatment of advanced RCC following PD-1/PD-L1 and VEGF-TKI therapies.

# Motivation

We used single-cell RNA-seq (scRNA-seq) data (GSE269826_Cultured.rds.gz) from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269826) to investigate transcriptional changes induced by belzutifan treatment. We aim to characterize the cellular and molecular consequences of HIF-2α inhibition and gain deeper insight into VHL-associated biological processes in renal cell carcinoma.

# Repository Structure
-`scRNA_data_processing_1.R` : Data pre-processing, normalization, dimensionality reduction, clustering, and differential expression after pseudo-time trajectory 

-`scRNA_data_TF_GRN_inference.py` : Regulon activity inference after pruning, visualize regulon activity along pseudo-time

-`belzutifan_RCC_V1.pdf` : Visualized result data
