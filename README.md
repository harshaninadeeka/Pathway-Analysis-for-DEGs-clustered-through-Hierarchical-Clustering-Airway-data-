# Pathway-Analysis-for-DEGs-clustered-through-Hierarchical-Clustering-Airway-data-
R-based bioinformatics pipeline reproducing key findings from the Himes et al. (2014) airway dataset. Features hierarchical clustering of DEGs in dexamethasone-treated human Airway Smooth Muscle (hASM) cells, followed by Reactome-based Over-Representation Analysis (ORA) to map pathways regulating airway remodeling and cell adhesion
This repository contains a complete, end-to-end transcriptomics workflow that analyzes the effects of Dexamethasone (a potent synthetic glucocorticoid) on Human Airway Smooth Muscle (hASM) cells. 
By leveraging the robust RNA-Seq dataset originally published by Himes et al. (2014), this project demonstrates how computational pipelines can dissect the molecular mechanisms behind steroid-mediated therapies for asthma and chronic inflammatory airway diseases. 

### 🔬 The Science Behind the Project 
When hASM cells are treated with 1µM dexamethasone for 18 hours, a complex cascade of gene expression changes occurs. To make sense of these differentially expressed genes (DEGs), this pipeline performs hierarchical clustering to group genes by expression pattern, followed by Over-Representation Analysis (ORA)using the ReactomePA package in R. 
This analysis successfully maps these gene clusters to three major biological themes that mirror the core findings of the original paper: 1. 

**Cluster 1 (Neuronal & Cell Adhesion): Highlighting NCAM1-mediated adhesion pathways that modulate neuromuscular crosstalk to reduce airway hyperresponsiveness. 
**Cluster 2 (Metal Homeostasis & Stress Response): Upregulation of metallothionein-driven antioxidant responses that protect hASM cells against oxidative stress. 
**Cluster 3 (Extracellular Matrix Remodeling): Pinpointing the structural shifts, ECM degradation, and O-linked glycosylation pathways that dexamethasone targets to arrest pathological airway scarring and tissue remodeling. 

###Visualization



##cluster 1-Dotplot and Bar plot


<img width="631" height="507" alt="Dotplot ORA-Cluster 01" src="https://github.com/user-attachments/assets/318f3a55-b085-46db-b4d0-ef2a6a1c5c2d" />

<img width="631" height="507" alt="Barplot ORA-Cluster 1" src="https://github.com/user-attachments/assets/a40631b0-2551-4484-9854-0ae9bd17be96" />


##cluster 2-Dotplot and Bar plot

<img width="631" height="507" alt="Dotplot ORA -Cluster 2" src="https://github.com/user-attachments/assets/a97c5411-c18c-42f7-af98-1c0715d8ff68" />


<img width="631" height="507" alt="Barplot-Cluster 02" src="https://github.com/user-attachments/assets/12d643f0-0ed3-405c-a511-604a7c47966c" />


##cluster 3-Dotplot and Bar plot

<img width="631" height="507" alt="Dotplot-Cluster 03" src="https://github.com/user-attachments/assets/13deea88-ee08-4bad-a439-124604afb709" />


<img width="631" height="507" alt="Barplot-Cluster 03" src="https://github.com/user-attachments/assets/f1a552de-e80b-4f82-ad3c-e743d417ea63" />




## 🛠️ Bioinformatics Pipeline & Methodology 
* **Differential Expression Analysis: Utilizing standard Bioconductor methods to isolate high-confidence DEGs. Identified significant DEGs p. adj < 0.05 from the raw count matrix.
* **Hierarchical Clustering: Grouping DEGs by normalized z-scores to isolate distinct therapeutic response profiles. 
* ** Pathway Enrichment: Performing ORA against the Reactome database to identify biologically relevant pathways. (0.13 < p. adj < 0.22).


Leveraging my molecular biotechnology background, this project helped me to improve my skills in R-based transcriptomics, from identifying DEGs from RNA seq data, hierarchical clustering of differentially expressed genes to Reactome pathway enrichment analysis (ORA) of glucocorticoid-treated airway cells.

I will upload my R script, cluster files and visualized plots.
