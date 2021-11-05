# scRNA-zebrafish

Code to perform the analysis presented in the manuscript "Single cell transcriptomics reveals critical molecular programming of progenitors for motor neurons and oligodendrocytes in zebrafish".


Required software:
Python3 ,R (3.6.0),

Description
1.main.R: Filter out low-quality cells and non-expressed genes. Batch effect correct. Seperated cells into differential clusters based on genes' expression. Calculate G2M score; plot cell ratio between time groups. Plot cells in UMAP dimension; Find markers and plot;Functional enrichment of marker genes or DEGs; 
2. paga.all.py : PAGA analysis for all the cell types found in this research
3. lineage.R : lineage analysis by R monocle package
4. find.priOPC.r: script for detect pri opc cell type;
5. Metabolic.pmn.R: data analysis about metabolic related genes expressed in PMN sub celltype;
6. Lineage.pmn.R:lineage analysis of PMN sub celltype by R monocle package
