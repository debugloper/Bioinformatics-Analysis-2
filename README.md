# Bioinformatics-Analysis-2
Typical Microarray statistical analysis and comparative studies

The title of the study is: "Multivariate analysis of microarray data analysis".

This script was produced in an attempt to study about microarray gene expression of Taiwanese Lung Cancer Patients.
The data was retrieved from GEO database which is searchable at: http://www.ncbi.nlm.nih.gov/sites/GDSbrowser/. 
GDS3837 data is consisting of 60 normal patient samples and 60 cancer patient samples, with total of 54,675 gene expression by
120 columns of samples. 

The project aims to study about microarray and statistical analysis method by using microarray data which is numerical value
of gene expression. It is an analysis about the GDS3837 data, which is paired tumor and normal tissues for gene expression
microarray data of non-smoking female lung cancer in Taiwan. GDS3837 data is provided by gene expression omnibus site GEO.
Before analysis, microarray data was normalized in order to remove bias according to factor except biological factor. The script
performs paired t-test for confirming difference between normal and tumor tissue using the normalization data. We are able to find
significant genes from the p-value of paired t-test and significant genes are used to do multivariate analysis such as 
principal component analysis (PCA), independent component analysis, random forest clustering and hierarchical clustering. We also 
find differentially expressed genes from the whole genes by implementing the multiple comparison method.
