# Sooyoun Lee (Saxophone) BF528 Individual Project - Microarray Based Tumor Classification 

## Contributor: Sooyoun Lee leesu@bu.edu

The purpose of this project is to compare and identify the C3 and C4 tumor subtypes by using the 134 tumor samples by performing noise filtering, dimensionality reduction, hierarchical clustering, subtype clustering, and gene set analysis. 

The following paper was used as a reference:  Marisa et al. Gene Expression Classification of Colon Cancer into Molecular Subtypes: Characterization, Validation, and Prognostic Value. PLoS Medicine, May 2013. PMID: 23700391


In this project our goal is the following:
- Employ noise filtering techniques to reduce data dimensionality.
- Perform data-driven analyses, such as hierarchical clustering, to discover novel relationships among samples in a given dataset.

# Repository Contents
#### 1. data.table
The data.table library called since this library helps with subsetting the rows, and compute on columns and also perform aggregations by the groups. 

#### 2. hclust
The hclust (Hierarchical clustering) function was used to analyze the set of dissimilarities and methods. 

#### 3. cutree
A cutree function was used to cut hierarchical models. This function allows us to cut the tree based on a certain height h or by a certain number of clusters k. 

#### 4. p.adjust
A p.adjust function was used to adjust p-values from a set of an un-adjusted p-values, using a number of adjustment procedures. 

#### 5. BiocManager 
The BiocManager package was installed since this function focuses on the statistical analysis and comprehension of high-throughput genomic data. 

#### 6. GSEABase
The GSEABase function is a universal gene set enrichment analysis tools that was installed in this project. 
