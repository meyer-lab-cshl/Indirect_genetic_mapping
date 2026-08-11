# Mapping Indirect Genetic Associations in Interacting Cell Types
This contains information for developing a method for analyzing the associations between genetic variants and their effect in interacting cell types. This is especially relevant to help delineate the interactions between cell types in human health and disease contexts. 

## Project Background
To characterize genetic variants that are involved in changing the interaction between cell types, I will be expanding upon the conventional expression quanititbve eQTL mapping framework. A genetic variant is considered an eQTL if the change in genotype is associated with a change in gene expression. Before implementing modifications to the framework, it is important to establish the baseline power and calibration of the conventional eQTL-mapping model. Thus, the first part of the project was to determine the baseline power and calibration of the method under a range of parameters.

After validating the power and calibration of the model, I extending the framework to include a term for a test statistic to model the impact of an upstream cell type to a downstream cell typoe. 


## Project Workflow
### Conventional eQTL model benchmarking
To benchmark the results of the conventional eQTL model, I used RESHAPE to create a synthetic genetic reference dataset from a public genetic dataset. Then, I inputted the synthetic genetic reference data into splatPop to generate simulated gene expression data with known eQTLs. I inputted a range of parameters from reviewing relevant biological literature to determine the cell size, sample size, and eQTL effect sizes. Then, I used the eQTL mapping software limixQTL to map eQTLs.
<img width="1131" height="407" alt="method_flowchart" src="https://github.com/user-attachments/assets/7e844f35-fc48-4f8e-89fe-8ff5895a07e4" />

### Simulated Indirect Genetic Effects

## Files Included in Repository
