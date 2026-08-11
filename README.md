# Mapping Indirect Genetic Associations in Interacting Cell Types
Within the immune systems, cell are constantly interact to work together and regulate each other to help protect the body from invaders. Thus, the effects of a genetic variation that act in one cell type can have downstream impacts on another cell type. This project builds a method to detect these "indirect" genetic effects to help delineate important pathways in human health and disease-related contexts. 

## Project Background
To characterize genetic variants that are involved in changing the interaction between cell types, I will be expanding upon the conventional expression quanititbve eQTL mapping framework. A genetic variant is considered an eQTL if the change in genotype is associated with a change in gene expression. Before implementing modifications to the framework, it is important to establish the baseline power and calibration of the conventional eQTL-mapping model. Thus, the first part of the project was to determine the baseline power and calibration of the method under a range of parameters.

After validating the power and calibration of the model, I extended the framework to include a term for a test statistic to model the impact of an upstream cell type to a downstream cell type.


## Project Workflow
### Conventional eQTL Model Benchmarking
To benchmark the results of the conventional eQTL model, I used RESHAPE to create a synthetic genetic reference dataset from a public genetic dataset. Then, I inputted the synthetic genetic reference data into splatPop to generate simulated gene expression data with known eQTLs. I inputted a range of parameters from reviewing relevant biological literature to determine the cell size, sample size, and eQTL effect sizes. Then, I used the eQTL mapping software limixQTL to map eQTLs.
<img width="1131" height="407" alt="method_flowchart" src="https://github.com/user-attachments/assets/7e844f35-fc48-4f8e-89fe-8ff5895a07e4" />

### Indirect Genetic Associations
After benchmarking the conventional eQTL model, I extended it to detect indirect genetic
associations due to downstream impacts of changes in one cell type due to genetic variation. Similar to benchmarking the conventional eQTL, I simulated genetic data and gene expression data with RESHAPE and splatPop respectively.

I simulated this in two ways, based on two hypotheses of how these effects might arise:
- **Overall model**: many small genetic effects in one cell type add up to shift gene expression
  in another cell type.
- **Pathway model**: a single variant with a strong effect is enough to cause downstream changes
  on its own.
To detect these effects, I summarize the genetic activity in the "upstream" cell type into a
single molecular state score per sample, then add that score into the eQTL model for the
"downstream" cell type. I tested this approach on simulated data with a known ground truth to
see how well it could recover the indirect effects under each hypothesis.

Both of these projects are implemented as Snakemake pipelines to ensure reproducibility, quality, and easily extensible. 

## Files Included in Repository
