### TMMS-GUI is an application for intergrating and processing transcriptomics and metabolomics data
There are so far only a handfull of open source user friendly tools available for intergrating transcriptomics and metabolimic datasets. The aim of TMMS-GUI is to help bridge this gap.
The app itself is relatively self explainatory; however, there are a few rules that you need to follow.

1. Please make sure the sample naming convention is consistent across the transcriptomics and metabolomics datasets
2. Please make sure you have a unique identifier for your gene name and metabolomics compound name (You can use gene name and compound name directly, providing they are unique)
3. For the transcriptomics dataset, the metabolomics dataset and sample information, please ensure sample names are in rows and the genes, compound ID and sample observations are in columns. Example files (flolai/TMMS-GUI/Example_Data_shiny_app_testing) can be downloaded together with the application and some samples of the data are shown below.
4. Data must be in .csv format
5. User data can be either pre-normalised before data input or to the log2 normalisation method within the application can be used.
6. Once gene expression and metabolomics (including gene names and metabolite names) are inputted into the Data Input tab, please continue in the Plots tab for the input of sample information (treated/control etc.)
7. Choose the relevant plot selection according to the normalisation method and model you need. (Here, the data will be merged and then centred and autoscaled for modelling)


**Please be patient, with large datasets (over 10000 genes and metabolites) the application may take a few minutes to generate the plot. While the interactive plot is being generated, a basic plot will be displayed in your rstudio plot section.**

### Example transcriptomics data showing samples names in rows and genes in columns (same for metabolomics)
```sh
sample_name	gene_Y880_RS00170	gene_Y880_RS00950	gene_Y880_RS00960	gene_Y880_RS01375
P6-2-24-P	       6	              301	             19	                 130
P6-3-24-P	       9	              240	             16	                 137
P6-2-1-P	       9                  191              	 26	                 239
P6-3-1-P	       11                 90	             9                 	 145
```


### Example transcriptomics gene identifier file (same for metabolomics)

```sh
gene_ID	            gene_name
gene_Y880_RS13925	zwf
gene_Y880_RS01680	ygbF
gene_Y880_RS22710	yciH
gene_Y880_RS28945	yccS
gene_Y880_RS01650	ybgC
```

### The statistics
There are two models that can be computed using this app. The first model is Principal Component Analysis (PCA) and the second one is Orthogonal Projections to Latent Structures Discriminant Analysis (OPLS-DA). Both of these model uses the ropls package in R [1]. A very good tutorial for using this package can be found here - https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html (use of command line only)

## PCA

<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/plot_area_1.png?raw=true" width="1000">


A good explaination of PCA and the scores plot and loadings plot can be found in the links below

- https://blog.umetrics.com/what-is-principal-component-analysis-pca-and-how-it-is-used
- https://blog.bioturing.com/2018/06/18/how-to-read-pca-biplots-and-scree-plots

In a nutshell, when a group of samples are clustered together in a scores plot, it suggests that there are some similarities between them. In the example data above, on the left you have sample type P6 at different time points and different treatments so it makes sense for them to be clustered together.  The positioning of genes and metabolites in the loadings plot correlates to the positioning of samples in the scores plot. As there are likely to be many variables in a multiomic dataset, the name/ID in the loadings plot will not be shown unless you hover over to that specific gene/metabolite. You can use the group selection settings (highlighted in red) available to select a section from the loadings plot and these can be copied to a text file or csv file for further interrogation if needed.

## OPLS-DA

<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/oplsda_plot_explain2.png?raw=true" width="1000">


A very simple explaination of OPLS-DA can be found in the link below

- https://blog.umetrics.com/explaining-differences-or-grouping-data-opls-da-vs-pca-data-analysis

OPLS-DA is a discrimination method where the model will try to find the maximun separation between the two groups. In the above example, there is a good separation between the two sample types as shown in the OPLS-DA scores plot. However, like any supervised model, it is very important to validate the model. Here, the permutation plot is shown together with the model's scores plot. The purpose of the permutation test is to determine whether the separation seen in the scores plot happened by chance. The dots on the left of the permutation plots are the scrambled observation (default set to 100), dots at x-axis position 1 are the original model data. The scambled data should never be better than the original points; also, pR2Y and pQ2 should be less than 0.05 or the model is invalid. For further model validation explaination please see reference [2].

Reference:

[1]Thevenot EA, Roux A, Xu Y, Ezan E, Junot C (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. Journal of Proteome Research, 14, 3322-3335. http://pubs.acs.org/doi/full/10.1021/acs.jproteome.5b00354.

[2]Szymańska, Ewa & Saccenti, Edoardo & Smilde, Age & Westerhuis, Johan. (2012). Double-check: Validation of diagnostic statistics for PLS-DA models in metabolomics studies. Metabolomics : Official journal of the Metabolomic Society. 8. 3-16. 10.1007/s11306-011-0330-3. 


