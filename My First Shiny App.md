### Simple app for intergrating transcriptomics and metabolomics data
There are so far only a handfull of open source user friendly tools available for intergrating transcriptomics and metabolimic dataset. It is hope that the easy to use TMMS-GUI will help this growing field.
The app itself is relatively self explainatory; however, there are some simple rules that you so need to follow.

1. Please make sure sample naming convention is consistance across transcriptomics and metabolomics data set
2. Please make sure you have a unique identifier for your gene name and metabolomics name (You can use gene name and compound name directly, providing they are unique)
3. For both transcriptomics and metabolomics data set, please ensure sample names are in rows, genes and compound ID as columns. Example files can be downloaded together with the APP and example snapshots are shown below.
4. data must be in .csv format
4. User data can be either pre-normalised before data input or to use the log2 normalisation method within the application.

**Please be patient, large data set (over 10000 genes and metabolites) may take a few minutes to generate the plot. If you are in a rush, basic model graphic will be display in your rstudio plot section.**

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
gene_ID	            gene
gene_Y880_RS13925	zwf
gene_Y880_RS01680	ygbF
gene_Y880_RS22710	yciH
gene_Y880_RS28945	yccS
gene_Y880_RS01650	ybgC
```

### The statistics
There are two models that can be computed using this app. The first model is Principal component analysis (PCA) and Orthogonal Projections to Latent Structures Discriminant Analysis (OPLS-DA). Both of these model uses the ropls package in R[1]. A very good tutorial for using this package can be found here - https://www.bioconductor.org/packages/devel/bioc/vignettes/ropls/inst/doc/ropls-vignette.html (use of command line only)

## PCA

<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/plot_area_1.png?raw=true" width="1000">


A good explaination of PCA scores plot and loadings plot explaination can be found in the website below

- https://blog.umetrics.com/what-is-principal-component-analysis-pca-and-how-it-is-used
- https://blog.bioturing.com/2018/06/18/how-to-read-pca-biplots-and-scree-plots

Briefly, when a group of samples clustered together in a scores plot, it is likely that there are some similarities between them. In the example data above, on the left you have sample type P6 at different time points and different treatment so it make sense for them to be clustered together.  The positioning of genes and metabolites in the loadings plot correlates to the positioning of samples in the scores plot. As there are likely to be many variables in a multiomic dataset, the name/ID in the loadings plot will not be shown unless you hover over to that specific gene/metabolite. You can use the group selection settings (highlighted in red) available to select a section from the loadings plot and these can be copied to a text file or csv file for further interrogation if needed.

## OPLS-DA

<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/oplsda_plot_GUI.png?raw=true" width="1000">


A very simple explaination of OPLS-Da model can be found in the link below
- https://blog.umetrics.com/explaining-differences-or-grouping-data-opls-da-vs-pca-data-analysis

OPLS-DA is a discrimination method where the model will try to find the maximun separation between the two groups. In the above example, there is a good separation between the two sample tpye is shown in the OPLS-DA scores plot. However, like any supervised models, it is very important to validate the model. In here, the permutation plot are shown together with the model scores plot. The purpose of the permutation test is to determine whether the separation seen in the scores plot can be happened by chance. The 'blue' dots in the perputation plots are the scrambled observation labels (default set to 100) amd the 'red' dots and the 'blue' dots at 1 in the x-axis are the original points. The scambled data should never be better than the original points or the model is invalid. For further model validation explaination please find in reference [2].

Reference:

[1]Thevenot EA, Roux A, Xu Y, Ezan E, Junot C (2015). Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. Journal of Proteome Research, 14, 3322-3335. http://pubs.acs.org/doi/full/10.1021/acs.jproteome.5b00354.

[2]Szymańska, Ewa & Saccenti, Edoardo & Smilde, Age & Westerhuis, Johan. (2012). Double-check: Validation of diagnostic statistics for PLS-DA models in metabolomics studies. Metabolomics : Official journal of the Metabolomic Society. 8. 3-16. 10.1007/s11306-011-0330-3. 


