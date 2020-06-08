### Simple app for intergrating transcriptomics and metabolomics data
There are so far only a handfull of open source user friendly tools available for intergrating transcriptomics and metabolimic data set. It is hope that the easy to use TMMS-GUI will help this growing field.
The app itself is relatively self explainatory, but there are some simple rule (mainly common sence) that you so need to follow.

1. Please make sure sample naming convention is consistance across transcriptomics and metabolomics data set
2. Please make sure you have a unique identifier for your gene name and metabolomics name (You can use gene name and compound name directly, providing they are unique)
3. For both transcriptomics and metabolomics data set, please ensure sample names are in rows, genes and compound ID as columns. Example files can be downloaded together with the APP and example snapshots are listed below.


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

Add csv file
add sample details
brief explaination of what to expect

Need to add titles to plot and instruction in usage


