
# TMMS-GUI Shiny
### Transcriptomics and metabolomics data multivariate statistical analysis
- This application is a concatenation-based multivariate analytical tool for integrating transcriptomic and metabolomic dataset.
- The analytical algorithm used in this application is based on the ropls package[1]
- The example data used during the development of this appilcation based reference [2] and the subset of the data can be found in the folder 'Example_Data_shiny_app_testing' 
- More detailed user guide can be found within the application
###### Please use the examples files in /TMMS-GUI/Example_Data_shiny_app_testing/ as a reference for file format

## Installation instruction:

```sh
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ropls")
install.packages("shiny")
shiny::runGitHub("TMMS-GUI", "flolai")
    
```    

### PCA in TMMS-GUI
<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/plot_area_1.png?raw=true" width="600">

### OPLS-DA in TMMS-GUI
<img src="https://github.com/flolai/TMMS-GUI/blob/master/app_graphics/oplsda_plot_GUI.png?raw=true" width="600">

### Contact
Florence : clai04@mail.bbk.ac.uk

Reference:

[1]Analysis of the human adult urinary metabolome variations with age, body mass index and gender by implementing a comprehensive workflow for univariate and OPLS statistical analyses. Thevenot EA, Roux A, Xu Y, Ezan E, Junot C (2015). Journal of Proteome Research, 14, 3322-3335. http://pubs.acs.org/doi/full/10.1021/acs.jproteome.5b00354.

[2]Comparative Metabolomics and Transcriptomics Reveal Multiple Pathways Associated with Polymyxin Killing in Pseudomonas aeruginosa
Mei-Ling Han, Yan Zhu, Darren J. Creek, Yu-Wei Lin, Alina D. Gutu, Paul Hertzog, Tony Purcell, Hsin-Hui Shen, Samuel M. Moskowitz, Tony Velkov, Jian Li
mSystems Jan 2019, 4 (1) e00149-18; DOI: 10.1128/mSystems.00149-18
