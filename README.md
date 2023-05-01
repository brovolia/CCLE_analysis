# CCLE_analysis
Here I worked with dataset from Cancer cell line encyclopedia. I subsetted the z-score matrix according to different genes related to AhR and then performed the trajectory analysis with dynverse package. The overall analysis is memory-consuming, that's why I made it using server. 
So in this analysis this "preparation for server" step was the most important for me. 
I wrote the script in R-studio, then I attached to cluster with Putty, created a repository for a project in my home directory and transfered script and input files from local machine to server via FileZilla.
I should also created in my home directory a folder for R packages \
`mkdir R_packages`\
And then in R script I specified that I want install the packages library in this directory\
`options(repos = "https://cran.rstudio.com/")
install.packages('devtools', lib = '/home/o515a/R_packages')
library(devtools)
devtools::install_github('dynverse/dyno', lib = '/home/o515a/R_packages')`\
