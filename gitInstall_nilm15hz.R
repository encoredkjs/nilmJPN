#! /usr/bin/Rscript
options(unzip = 'internal')

install.packages('devtools', dependencies = TRUE)
library(devtools)

### Set basic information for the installation of private repository
setPAT = '640e2849d4b51ab2de08486e28dbc2985008ee8b'

setPackage = 'encoredDataETL'
install_github('EncoredTech/research_R_package', subdir = setPackage, auth_token = setPAT)

setPackage = 'encoredLog'
install_github('EncoredTech/research_R_package', subdir = setPackage, auth_token = setPAT)

setPackage = 'nilm15hz'
install_github('EncoredTech/research_R_package', subdir = setPackage, auth_token = setPAT)
