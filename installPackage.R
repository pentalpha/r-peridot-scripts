installFromBioLite <- function(packName){
  source("https://bioconductor.org/biocLite.R")
  biocLite()
  biocLite(packName)
}

installFromCran <- function(packName){
  install.packages(packName, dependencies=TRUE, repos="https://cloud.r-project.org/")
}

installFromPeridotRepo <- function(packName, repo){
  if(repo != "0.0.0.0"){
    install.packages(packName, dependencies=TRUE, repos=repo)
  }
}

packIsInstalled <- function(s){
  return(s %in% rownames(installed.packages()))
}

args = commandArgs(trailingOnly = F)
packageToInstall <- args[length(args)-1]
repo <- args[length(args)]

installFromPeridotRepo(packageToInstall, repo)

if(packIsInstalled(packageToInstall)){
  "Package successfully installed."
}else{
  "Could not install from R-Peridot repositories, trying to download from a CRAN mirror."
  installFromCran(packageToInstall)
  if(packIsInstalled(packageToInstall)){
    "Package successfully installed."
  }else{
    "Could not install from CRAN repository, trying to download from Bioconductor."
    installFromBioLite(packageToInstall)
    if(packIsInstalled(packageToInstall)){
      "Package successfully installed."
    }else{
      "The package could not be installed successfully."
    }
  }
}