installFromPeridotRepo <- function(packName, repo){
  install.packages(packName, dependencies=TRUE, repos=repo)
}

packIsInstalled <- function(s){
  return(s %in% rownames(installed.packages()))
}

args = commandArgs(trailingOnly = F)
packageToInstall <- args[length(args)-1]
repo <- args[length(args)]

if(packIsInstalled(packageToInstall)){
  "Package already installed. Skiping installation."
}else{
  installFromPeridotRepo(packageToInstall, repo)
  if(packIsInstalled(packageToInstall)){
    "Package successfully installed."
  }else{
    "The package could not be installed successfully."
  }
}