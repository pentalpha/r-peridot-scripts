R.Version()

Packages = installed.packages()

Packages = as.data.frame(Packages, row.names = FALSE)

onlyNamesAndVersion = Packages[,c(1,3)]

onlyNamesAndVersion
