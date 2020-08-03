##============================================================================##
## Variables used
##============================================================================##

Sys.setenv(TZ='Etc/GMT-1')
dir.raw <- "D:/PhD/HOAL/raw_data/"
dir.processed <- "D:/PhD/HOAL/processed_data/"
dir.baro <- "D:/PhD/HOAL/data_lovrenc/baro/"
dir.GW <- "D:/PhD/HOAL/data_lovrenc/Groundwater/"
dir.GIS <- "D:/PhD/HOAL/GIS/"

#==============================================================================#
## HOAL zones ####
HOALzone <- PiezoInstall()[,c("Station","Group")]
HOALzone$Group <- as.factor(HOALzone$Group)
