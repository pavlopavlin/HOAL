# HOAL package
HOAL package is a collection of function mainly used for time-series manupulation and analysis.
This package was created in the course of my PhD in the Hydrological Open Air Laboratory (HOAL). As such varies in synthax and coding style.

## Install
To install the package from GitHub use:

```{r}
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("pavlopavlin/HOAL")

```

## Funcionalities

This package could be broken down to the following categories of functions:

- time-series (xts objects) manipulation 
- processing of raw groundwater data from vanEssen pressure transducers
- soil moisture data filtering and aggregation
- plotting wrappers based on *dygraph* and *ggplot2* packages
- other


