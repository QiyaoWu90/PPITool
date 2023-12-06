# About

PPI Tool is for protein-protein interactor screening from mass spectrometry data.

It can be applied to both affinity-purification mass spectrometry data and proximity-labeling mass spectrometry data. PSM information is required but only label-free DDA data has been tested.

PPI tool started by modifying Fisher's exact test. We suppose the distribution of true and false interactors from PPI studies follows the same hypergeometric distribution. Generally, data with biological replicates are recommended but data without replicates are also suitable.

# Usage

Rstudio is recommended to launch the software

1. request the shiny package first:

``` bash
if(!require(shiny)){
  install.packages("shiny")
  library(shiny)
}
```

2. Run the following command to launch software:

``` bash
runGitHub( "PPITool", "QiyaoWu90")
```



# Input examples

Two examples can be found in input_example folder (three/four repeats dataset)

# How to use

Please check the "Documentation" Tab of the software for the usage.
