# ARGCompare App

This is an interactive dashboard for exploring and comparing antibiotic resistance gene (ARG) detection tools.

## Running the App using Pixi
Pixi is a package used to provide a reproducible environment (For this project, it includes R + system dependecies + R packages, as specified in `pixi.toml` & `pixi.lock`).

### Prerequisites
- Install [pixi](https://pixi.sh)
- Clone the GitHub repository

### 1. Clone the argcompare repo and go to the app folder
```bash
git clone https://github.com/BigDataBiology/arg_compare.git
cd arg_compare/shiny/app
```
### 2. Install the Pixi environment
```bash
pixi install
```
### 3. Run the App
Make sure to run this command in your current working directory:
```bash
pixi run R -e 'shiny::runApp(host="0.0.0.0", port =3838)'
```
Then you click on the url link it shows to view the app.

## Running the App in  Rstudio
You can also run the app directly in Rstudio if you prefer.

### 1. Clone the repository, open a project and set the working directory.
```
git clone https://github.com/BigDataBiology/arg_compare.git
```

### 2. Install required packages
Run this command in your console.
```
install.packages(c(
  "shiny",
  "bslib",
  "tidyverse",
  "shinyWidgets",
  "dplyr",
  "ggplot2",
  "gridExtra",
  "RColorBrewer",
  "ggpattern",
  "reactable",
  "Cairo",
  "ggalluvial",
  "cowplot",
  "scales",
  "magrittr",
  "shinycssloaders",
  "ragg",
  "qs",
  "ggrastr"
))
```
### 3. Run the app
There are two ways to do this:
1. Click **Run App** (top-right) of the Rstudio UI
2. Run this command from the console:
```r
shiny::runApp('shiny/app')
```
