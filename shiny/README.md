# ARG Detection Tools - Shiny app

An interactive Shiny dashboard for exploring and comparing antibiotic resistance gene (ARG) detection tools.


To run the app:
- install R
- Do in terminal: 
```bash
R -e "install.packages(c('shiny','bslib','tidyverse','shinyWidgets','Cairo','ggalluvial'), repos='https://cloud.r-project.org')"
```
- Then:
```bash
R -e "shiny::runApp()"
```