# Citizen Science Deepens the Ecological and Climatic Dimensions of Mosquito Surveillance

**ABSTRACT**

As mosquito-borne diseases continue to expand worldwide, integrating citizen science into vector surveillance presents untapped potential. This study compares ecological models of \textit{Aedes albopictus}, an invasive mosquito and global vector of dengue and other arboviruses, in Spain (2020–2022), using two contrasting data sources: traditional traps and citizen science. While both showed strong seasonal agreement, spatial discrepancies emerged in summer, particularly in northwestern Spain. These differences were linked to the broader environmental coverage provided by citizen observations, especially under extreme temperature and humidity conditions. Citizens can report mosquito activity in places and times where traps are absent or ineffective, capturing real host-seeking behavior and reducing sampling biases inherent to fixed trap locations. Incorporating citizen data into trap-based models improved monthly predictive performance, demonstrating its value as a complementary source. Combined, they enhance mosquito surveillance, support early warning systems, and inform decisions under climate-driven range expansion.
## Authors and Afiliations
Catuxa Cerecedo-Iglesias—Centre d'Estudis Avançats de Blanes (CEAB-CSIC), Blanes, Spain (catuxa.cerecedo@ceab.csic.es).<br>
John R.B. Palmer—Universitat Pompeu Fabra (UPF), Barcelona, Spain (john.palmer@upf.edu).<br>
Frederic Bartumeus—Centre d'Estudis Avançats de Blanes (CEAB-CSIC), Blanes, Spain; Centre de Recerca Ecològica i Aplicacions Forestals (CREAF), Cerdanyola del Vallès, Barcelona, Spain; Institut Català de Recerca i Estudis Avançats (ICREA), Barcelona, Spain (fbartu@ceab.csic.es).

*Corresponding author: Catuxa Cerecedo-Iglesias

**CONTENTS**

This repository contains the scripts for the primary statistical analyses presented in the manuscript titled *Citizen Science Deepens the Ecological and Climatic Dimensions of Mosquito Surveillance* submitted on *Nature Ecology & Evolution*. In some cases, it also provides details on how the figures included in the manuscript were generated.

The `DATA` directory contains cleaned datasets derived from both traditional sampling and citizen science observations, along with the covariates used to build the models. It also includes additional datasets related to land use and other spatial variables required for georeferencing the analyses and figures.

Scripts are organized sequentially (from `R1` to `R6`) to reflect the workflow followed during the development of this study.

**PREREQUISITES**
- R version(4.4.2.)
- Requieres R packages: `brms`, `tidyverse`, `rstanarm`, `loo`, `cmdstanr`, `tidyverse`, `sf`, `parallel`, `dplyr`, `data.table`, `terra`, `patchwork`, `ggplot2`, `cowplot`, `glmmTMB`, `DHARMa`, `MuMIn`, `caret`.




