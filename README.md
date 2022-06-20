BDNF DNA Methylation Map
================
Lacey W. Heinsberg

# Copyright information

Copyright 2022, University of Pittsburgh. All Rights Reserved. License:
[GPL-2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# Repository description

This repository contains R code used to create the BDNF DNA Methylation
(DNAm) Map, a web-based shiny application developed to synthesize the
results of a systematic review paper:

Treble-Barna A, Heinsberg LW, Breazeale S, Davis TS, Stec Z, Kesbhat AA,
Chattopadhyay A, VonVille HM, Ketchum AM, Yeates KO, Kochanek PM, Weeks
DE, Conley YP. Brain-Derived Neurotrophic Factor (BDNF) Epigenomic
Modifications and Brain-Related Phenotypes in Humans: A Systematic
Review. (LINK TO BE ADDED UPON PUBLICATION).

The purpose of this paper was to systematically review studies that
investigated BDNF epigenomic modifications in association with
brain-related phenotypes in humans and to present the findings in a
researcher-friendly format.

The BDNF DNAm Map was developed to support users in interactively
visualizing the specific positions of BDNF DNAm CpG sites investigated
across all studies for which these data were available. The application
displays a to-scale depiction of the BDNF gene and allows users to
interactively customize the figure to display DNAm sites customizable by
phenotype(s), tissue type(s), statistical significance, reference
identifier (RefID), and CpG position number.

The BDNF DNAm Map is available at
<https://lwheinsberg.shinyapps.io/BDNF_DNAmMap/> for public use.

The purpose of this repository is to more fully document the application
construction. The code is presented as a single file, `app.R`, which
calls two data files, `BDNF_BaseCoordinates.xlsx` and
`BDNF_DNAmMap_AKA_Long.xlsx`, located within this repository.

# BDNF DNAm Map development

The BDNF DNA Methylation Map was developed in R version 4.1.2 using a
suite of packages including shiny (Chang et al., 2015), ggplot2
(Wickham, 2016), and plotly (Sievert, 2020). The application runs on
data synthesized through the systematic review described above, curated
in long form.

In brief, we created a to-scale depiction of the BDNF gene (base plot),
including exons, promoter regions, and transcription start sites. Next,
CpG sites and results for a single study x CpG x phenotype x tissue were
dynamically layered on the figure and then presented within the shiny
application with customizable filters. The application uses reactive
functions designed to be responsive to user-defined input as well as
interactive functions that produce hover text and allow the user to zoom
in (user-defined plot). All figure data were mapped to human genome
assembly 38 (hg38).

For more complete details on the application configuration and
development, please refer to the full source code (`app.R`) located in
this repository.

# Dependencies

The R code, `app.R`, relies on calling `BDNF_BaseCoordinates.xlsx`,
`BDNF_DNAmMap_AKA_Long.xlsx`, and the R package dependencies listed
within the “Load libraries” sections of application source code.

# Execution

R Studio can be used to execute the .Rmd code using the `Run App`
function from the quick bar. Alternatively, users can interact with an
[online version](https://lwheinsberg.shinyapps.io/BDNF_DNAmMap/) that
requires no R software or code.

# Instructions for use

Instructions for use are located within the application.

# Developers

This application was created by Lacey Heinsberg with oversight from
Daniel Weeks. Note this application could not function without the
accompanying database curated by Amery Treble-Barna and her team which
included Lacey Heinsberg, Tara Davis, Stephen Breazeale, Zachary Stec,
Aboli Kesbhat, Julia Tefs, and Ansuman Chattopadhyay. This curation took
a great deal of collaborative effort and many many many hours of work.

# Contact information

If you have any questions or comments, please contact Lacey W.
Heinsberg, PhD, RN (<law145@pitt.edu>).

# References

Chang W, Cheng J, Allaire J, et al. shiny: Web Application Framework for
R. 2021. <https://cran.r-project.org/package=shiny>.

Wickham H. Ggplot2: Elegant Graphics for Data Analysis. New York:
Springer-Verlag; 2016. <https://ggplot2.tidyverse.org>.

Sievert, C., 2020. Interactive Web-Based Data Visualization with R,
plotly, and shiny. Chapman and Hall/CRC Florida.
