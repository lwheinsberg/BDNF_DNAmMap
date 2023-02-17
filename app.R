# BDNF DNA Methylation Brain Map: https://lwheinsberg.shinyapps.io/BDNF_DNAmMap/
# Application created by Lacey W. Heinsberg, PhD, RN with oversight from Daniel E. Weeks, PhD 
# Copyright 2022, University of Pittsburgh. All Rights Reserved.
# License: GPL-2 https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html

# Load libraries
library(openxlsx)
library(ggplot2)
library(tidyr)
library(plotly)
library(dplyr)
library(stringr) 
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(dplyr)
library(ggstance) # Vertical dodge command
library(viridis)
library(htmlwidgets) # Add clickable hyperlinks
library(rsconnect)
library(shinycssloaders) # Add progress bar for plots
library(DT)

############ PART 1: DATA LOADING AND PREP
# Load and clean database files  
# Read in study database 
df <- read.xlsx("BDNF_DNAmMap_AKA_Long.xlsx") # Long form 
# Trim white space
df$Position <- str_trim(df$Position, side = c("both"))
df$Stop <- str_trim(df$Stop, side = c("both"))
df$Tissue <- str_trim(df$Tissue, side = c("both")) 
df$Study_Probe_ID <- str_trim(df$Study_Probe_ID, side = c("both"))
# Ensure numbers are numeric 
df$Position <- as.numeric(as.character(df$Position))
df$Stop <- as.numeric(as.character(df$Stop))
# Create an Illumina-specific probe column (as this information may be of specific interest to some users)
df$Illumina <- df$Study_Probe_ID
df$Illumina[!startsWith(df$Illumina, "cg")] <- NA # Set anything that doesn't start with "cg" (Illumina annotation) to missing 

# Read in base gene information (data base created using Pruunsild, P., Kazantseval, A., Aid, T., Palm, K., Timmusk, T., 2007. Dissecting the human BDNF locus: Bidirectional transcription, complex splicing, and multiple promoters. Genomics 90. https://doi.org/10.1016/j.ygeno.2007.05.004)
base <- read.xlsx("BDNF_BaseCoordinates.xlsx")

# Create CpG only database (for laying segments on initial base plot)
# First drop any multi CpG sites 
cpgs3 <- subset(df, is.na(Stop))
myvars <- c("Position", "Study_Probe_ID")
cpgs2 <- cpgs3[myvars]
# Drop any non-illuimna probe IDs from this dataframe 
cpgs2$Study_Probe_ID[!startsWith(cpgs2$Study_Probe_ID, "cg")] <- NA 
names(cpgs2)[2] <- "Illumina"
# De-duplicate 
cpgs <- unique(cpgs2)

# Prepare df (for laying CpG segments on the plots)
# Set heights of CpGs
df$x <- df$Position 
df$xend <- df$Position
df$y <- -2.5
df$yend <- -1.7

# Set y-value for initial height of first study-specific dot 
# (used along with position_dodgev() so dots don't lay on top of each other)
df$y2 <- -1.5 # Single CpGs
df$y2_multi <- -2.7 # Mean methylation or methylation ratio 

# Prepare cpgs (for laying segments on initial base plot) 
cpgs$x <- cpgs$Position
cpgs$xend <- cpgs$Position
cpgs$y <- -2.5
cpgs$yend <- -1.7

# Prepare TSS (for laying TSS segments on base plot)
temp <- which(base$name=="ts1")
tss <- base[temp:nrow(base),]
tss$short_label <- NULL
tss$label <- NULL 
tss$start <- NULL 
tss$x <- tss$stop
tss$xend <- tss$stop
tss$y <- -2.4
tss$yend <- -1.8

# Create list of phenotypes (for use in Shiny App)
phenotypes <- unique(df$Phenotype)
# Order as alphabetical 
phenotypes <- sort(phenotypes)

# Create list of tissues (for use in Shiny App)
tissues <- unique(df$Tissue)
tissues <- sort(tissues)

# Create list of refIDs (for us in Shiny App)
refIDs <- unique(df$RefID)
refIDs <- sort(refIDs)

# Create list of CpGs (for us in Shiny App)
cpgs_names <- cpgs$Position
cpgs_names <- sort(cpgs_names)

# Min and max CpG locations 
min.cpg <- min(df$Position)
max.cpg <- max(df$Position)

############ Draw figure (BDNF base plot)
# Store start/stop positions to use as coordinates in geom_rect below
gene <- which(base$name=="gene")
exon_1 <- which(base$name=="exon_I")
exon_2 <- which(base$name=="exon_II")
exon_3 <- which(base$name=="exon_III")
exon_4 <- which(base$name=="exon_IV")
exon_5 <- which(base$name=="exon_V")
exon_5h <- which(base$name=="exon_Vh")
exon_6 <- which(base$name=="exon_VI")
exon_7 <- which(base$name=="exon_VII")
exon_8 <- which(base$name=="exon_VIII")
exon_8h <- which(base$name=="exon_VIIIh")
exon_9 <- which(base$name=="exon_IX")
promoter_1 <- which(base$name=="pr1")
promoter_2 <- which(base$name=="pr2")
promoter_3 <- which(base$name=="pr3")
promoter_4 <- which(base$name=="pr4")
promoter_5 <- which(base$name=="pr5")
promoter_6 <- which(base$name=="pr6")

# Set limits of plot
my.theme <- list(theme(panel.border = element_rect(colour = "black", fill=NA, size=1)))
p <- ggplot() + ylim(-4, 2) + xlim(base$start[gene]-3509, base$start[gene]+29991) + theme_classic()  + my.theme +  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())

# Add gene body 
p <- p + geom_rect(data=base, aes(xmin=start[gene],
                                  xmax=stop[gene], ymin=-2.15, ymax=-2.00),
                   colour = "black", fill = "white")

# Add exons
# Exon 1
p <- p + geom_rect(data=base, aes(xmin=start[exon_1], 
                                  xmax=stop[exon_1], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
# R hack to add in hover text (since hover text won't work with geom_rect)
p <- p + geom_segment(data = base, aes(x = start[exon_1], y = -2.30, xend = stop[exon_1], yend = -2.30, text=paste("Exon I")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_1], y = -2.1, xend = stop[exon_1], yend = -2.1, text=paste("Exon I")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_1], y = -1.90, xend = stop[exon_1], yend = -1.90, text=paste("Exon I")), col="transparent")
# Exon 2 
p <- p + geom_rect(data=base, aes(xmin=start[exon_2], 
                                  xmax=stop[exon_2], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_2], y = -2.30, xend = stop[exon_2], yend = -2.30, text=paste("Exon II")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_2], y = -2.1, xend = stop[exon_2], yend = -2.1, text=paste("Exon II")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_2], y = -1.90, xend = stop[exon_2], yend = -1.90, text=paste("Exon II")), col="transparent")
# Exon 3, etc. 
p <- p + geom_rect(data=base, aes(xmin=start[exon_3], 
                                  xmax=stop[exon_3], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_3], y = -2.30, xend = stop[exon_3], yend = -2.30, text=paste("Exon III")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_3], y = -2.1, xend = stop[exon_3], yend = -2.1, text=paste("Exon III")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_3], y = -1.90, xend = stop[exon_3], yend = -1.90, text=paste("Exon III")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_4], 
                                  xmax=stop[exon_4], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_4], y = -2.30, xend = stop[exon_4], yend = -2.30, text=paste("Exon IV")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_4], y = -2.1, xend = stop[exon_4], yend = -2.1, text=paste("Exon IV")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_4], y = -1.90, xend = stop[exon_4], yend = -1.90, text=paste("Exon IV")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_5], 
                                  xmax=stop[exon_5], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_5], y = -2.30, xend = stop[exon_5], yend = -2.30, text=paste("Exon V")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_5], y = -2.1, xend = stop[exon_5], yend = -2.1, text=paste("Exon V")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_5], y = -1.90, xend = stop[exon_5], yend = -1.90, text=paste("Exon V")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_5h], 
                                  xmax=stop[exon_5h], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_5h], y = -2.30, xend = stop[exon_5h], yend = -2.30, text=paste("Exon Vh")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_5h], y = -2.1, xend = stop[exon_5h], yend = -2.1, text=paste("Exon Vh")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_5h], y = -1.90, xend = stop[exon_5h], yend = -1.90, text=paste("Exon Vh")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_6], 
                                  xmax=stop[exon_6], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_6], y = -2.30, xend = stop[exon_6], yend = -2.30, text=paste("Exon VI")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_6], y = -2.1, xend = stop[exon_6], yend = -2.1, text=paste("Exon VI")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_6], y = -1.90, xend = stop[exon_6], yend = -1.90, text=paste("Exon VI")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_7], 
                                  xmax=stop[exon_7], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_7], y = -2.30, xend = stop[exon_7], yend = -2.30, text=paste("Exon VII")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_7], y = -2.1, xend = stop[exon_7], yend = -2.1, text=paste("Exon VII")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_7], y = -1.90, xend = stop[exon_7], yend = -1.90, text=paste("Exon VII")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_8], 
                                  xmax=stop[exon_8], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_8], y = -2.30, xend = stop[exon_8], yend = -2.30, text=paste("Exon VIII")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_8], y = -2.1, xend = stop[exon_8], yend = -2.1, text=paste("Exon VIII")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_8], y = -1.90, xend = stop[exon_8], yend = -1.90, text=paste("Exon VIII")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_8h], 
                                  xmax=stop[exon_8h], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_8h], y = -2.30, xend = stop[exon_8h], yend = -2.30, text=paste("Exon VIIIh")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_8h], y = -2.1, xend = stop[exon_8h], yend = -2.1, text=paste("Exon VIIIh")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_8h], y = -1.90, xend = stop[exon_8h], yend = -1.90, text=paste("Exon VIIIh")), col="transparent")
p <- p + geom_rect(data=base, aes(xmin=start[exon_9], 
                                  xmax=stop[exon_9], ymin=-2.30, ymax=-1.90),  
                   colour = "black", fill = "grey")
p <- p + geom_segment(data = base, aes(x = start[exon_9], y = -2.30, xend = stop[exon_9], yend = -2.30, text=paste("Exon IX")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_9], y = -2.1, xend = stop[exon_9], yend = -2.1, text=paste("Exon IX")), col="transparent")
p <- p + geom_segment(data = base, aes(x = start[exon_9], y = -1.90, xend = stop[exon_9], yend = -1.90, text=paste("Exon IX")), col="transparent")

# Add promoters (note these are not specifically labeled in hover text)
p <- p + geom_rect(data=base, aes(xmin=start[promoter_1], 
                                  xmax=stop[promoter_1], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")
p <- p + geom_rect(data=base, aes(xmin=start[promoter_2], 
                                  xmax=stop[promoter_2], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")
p <- p + geom_rect(data=base, aes(xmin=start[promoter_3], 
                                  xmax=stop[promoter_3], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")
p <- p + geom_rect(data=base, aes(xmin=start[promoter_4], 
                                  xmax=stop[promoter_4], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")
p <- p + geom_rect(data=base, aes(xmin=start[promoter_5], 
                                  xmax=stop[promoter_5], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")
p <- p + geom_rect(data=base, aes(xmin=start[promoter_6], 
                                  xmax=stop[promoter_6], ymin=-2.25, ymax=-1.95),  
                   colour = "black", fill = "blue")

# Add TSS 
p <- p + geom_segment(data = tss, aes(x = x, y = y, xend = xend, yend = yend, text=paste("TSS", "\n", "hg38: ", xend)), col="orange")

# Add promoter/exon/TSS/CpG legend
p <- p + geom_rect(data=base, aes(xmin=start[gene]+1000, 
                                  xmax=start[gene]+4000, ymin=1, ymax=1.3, text="here is some text"),  
                   colour = "black", fill = "blue")
p <- p + geom_text(data = base, aes(x = start[gene]+8500, y = 1.15, label = "Promoter"),  hjust = 0, size = 4)
p <- p + geom_rect(data=base, aes(xmin=start[gene]+1000, 
                                  xmax=start[gene]+4000, ymin=0.5, ymax=0.8),  
                   colour = "black", fill = "grey")

# Add transcription direction arrow 
p <- p + geom_segment(data = base, aes(x = start[gene]+30000, y = -3.6, xend = start[gene]+27000, yend = -3.6))
p <- p + geom_segment(data = base, aes(x = start[gene]+27000, y = -3.6, xend = start[gene]+27500, yend = -3.5))
p <- p + geom_segment(data = base, aes(x = start[gene]+27000, y = -3.6, xend = start[gene]+27500, yend = -3.7))
p <- p + geom_text(data = base, aes(x = start[gene]+29000, y = -3.9, label = "Transcription direction (R to L)"),  hjust = 0, size = 4)

# Other info (legend)
p <- p + geom_text(data = base, aes(x = start[gene]+8500, y = 0.65, label = "Exon     "),  hjust = 0, size = 4)
p <- p + geom_segment(data = base, aes(x = start[gene]+20000, y = 1, xend = start[gene]+20000, yend = 1.3))
p <- p + geom_segment(data = base, aes(x = start[gene]+20000, y = 0.5, xend = start[gene]+20000, yend = 0.8), col="orange")
p <- p + geom_text(data = base, aes(x = start[gene]+25000, y = 1.15, label = "CpG site"),  hjust = 0, size = 4)
p <- p + geom_text(data = base, aes(x = start[gene]+25000, y = 0.65, label = "TSS     "),  hjust = 0, size = 4)

# Set axis/labels
p_base <- p + labs(y=" ", x = "Genome Reference Consortium Human Build 38 (hg38) Chr11 position") + scale_x_continuous(breaks = round(seq(base$start[gene]-3509, base$stop[gene]+29991, by = 10000),1)) 

# Add CpG segments for all CpGs that have been identified 
p2 <- p_base + geom_segment(data = cpgs, aes(x = x, y = y, xend = xend, yend = yend, text=paste("hg38: ", xend, "\n", "Illumina probe: ", Illumina)), col="black") 

# Wrap in ggplotly 
p_base_all <- ggplotly(p2, tooltip="text")

# Plot 1 = p_base: Plot with only base drawing (no CpGs) (layers added for user-defined figure on page 3 of app)
# Plot 2 = p_base_all: Plots with all CpGs with known positions identified as part of review (used on page 2 of app)

############ CREATE SEARCHABLE TABLE
# Modify the long data base so that it is prettier/searchable by the user
df2 <- df
x <- which(names(df2)=="Position")
names(df2)[x] <- "Position (hg38)"
x <- which(names(df2)==("Stop"))
names(df2)[x] <- "Stop Position (hg38)"
x <- which(names(df2)=="Tissue_celltype")
names(df2)[x] <- "Tissue_details"
df2$Illumina <- NULL
df2$x <- NULL
df2$y <- NULL
df2$xend <- NULL
df2$yend <- NULL
df2$y2 <- NULL
df2$y2_multi <- NULL
df2$Old_RefID <- NULL
x <- which(names(df2)=="hg19.Position")
names(df2)[x] <- "Position (hg19)"
x <- which(names(df2)=="hg19.Stop")
names(df2)[x] <- "Stop Position (hg19)"
x <- which(names(df2)=="First.Author")
names(df2)[x] <- "First Author (Year)"
# Change code value from 1/2/3 to Significant, Non-significant, NA
df2$Code[df2$Code==1] <- "NS"
df2$Code[df2$Code==2] <- "Sig"
df2$Code[df2$Code==3] <- "Not examined individually"
head(df2)
table_search <- df2 
table_search$RefID <- as.factor(table_search$RefID)
rm(df2)


# PART 2: CONSTRUCT BDNF DNA Methylation (DNAm) Brain Map 
## ui
ui <- dashboardPage(
  
  # Title ----
  dashboardHeader(title = "BDNF DNAm Map"),
  
  # Sidebar ---- 
  dashboardSidebar(
    sidebarMenu(id = "sidebarid",
                menuItem("Overview", tabName = "page1"),
                menuItem("Base plot", tabName = "page2"),
                menuItem("User-defined plot", tabName = "page3"),
                menuItem("Table", tabName = "page4"),
                menuItem("Developers", tabName = "page5"),
                # Add conditional panel with user input options (as these are only applicable 
                # to the user-defined plot on page 3)
                conditionalPanel(
                  'input.sidebarid == "page3"',
                  pickerInput("ref","Select the refID(s) of interest.", choices=refIDs, selected=refIDs, options = list(`actions-box` = TRUE), multiple = T),
                  pickerInput("cpg","Select the CpG(s) of interest by hg38 position.", choices=cpgs_names, selected=cpgs_names, options = list(`actions-box` = TRUE), multiple = T),
                  #checkboxGroupInput("phenotype", "Select the brain-related phenotype(s) of interest.", phenotypes, selected=phenotypes), 
                  pickerInput("phenotype","Select the brain-related phenotype(s) of interest.", choices=phenotypes, selected=phenotypes, options = list(`actions-box` = TRUE), multiple = T),
                  #checkboxGroupInput("tissue", "Select the tissue(s) of interest.", tissues, selected=tissues), 
                  pickerInput("tissue","Select the tissue(s) of interest.", choices=tissues, selected=tissues, options = list(`actions-box` = TRUE), multiple = T),
                  radioButtons("color_by", label="Choose 'color by' factor.", choices=list("Significance" = 1, "Tissue" = 2, "Phenotype" = 3), selected=1), 
                  radioButtons("sig_cpgs", label="Include nonsignificant CpGs on figure.", choices=list("Yes" = 1, "No" = 2), selected=1)
                )
    )
  ),
  
  # Body of application ----
  dashboardBody(
    tabItems(
      # Page 1 ----
      tabItem(tabName = "page1", 
              p(a("BDNF DNA Methylation Map"), style = "font-size:30px"),
              p(a("Overview"),style = "font-size:27px"),
              p(a("Introduction"),style = "font-size:25px"),
              p("The BDNF DNA Methylation Map is designed to interactively visualize the specific positions of BDNF DNA methylation (DNAm) CpG sites investigated in association with brain-related phenotypes in humans. The data for this application were curated by Treble-Barna et al as part of a manuscript entitled Brain-Derived Neurotrophic Factor (BDNF) Epigenomic Modifications and Brain-Related Phenotypes in Humans: A Systematic Review. The paper is available at https://pubmed.ncbi.nlm.nih.gov/36764636/. Source code and additional information for this application are available via the BDNF_DNAmMap GitHub repository [https://github.com/lwheinsberg/BDNF_DNAmMap]. This application is available under license GPL-2 (https://creativecommons.org/licenses/by-sa/3.0/). Copyright 2022, University of Pittsburgh. All Rights Reserved.",style = "font-size:15px"),               
              p(a("Application usage"),style = "font-size:25px"),
              p("This application can be navigated using the menu bar to the left, with pages described below.",style = "font-size:15px"), 
              p("Overview: Introduces the application and provides an overview of its use.", style = "font-size:15px"), 
              p("Base plot: Provides a to-scale depiction of the basic BDNF gene structure (i.e., exons, promoters, transcription start sites) as well as all BDNF-associated CpG sites examined in studies included in our systematic review for which positions could be identified.", style = "font-size:15px"), 
              p("User-defined plot: Provides a customizable and interactive plot synthesizing the results of the systematic review. Customizable settings include: (1) article reference ID (corresponding to the systematic review); (2) CpG position (Genome Reference Consortium Human Build 38 [hg38]); (3) broad brain-related phenotype; (4) tissue type; (5) ‘color-by’ factor allowing the user to color the study-specific findings by statistical significance, tissue type, or broad phenotype; and (6) an option to display only CpGs that were statistically significant in association with brain-related phenotypes per the study authors’ definition.",style = "font-size:15px"), 
              p("Table: Provides a searchable data base that corresponds to the systematic review findings and user-defined plot.",style = "font-size:15px"), 
              p("Developers: Provides names and contact information for the application developers.",style = "font-size:15px"),
              p("More details can be found within the application pages. Please note the figures in this application were created using ggplotly() which allows users to zoom into different sections and hover over objects to display text with exon label or CpG site position.",style = "font-size:15px"),
              p(a("Video instructions"),style = "font-size:25px"),
              p("A brief video tutorial has been created and is available at https://www.youtube.com/watch?v=1qEwysAbmCQ.",style = "font-size:15px"),               
      ),
      # Page 2 ----
      tabItem(tabName = "page2", 
              p(a("BDNF DNA Methylation Map"),style = "font-size:30px"),
              p(a("Base plot"),style = "font-size:27px"),
              p("The plot below displays the positions of each CpG site (black vertical lines) examined within each study for which this piece of data was available. CpG sites are shown in relation to 11 BDNF exons (gray boxes), 6 promoter regions (blue boxes), and 32 transcription start sites (gold vertical lines) according to Pruunsild and colleagues (Pruunsild, P., Kazantseval, A., Aid, T., Palm, K., Timmusk, T., 2007. Dissecting the human BDNF locus: Bidirectional transcription, complex splicing, and multiple promoters. Genomics 90. https://doi.org/10.1016/j.ygeno.2007.05.004). All data are mapped to the Genome Reference Consortium Human Build 38 (hg38). Please note the plot below was created using ggplotly() which allows users to zoom into different sections and hover over objects to display text with exon label or CpG site hg38 chr11 position. Promoter labels are not listed. Direction of transcription is right to left.", style="font-size:15px"),
              br(), br(),
              # Add progress spinner
              withSpinner(
                plotlyOutput("BDNF_plot1"))),
      # Page 3 ----
      tabItem(tabName = "page3", 
              p(a("BDNF DNA Methylation Map"),style = "font-size:30px"),
              p(a("User-defined plot"),style = "font-size:27px"),
              p("The plot below allows the user to customize a figure using the navigation bar to the left. To customize a plot, select the RefID(s), CpG(s), phenotype(s), and tissue(s) of interest, as well as a 'color by' and CpG statistical significance preference. The plot will automatically update based on the user-defined inputs. Users can zoom in to view specific regions of the figure in more detail. Please note there may be some delays depending on internet speed/selection.", style = "font-size:15px"),
              p("When the 'color by' option is set to 'CpG significance', study x CpG x phenotype x tissue-specific results are represented by a black circle (non-significant in association with brain-related phenotype) or red triangle (significant in association with brain-related phenotype). Shapes located above CpG sites (vertical black line) represent results for an association test between a single CpG site and a brain-related phenotype. Shapes located below CpG sites represent results for an association test between the mean value of many CpG sites and a brain-related phenotype, with the shape positioned at the center position of CpGs examined. Taller towers of shapes represent a greater number of studies that examined a specific CpG site or region. Shapes can also be colored by either phenotype or tissue type. Users can hover over shapes to view a pop-up box containing the hg38 chromosome 11 position(s) (e.g. 27,723,268), broad brain-related phenotype, more specific phenotype as applicable, tissue and cell type from which DNA was extracted, the reference number for the study in Tables 2-14 of our systematic review manuscript, and the label used to refer to the CpG site (e.g. CpG1) by the authors of the study, if applicable. Clicking on each shape will navigate the user to the website containing the abstract and/or full text of the study. Direction of transcription is right to left.", style = "font-size:15px"),
              br(), br(),
              # Add progress spinner
              withSpinner(
                plotlyOutput("BDNF_plot"))),  
      # Page 4 ----
      tabItem(tabName = "page4",
              p(a("BDNF DNA Methylation Map"),style = "font-size:30px"),
              p(a("Table"),style = "font-size:27px"),
              p("The table below displays all data curated as part of the systematic review for which CpG positions could be identified. This table is free-text searchable (e.g., search by hg19 or hg38 CpG position, probe number, phenotype, etc) and filters by row.", style = "font-size:15px"),
              p("RefID=reference identifier, corresponding to the systematic review which accompanies this application; Study_Probe_ID=study-specific CpG identifier; Position (hg38)=CpG position in Genome Reference Consortium Human Build 38; Stop Position (hg38)=applicable only for studies that looked at mean DNA methylation across many sites (vs. a single site) and indicates the stopping point of the range of CpGs examined (reported in Genome Reference Consortium Human Build 38); Position (hg19)=CpG position in Genome Reference Consortium Human Build 37 (common synonym for this build is hg19); Stop Position (hg19)=mean DNA methylation stopping point in Genome Reference Consortium Human Build 37; Tissue=study-specific tissue examined; Phenotype=broad phenotype category examined (corresponds to Tables 2-14 of the systematic review); Phenotype_details=more nuanced phenotype description; Code=study significance based on the author's definition; URL=link to reference.", style = "font-size:15px"),
              br(), br(),
              # Add progress spinner
              withSpinner(
                # ORIGINAL : DT::dataTableOutput("mytable"))),  
                fluidRow(
                  column(width = 12, align = "right", DT::dataTableOutput("mytable"),style = "height:500px; overflow-y: scroll;overflow-x: scroll;")
                ))),
      # Page 5 ----
      tabItem(tabName = "page5", 
              p(a("BDNF DNA Methylation Map"),style = "font-size:30px"),
              p(a("Developers"),style = "font-size:27px"),
              p(a("Lacey W. Heinsberg, PhD, RN", href="https://github.com/lwheinsberg", target="_blank"),style = "font-size:18px"),
              p(a("Daniel E. Weeks, PhD", href="https://publichealth.pitt.edu/home/directory/daniel-e-weeks", target="_blank"),style = "font-size:18px"),
              p("For additional help or to submit feedback or bug reports, please contact Lacey Heinsberg at law145@pitt.edu.", style = "font-size:15px"),
              p("Note this application could not function without the accompanying database curated by Amery Treble-Barna and her team which included Lacey Heinsberg, Tara Davis, Stephen Breazeale, Zachary Stec, Aboli Kesbhat, Julia Tefs, and Ansuman Chattopadhyay. This curation took a great deal of collaborative effort and many many many hours of work. The paper that summarizes this effort can be found at https://www.sciencedirect.com/journal/neuroscience-and-biobehavioral-reviews.", style = "font-size:15px"),
      )
    )
  )
)

## Server 
server <- function(input, output) {
  
  # Page 2 ----
  # Base plot with all CpGs identified in systematic review 
  output$BDNF_plot1 <- renderPlotly({
    p_base_all
  })
  
  # Page 3 ----
  # User-defined plot 
  output$BDNF_plot <- renderPlotly({
    
    # Store user-defined input options as reactive functions
    # -------     
    # RefID(s) of interest
    refid <- reactive({
      input$ref
    }) %>% debounce(4000) # Use debounce to add a 4 second delay 
    # The debounce is helpful because it transforms a reactive expression by
    # preventing its invalidation signals from being sent unnecessarily often
    # This lets you ignore a very "chatty" reactive expression until it 
    # becomes idle, which is useful when the intermediate values don't 
    # matter as much as the final value, and the downstream calculations
    # that depend on the reactive expression take a long time
    
    # CpG(s) of interest
    cpg_s <- reactive({
      input$cpg
    }) %>% debounce(4000)
    
    # Phenotype(s) of interest
    colm <- reactive({
      input$phenotype
    }) %>% debounce(4000) 
    
    # Tissue(s) of interest
    tissue <- reactive({
      input$tissue
    }) %>% debounce(4000)
    
    # Color by preference 
    color_by <- reactive({
      input$color_by
    }) %>% debounce(4000)
    
    # Do you want to display non-significant CpGs on the plot? 
    sig_cpgs <- reactive({
      input$sig_cpgs
    }) %>% debounce(4000) 
    
    # ------- 
    if(sig_cpgs()==1){ # Option 1, loop logic used if the user wants to 
      # include non-significant CpGs on the plot 
      
      # Subset df (by row) for refIDs(s) selected by user 
      df2.1 <- df %>% filter_at(vars(RefID), any_vars(. %in% refid()))
      # Subset for CpG(s) selected by user 
      df2.2 <- df2.1 %>% filter_at(vars(Position), any_vars(. %in% cpg_s()))
      # Subset for phenotype(s) selected by user 
      df2 <- df2.2 %>% filter_at(vars(Phenotype), any_vars(. %in% colm()))
      # Subset for tissue(s) selected by user 
      df3 <- df2 %>% filter_at(vars(Tissue), any_vars(. %in% tissue())) 
      
      # Null plot
      # If df3 is empty, then print base plot with a message to user saying no CpGs meet the criteria they selected 
      if (dim(df3)[1]==0) {
        p_null <- p_base + annotate("text", x=27683893, y=0, label= "No CpGs meet search criteria") 
        p_3_phen <- ggplotly(p_null)
        
      } else {
        # -------        
        # Note: The data base for this application was curated using three codes
        # Code==1 means a CpG site was examined, but not significantly associated with a phenotype 
        # Code==2 means a CpG site was examined, and was significantly associated with a phenotype 
        # Code==3 means a CpG site was included in a study, but was not examined for association with
        # a phenotype; instead, code 3 sites were examined in combination with other sites as 
        # mean methylation or methylation ratio values
        
        # In the case of a mean/ratio study, we still want to plot the individual CpG sites on the 
        # graph, but we only want to include a study-specific dot for the multi-CpG results
        # To distinguish this difference on the plot, I will plot single site CpG results ABOVE
        # the gene, and multi-site resutls BELOW the gene
        # This will be done by layering the results on the figure using different approaches 
        
        # Subset table (by row) based on (1) single sites only and (2) mean/ratio sites 
        df4_single <- subset(df3, is.na(Stop)) # Data frame for single CpG sites
        df4_multi <- subset(df3, !is.na(Stop)) # Data frame of mean/ratio sites
        
        # Add CpG segments to the plot if df4_single is not empty
        if (dim(df4_single)[1] != 0) {
          # ------- Base plot with CpGs 
          
          # This builds a base plot with CpG segments of interest (both single sites and mean/ratio sites)
          # I will add study-specific dots to this plot below 
          p_phen_base <- p_base + geom_segment(data = df4_single, aes(x = x, y = y, xend = xend, yend = yend, text=paste("hg38: ", xend, "\n", "Illumina ID: ", Illumina))) 
          
          # Now that we have layered the CpG onto the plot, are done with rows Coded as 3
          # Change values of 3 to 0 
          df4_single$Code[df4_single$Code==3] <- 0
          # Filter out 0 Codes 
          df5_single2 <- df4_single %>% filter(Code != 0)
          # Index the CpGs for vertical placement/organization on the plot, creating an 'n' variable 
          # that is indexes the CpG position 
          if (dim(df5_single2)[1] != 0) {
            df5_single <- df5_single2 %>% group_by(Position) %>% mutate(n = factor(1:n()))
          } else {
            df5_single <- df5_single2
          }
          
          # -------   Layer on study-specific dots (study x CpG x tissue x phenotype) 
          
          # Recode 'Code' as a factor with a label for use in ggplot
          df5_single$Code <- factor(df5_single$Code, levels = c("1", "2"), 
                                    labels = c("Nonsignificant", "Significant"))
          
          # ------    # Add conditional logic for user-defined color codes         
          # -----     # Color by significance 
          if (color_by()==1) {
            # Single site results
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, shape= Code, color= Code, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_manual(values=c("Nonsignificant"="Black", "Significant"="Red")) +scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17)) + labs(color="Statistical\nSignificance", shape="Statistical\nSignificance") 
            } else {
              p_phen_f1 <- p_phen_base
            }
            
            # Multi site results 
            if (dim(df4_multi)[1] != 0) {
              # Create a midpoint value for placement of mean/ratio dot on plot (below gene)
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              # Index the CpGs for vertical placement/organization on the plot 
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              # Recode as factor
              df5_multi$Code <- factor(df5_multi$Code, levels = c("1", "2"), 
                                       labels = c("Nonsignificant", "Significant"))
              # Layer dots on plot 
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, shape= Code, color= Code,customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_manual(values=c("Nonsignificant"="Black", "Significant"="Red")) +scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17))  + labs(color="Statistical\nSignificance", shape="Statistical\nSignificance") + guides(size = FALSE)
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
          
          # -----   # Color by Tissue
          if (color_by()==2) {
            # Single site results
            df5_single$Tissue <- as.factor(df5_single$Tissue)
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, shape= Code, color= Tissue, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17)) + labs(color="Tissue") + guides(size = FALSE, shape=FALSE)
            } else {
              p_phen_f1 <- p_phen_base
            }
            # Multi site results 
            if (dim(df4_multi)[1] != 0) {
              # Create midpoint for placement of dot 
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              # Index the CpGs for vertical placement/organization on the plot 
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              df5_multi$Code <- factor(df5_multi$Code, levels = c("1", "2"), 
                                       labels = c("Nonsignificant", "Significant"))
              df5_multi$Tissue <- as.factor(df5_multi$Tissue)
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, shape= Code, color= Tissue, customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17)) + labs(color="Tissue") + guides(size = FALSE, shape=FALSE) 
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
          
          # -----   # Color by Phenotype
          if (color_by()==3) {
            # Single site results
            df5_single$Phenotype <- as.factor(df5_single$Phenotype)
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, shape= Code, color= Phenotype, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17)) + labs(color="Phenotype") + guides(shape = FALSE)
            } else {
              p_phen_f1 <- p_phen_base
            }
            # Multi site results 
            if (dim(df4_multi)[1] != 0) {
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              # Index the CpGs for vertical placement/organization on the plot 
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              df5_multi$Code <- factor(df5_multi$Code, levels = c("1", "2"), 
                                       labels = c("Nonsignificant", "Significant"))
              df5_multi$Phenotype <- as.factor(df5_multi$Phenotype)
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, shape= Code, color= Phenotype, customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_shape_manual(values=c("Nonsignificant"=16, "Significant"=17)) + labs(color="Phenotype") + guides(size = FALSE, shape=FALSE) #, size="Statistical\nSignificance") 
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
        }
        
        # ------
        # Wrap in plotly 
        p_3_phen <- ggplotly(p_phen_fin, tooltip="text")
        
        # ------
        #########################
        # Function for opening URL
        p_3_phen <- htmlwidgets::onRender(p_3_phen, "
          function(el, x) {
            el.on('plotly_click', function(d) {
              var URL = d.points[0].customdata;
              //URL
              window.open(URL);
            });
          }
        ")
        # End function for opening URL
        #########################
      } 
      
    } # End condition cpg_sig==1 loop 
    
    # ------     
    else { # Option 2, if the user does not want to include non-significant CpGs on the plot 
      
      # Repeat above logic, but pull out only significant studies 
      # This section of the server and below mimics above, so will have fewer annotations
      # For annotation details see above loop
      df2.1 <- df %>% filter_at(vars(RefID), any_vars(. %in% refid()))
      df2.2 <- df2.1 %>% filter_at(vars(Position), any_vars(. %in% cpg_s()))
      df2 <- df2.2 %>% filter_at(vars(Phenotype), any_vars(. %in% colm()))
      df3 <- df2 %>% filter_at(vars(Tissue), any_vars(. %in% tissue())) 
      df4 <- df3 %>% filter(Code==2|Code==3) # Select only significant CpG
      
      if (dim(df4)[1]==0) {
        p_null <- p_base + annotate("text", x=27683893, y=0, label= "No CpGs meet search criteria") 
        p_3_phen <- ggplotly(p_null)
        
      } else {
        
        df4_single <- subset(df4, is.na(Stop)) # Data frame for single CpG sites
        df4_multi <- subset(df4, !is.na(Stop)) # Data frame of mean methylation 
        
        # ------ # Add CpG segments to the base plot 
        if (dim(df4_single)[1] != 0) {
          p_phen_base <- p_base + geom_segment(data = df4_single, aes(x = x, y = y, xend = xend, yend = yend, text=paste("hg38: ", xend, "\n", "Illumina ID: ", Illumina))) 
          
          df5_single2 <- df4_single %>% filter(Code==2) # Select only significant CpG (ignoring codes of 3)
          if (dim(df5_single2)[1] != 0) {
            df5_single <- df5_single2 %>% group_by(Position) %>% mutate(n = factor(1:n()))
          } else {
            df5_single <- df5_single2
          }
          df5_single$Code <- factor(df5_single$Code, levels = c("2"), 
                                    labels = c("Significant"))
          
          # ------           
          if (color_by()==1) {
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, color= Code, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_manual(values=c("Red")) +  labs(color="Statistical\nSignificance")
            } else {
              p_phen_f1 <- p_phen_base
            }
            if (dim(df4_multi)[1] != 0) {
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              df5_multi$Code <- factor(df5_multi$Code, levels = c("2"), 
                                       labels = c("Significant"))
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, color= Code, customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_manual(values=c("Red"))  + labs(color="Statistical\nSignificance") + guides(size=FALSE)
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
          
          if (color_by()==2) {
            df5_single$Tissue <- as.factor(df5_single$Tissue)
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, color= Tissue, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + labs(color="Tissue")
            } else {
              p_phen_f1 <- p_phen_base
            }
            if (dim(df4_multi)[1] != 0) {
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              df5_multi$Code <- factor(df5_multi$Code, levels = c("2"), 
                                       labels = c("Significant"))
              df5_multi$Tissue <- as.factor(df5_multi$Tissue)
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, color= Tissue, customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ",Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + labs(color="Tissue")+ guides( size = FALSE) 
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
          
          if (color_by()==3) {
            df5_single$Phenotype <- as.factor(df5_single$Phenotype)
            if (dim(df5_single2)[1] != 0) {
              p_phen_f1 <- p_phen_base + geom_point(data=df5_single, aes(x=Position, y=y2, group=n, color= Phenotype, customdata=URL, text=paste("hg38: ", Position, "\n", "Phenotype: ", Phenotype, "\n", "Tissue: ", Tissue, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + scale_color_viridis(discrete = TRUE, option = "D") + labs(color="Phenotype") 
            } else {
              p_phen_f1 <- p_phen_base
            }
            if (dim(df4_multi)[1] != 0) {
              df4_multi$mid <- (df4_multi$Position + df4_multi$Stop)/2
              df5_multi <- df4_multi %>% group_by(mid) %>% mutate(n = factor(1:n()))
              df5_multi$Code <- factor(df5_multi$Code, levels = c("2"), 
                                       labels = c("Significant"))
              df5_multi$Phenotype <- as.factor(df5_multi$Phenotype)
              p_phen_fin <- p_phen_f1 + geom_point(data=df5_multi, aes(x=mid, y=y2_multi, group=n, size=10, color=Phenotype, customdata=URL, text=paste("Multiple CpGs \n", "hg38: ", Position, "-", Stop, "\n", "Phenotype: ", Phenotype, "\n", "Phenotype details: ", Phenotype_details, "\n", "Tissue: ", Tissue, "\n", "Tissue details: ", Tissue_celltype, "\n", "RefID: ", RefID, "\n", "First Author (Year): ", First.Author, "\n", "Study-specific PID: ", Study_Probe_ID)), position=ggstance::position_dodgev(height=0.3)) + labs(color="Phenotype")+ guides( size = FALSE)
            } else {
              p_phen_fin <- p_phen_f1
            }
          }
          
        }
        
        p_3_phen <- ggplotly(p_phen_fin, tooltip="text")
        
        #########################
        # Function for opening URL
        p_3_phen <- htmlwidgets::onRender(p_3_phen, "
          function(el, x) {
            el.on('plotly_click', function(d) {
              var URL = d.points[0].customdata;
              //URL
              window.open(URL);
            });
          }
        ")
        # End function for opening URL
        #########################
        
        #p_3_phen <- p_3_phen %>% toWebGL() # To speed up 
        #p_3_phen <- p_3_phen %>% partial_bundle() # To speed up
        
      } 
    } # End condition cpg_sig==2 loop
    
    # Page 3 ---- 
    # User-defined plot 
    p_3_phen
  }) # End render plotly element 
  
  # Page 4 ----
  # Searchable table
  output$mytable <- DT::renderDataTable({
    DT::datatable(table_search, filter="top", rownames=FALSE)  
  })
} # Close server 

## Deploy application 
shinyApp(ui, server)

