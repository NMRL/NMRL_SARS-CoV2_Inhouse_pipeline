#!/mnt/home/groups/nmrl/cov_analysis/SARS-CoV2_assembly/tools/rbase_env/bin/Rscript

#CHECK IF PACKAGES ARE INSTALLED, IF NOT - INSTALL
packageList <- c("scales","lubridate","tidyverse","ggplot2","dplyr","grDevices", "plotly", "htmlwidgets")
newPackages <- packageList[!(packageList %in% installed.packages()[,"Package"])]
if(length(newPackages)) install.packages(newPackages, repos = c("http://cran.us.r-project.org"))
 
#LOAD PACKAGES
library(grDevices)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(lubridate)
library(scales)
library(tsibble)
library(plotly)
library(htmlwidgets)


#COMMAND-LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)

if (length(args)<1) {
  stop("
  First argument - path to filtered_mutations.csv
  ", call.=FALSE)
} 

#PATHS & PRE-PROCESSING
path <- args[1]
excel_data <- read.csv(path,  encoding = "UTF-8")#SUMMARY DATA
excel_data$sampling_date <- as.Date(excel_data$sampling_date, format="%Y-%m-%d", origin="1970-01-01") #CASTING DATE FROM CHAR TO DATE.OBJECT
mutation_frequency <- excel_data$Fill #RENAMING MUTATION FREQUENCY TOOLTIP
mutation <- excel_data$Mutation #RENAMING FOR TOOLTIP
excel_data$n <- factor(excel_data$sampling_date) #TO REORDER BY SAMPLING DATE
title <- paste("Sars-CoV2 S gēnā novēroto mutāciju biežums pēc paraugu paņemšanas datuma.\nApkopoti dati par mutācijām, kuras novērotas vismaz divos secīgos paraugu paņemšanas datumos un sastopamas iekš 65% lasījumu.")

#GENERATING HEATMAP PLOT
currentDate <- Sys.Date()
data <- excel_data %>% mutate(Mutation = fct_reorder(Mutation, desc(n)), text = paste0("datums: ", sampling_date, "\n", "mutācija: ", Mutation, "\n", "paraugu daļa: ",100*round(mutation_frequency,2),"%\nparaugu kopskaits: ", Count)) #ADDING TOOLTIP COLUMN AS DEFINED IN https://r-graph-gallery.com/79-levelplot-with-ggplot2.html
# grDevices::cairo_pdf(paste('./mut_heatmap_data/mutāciju_apkopojums_',currentDate,".pdf",sep=""), width=50, height=25*length(excel_data['Mutation'])) #TO PRODUCE STATIC PDF
heatmap <- data %>% ggplot(aes(sampling_date, Mutation, fill=mutation_frequency, text=text)) + 
  geom_tile(colour = "grey") + #OUTLINE COLOR OF HEATMAP CELLS
  labs(title=title,
        x = "Paraugu paņemšanas nedēļa", #X-AXIS LABEL
        y=element_blank(), #NO LABEL
        fill = element_blank()) + #LEGEND LABEL
        theme(
          plot.title = element_text(hjust = 0.5, vjust=0), #PLOT TITLE
          legend.position='none',
          axis.text.x = element_text( #FORMAT X-AXIS LABEL TEXT
            angle = 90, #ROTATE COUNTER-CLOCKWISE
            vjust=0.2, #OPTIMAL VERTICAL POSITIONING RELATIVE TO AXIS TICK
            hjust = 0.95,#OPTIMAL HORIZONTAL POSITIONING RELATIVE TO AXIS TICK
            size = 8), #FONT SIZE
            panel.grid.major.y = element_line(colour = 'grey'), #ADD COLOUR TO THE GRID LINES
            panel.grid.major.x = element_blank(), #REMOVE VERTICAL GRID LINES
            panel.background = element_blank(), #REMOVE BACKGROUND GRID PANEL COLOUR
            axis.line = element_line(colour = "black"),#CHANGE AXIS COLORS
            axis.title.x = element_text(vjust=-1.5), #ADD SPACE BETWEEN TICK LABELS AND X AXIS LABEL
            axis.title.y = element_text(vjust=3), #ADD SPACE BETWEEN TICK LABELS AND Y AXIS LABEL
            plot.margin = unit(c(2,0.5,1.5,0.5), "cm"), #ADD MARGINS TO THE PLOT
          ) + scale_fill_gradient(
            low="grey", high="#EB8955", #CONTROL OVER GRADIENT COLOR
            guide = guide_colorbar(
              label = TRUE,
              draw.ulim = TRUE, 
              draw.llim = TRUE,
              frame.colour = "black",
              ticks.colour = "black",
              ticks = TRUE, 
              nbin = 10))+
        scale_x_date(
          date_breaks = "1 week", #MAYOR BRAKES ON DATE-BASED X-AXIS
          date_minor_breaks = "1 week", #MINOR BRAKES ON DATE-BASED X-AXIS
          date_labels = "%Wned. /%Y") #TICK LABEL FORMAT ON DATE-BASED X-AXIS
# dev.off()
p <- ggplotly(heatmap, tooltip='text', height=ceiling(0.3*nrow(data)))#INTERACTIVITY + RESIZING: https://stackoverflow.com/questions/37241985/resize-plotly-r-ggplotly
saveWidget(p, paste('/mnt/home/groups/nmrl/cov_analysis/mutation_heatmap/mut_heatmap_datamutāciju_apkopojums_',currentDate,".html",sep="")) #SAVE INTERACTIVE PLOT