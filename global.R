##################################################
#                    APUAMA                      #
# a software tool for reaction rate calculations #
##################################################

# Load R packages
library(shiny)
library(shinyWidgets)
library(shinythemes)
library(shinyjs)
library(ggplot2)
library(dplyr)
library(DT)
library(Rcpp)
library(readODS)
library(htmlTable)
sourceCpp("apuama.cpp")
#update.packages(ask = FALSE, checkBuilt = TRUE) # to update current packages and also restore broken ones
