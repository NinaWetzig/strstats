#!/usr/bin/env Rscript

source("source/read_data.R")


library("optparse")

option_list = list(
  make_option(c("-p", "--profile"), type="character", default=NULL, 
              help="output = data frame", metavar="character"),
  make_option(c("-j", "--jinput"), type="character", default=NULL, 
              help="invisible output = data frame", metavar="character"),
  make_option(c("-e", "--jinput1"), type="character", default=NULL, 
              help="combines two data frames and adds columns 'j1' ans 'j2'. output = csv Table", metavar="character"),
  make_option(c("-z", "--jinput2"), type="character", default=NULL, 
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-r", "--rawinput"), type="character", default=NULL, 
              help="aggregates data frames and sums up the column 'Reads' per 'Marker' and per 'Marker' and 'Allele'.output = data frame", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

print("Lese datei:")

#Tabelle mit den Profilen
if (!is.null(opt$profile)){
  print("Profiles:")
  print(read_profile(opt$profile))
} else {}

# Input-Tabelle zu beliebigem Sample
if (!is.null(opt$jinput)){ 
  print("J_input")
  read_Jinput(opt$jinput)
} else {}

#Kombiniert zwei Samples in einer csv-Tabelle
if (!is.null(opt$jinput1) & !is.null(opt$jinput2)){ 
  print("csv")
  create_csv(opt$jinput1, opt$jinput2)
} else {}

#Berechnet Read-Summen
if (!is.null(opt$rawinput)){
  print("Tabellen:")
  print(read_raw_input(opt$rawinput))
} else{}
