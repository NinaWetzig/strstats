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
              help="aggregates data frames and sums up the column 'Reads' per 'Marker' and per 'Marker' and 'Allele'.output = data frame", metavar="character"),
  make_option(c("-u", "--Run01"), type="character", default=NULL, 
              help="tabel with samples, marker", metavar="character"),
  make_option(c("-f", "--profilefile"), type="character", default=NULL, 
              help="profiles tabel with patients", metavar="character"),
  make_option(c("-g", "--global_sample"), type="character", default=NULL, 
              help="table with patients and samples", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)

if(FALSE){
#Tabelle mit den Profilen
#if (!is.null(opt$profile)){
  #print("Profiles:")
  #print(read_profile(opt$profile))
#} else {}

# Input-Tabelle zu beliebigem Sample
#if (!is.null(opt$jinput)){
  #print("J_input")
  #read_Jinput(opt$jinput)
#} else {}

#Kombiniert zwei Samples in einer csv-Tabelle
#if (!is.null(opt$jinput1) & !is.null(opt$jinput2)){ 
  #print("csv")
  #create_csv(opt$jinput1, opt$jinput2)
#} else {}

#Berechnet Read-Summen
#if (!is.null(opt$rawinput)){
  #print("Tabellen:")
  #print(read_raw_input(opt$rawinput))
#} else{}

#source/strstats.R -p resources/profiles.csv  -j data/Run_01/variant_calling/J1_STR.out -e data/Run_01/variant_calling/J1_STR.out -z data/Run_01/variant_calling/J2_STR.out -r /home/nina/Downloads/Run_01_short.out
}

#Für die Funktion zur Zuordnung der Marker, Samples und Patienten und zur Abgleichung der callAllele mit den korrekten Allelen

#erstellt combined-Table aus der Run Tabelle und der global-Sample Tabelle
if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample) & !is.null(opt$output)){
  out_dir <- paste(opt$output, "/", sep = "")
  dir.create(out_dir, recursive = TRUE)
  #combine_Allele(opt$Run01, opt$profilefile, opt$global_sample)
  combined_table <- combine_Allele(opt$Run01, opt$profilefile, opt$global_sample, out_dir)
  #print("combine_Allele")
  #compare_Allele(combined_table)
  SN_LN_stutter_true_table <- compare_Allele(combined_table, out_dir)
  #print("compare_Allele")
  Doc_SCR_table <- DoC_SCR_StR(combined_table, SN_LN_stutter_true_table, out_dir)
  #print("DoC_SCR_StR")
  #ACR_function(SN_LN_stutter_true_table)
  ACR_table <- ACR_function(SN_LN_stutter_true_table, out_dir)
  #print("ACR")
  create_plots(ACR_table, combined_table, out_dir)
  #print("plots")
  create_pies(SN_LN_stutter_true_table, out_dir)
  #print("pies")
  
  print("combineAllele")
  #combine_Allele(opt$Run01, opt$profilefile, opt$global_sample)
} else {}

if(FALSE) {
  #Vergleicht die call Allele mit den korrekten Allelen und schreibt die Resultate in eine neue Spalte
#if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample)){ 
#  print("compareAllele")
#  compare_Allele(opt$Run01, opt$profilefile, opt$global_sample)
#} else {}

#Berechnet DoC und SCR und StR
#if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample)){ 
#  print("DoC-SCR")
#  DoC_SCR_StR(opt$Run01, opt$profilefile, opt$global_sample)
#} else {}

#Berechnet ACR
#if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample)){ 
  #print("ACR")
  #ACR_function(opt$Run01, opt$profilefile, opt$global_sample)
#} else {}

#erstellt plots
#if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample)){ 
  #print("plots")
 # create_plots(opt$Run01, opt$profilefile, opt$global_sample)
#} else {}

#erstellt pie_charts
#if (!is.null(opt$Run01) & !is.null(opt$profilefile) & !is.null(opt$global_sample)){ 
 # print("pies")
 # create_pies(opt$Run01, opt$profilefile, opt$global_sample)
#} else {}

#source/strstats.R -u /home/nina/projects/strstats/data/Run_01.out -f resources/profiles.csv -g resources/global_sample_overview_forensic.csv
}