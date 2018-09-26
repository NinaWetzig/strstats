#!/usr/bin/env Rscript

source("source/read_data.R")


library("optparse")

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="tabel with samples, marker", metavar="character"),
  make_option(c("-p", "--profilefile"), type="character", default=NULL, 
              help="profiles tabel with patients", metavar="character"),
  make_option(c("-g", "--global_sample"), type="character", default=NULL, 
              help="table with patients and samples", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output directory", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
print(opt)


#Für die Funktion zur Zuordnung der Marker, Samples und Patienten und zur Abgleichung der callAllele mit den korrekten Allelen

#erstellt combined-Table aus der Run Tabelle und der global-Sample Tabelle
if (!is.null(opt$input) & !is.null(opt$profilefile) & !is.null(opt$global_sample) & !is.null(opt$output)){
  out_dir <- paste(opt$output, "/", sep = "")
  dir.create(out_dir, recursive = TRUE)
  combined_table <- combine_Allele(opt$input, opt$profilefile, opt$global_sample, out_dir)
  SN_LN_stutter_true_table <- compare_Allele(combined_table, out_dir)
  Doc_SCR_table <- DoC_SCR_StR(combined_table, SN_LN_stutter_true_table, out_dir)
  ACR_table <- ACR_function(SN_LN_stutter_true_table, out_dir)
  create_plots(ACR_table, combined_table, out_dir)
  create_pies(SN_LN_stutter_true_table, out_dir)
} else {}
