#!/usr/bin/env Rscript

source("source/read_data.R")

args = commandArgs(trailingOnly=TRUE)
print("Lese datei:")
print(args)
print("Profiles:")
print(read_profile(args[[1]]))
print("J_input")
read_Jinput(args[[2]])
print("csv")
create_csv(args[[2]], args[[3]])




