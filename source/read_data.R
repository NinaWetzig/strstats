library(reshape)
library(tidyr)
library(dplyr)

#' params: profile_file - a profile csv file
#' return: data.frame of the prfile file
read_profile <- function(profile_file){
  profiles.names <- c("Marker", "HeLa", "Glenn", "Sarah" , "Trine" , "Nicole")
  profiles <- read.csv(profile_file, na.strings=c("","NA"), sep =",", col.names = profiles.names, header = TRUE, stringsAsFactors = FALSE)
  # munging
  profiles$Marker <- gsub(pattern = "\n", replacement="", x = profiles$Marker)
  profiles$Marker <- gsub(pattern = " ", replacement="", x = profiles$Marker)
  return(profiles)
}

#Funktion zur Berechnung der Read-Summen je nach Marker und Allel
read_raw_input  <- function(raw_input_file){
  raw_names <- c("Run", "Sample", "Marker", "Allele", "Pattern", "Quality", "Reads", "Variants")
  #Einlesen der Tabelle
  raw_input <- read.csv(raw_input_file, na.strings=c("","NA"), sep =",", col.names = raw_names, header = TRUE, stringsAsFactors = FALSE)
  #Berechnet die Read-Summe pro Marker
  agg_input <- aggregate(raw_input$Reads, by = list(Run = raw_input$Run, Sample = raw_input$Sample, Marker = raw_input$Marker), FUN = function(x) sum = sum(x))
  agg_input <- do.call(data.frame, agg_input)
  #Berechnet die Read-Summe pro Marker und Allele
  agg_input_Allele <- aggregate(raw_input$Reads, by = list(Run = raw_input$Run, Sample = raw_input$Sample, Marker = raw_input$Marker, Allele = raw_input$Allele), FUN = function(x) sum(x))
  agg_input_Allele <- do.call(data.frame, agg_input_Allele)
  #neue Namen der Spalten
  colnames(agg_input) <- c("Run", "Sample", "Marker", "Sum")
  colnames(agg_input_Allele) <- c("Run", "Sample", "Marker", "Allel", "Sum")
  #result ist eine Liste der gewünschten Ausgaben, damit die Funktion mehrere Ausgaben gibt
  result = list(agg_input, agg_input_Allele)
  return(result)
}

#Einlesen der Run_01_short Tabelle
#raw_names <- c("Run", "Sample", "Marker", "Allele", "Pattern", "Quality", "Reads", "Variants")
#raw_input_file <- read.csv("/home/nina/Downloads/Run_01_short.out",na.strings=c("","NA"), sep =",", col.names = raw_names, header = TRUE, stringsAsFactors = FALSE)
#Berechenen der Read-Summe pro Marker
#agg <- aggregate(raw_input_file$Reads, by = list(run = raw_input_file$Run, raw_input_file$Sample, state = raw_input_file$Marker), FUN = function(x) sum=sum(x))

#Run_01_short Tabelle
#verbinde die Spalten Marker und Allele zu einer
#tab_combined <- unite(raw_input_file, newcol, c(Marker, Allele))
#Berechnen der Read-Summe pro Marker und Allel
#agg_Run_01_short <- aggregate(tab_combined$Reads, by = list(run = tab_combined$Run, tab_combined$Sample, state = tab_combined$newcol), FUN = function(x) sum=sum(x))


#Einlesen der Run_01 Tabelle
#raw_names <- c("Run", "Sample", "Marker", "Allele", "Pattern", "Quality", "Reads", "Variants")
#Run_01 <- read.csv("/home/nina/Downloads/Run_01.out",na.strings=c("","NA"), sep =",", col.names = raw_names, header = TRUE, stringsAsFactors = FALSE)
#Berechnen der Read-Summe pro Marker
#### NEW VERSION
#agg_marker_Run_01 <- aggregate(Run_01$Reads, by = list(Run = Run_01$Run, Sample = Run_01$Sample, Marker = Run_01$Marker), FUN = function(x) sum=sum(x))

#Run_01 Tabelle
#verbinde die Spalten Marker und Allele zu einer
#Run_01_combined <- unite(Run_01, newcol, c(Marker, Allele))
#Berechnen der Read-Summe pro Marker und Allel
#agg_allele_Run_01 <- aggregate(Run_01$Reads, by = list(Run = Run_01$Run, Sample = Run_01$Sample, Marker = Run_01$Marker, Allele = Run_01$Allele), FUN = function(x) sum=sum(x))


read_Jinput <- function(j_input_file){
  names2 <- c("Marker", "Allele" , "Pattern" , "Pattern_Quality" , "Reads", "Variants")
  j_input <- read.csv(j_input_file, na.strings=c("","NA"), sep = ",", col.names = names2, header = FALSE, stringsAsFactors = FALSE, colClasses = c("character", "character", "character", "numeric", "numeric", "character"))
  invisible(j_input)
}

read_data <- function(profile_file, j1in, j2in){
  profiles = read_profile(profile_file)
  J1_input = read_Jinput(j1in)
  J2_input = read_Jinput(j2in)
  return(profiles)
}

ignoreme <- function(){
  read_data("resources/profiles.csv", "data/Run_01/variant_calling/J1_STR (copy).csv", "data/Run_01/variant_calling/J2_STR.out")
  
  
  profiles.names <- c("Marker", "HeLa", "Glenn", "Sarah" , "Trine" , "Nicole")
  
  profiles <- read.csv("resources/profiles.csv", na.strings=c("","NA"), sep =",", col.names = profiles.names, header = TRUE, stringsAsFactors = FALSE)
  
  names2 <- c("Marker", "Allele" , "Pattern" , "Pattern_Quality" , "Reads", "Variants")
  
  J1_input <- read.csv("data/Run_01/variant_calling/J1_STR (copy).csv", na.strings=c("","NA"), sep = ",", col.names = names2, header = FALSE, stringsAsFactors = FALSE,colClasses = c("character", "character", "character", "numeric", "numeric", "character"))
  J2_input <- read.csv("data/Run_01/variant_calling/J2_STR.out", 
                       na.strings=c("","NA"), 
                       sep = ",", 
                       col.names = names2, 
                       header = FALSE, 
                       stringsAsFactors = FALSE,
                       colClasses = c("character", "character", "character", "character", "numeric", "character"))
  
  
  
  vWA <- J1_input[J1_input$Marker == "vWA", ]
  TH01 <- J1_input[J1_input$Marker == "TH01", ]
  c <- rbind(vWA, TH01)
  head(subset(J1_input, Marker %in% c("vWA", "TH01")))
  # mit subset kann die Tabelle leichter erstellt werden
}
  
create_csv <- function(j1in,j2in){
  J1_input = read_Jinput(j1in)
  J2_input = read_Jinput(j2in)
  J1_input$Sample <- "J1"
  J2_input$Sample <- "J2"
  #J1_input <- cbind(Sample = "J1" ,J1_input )
  #J2_input <- cbind(Sample = "J2" ,J2_input )
  #Beides fügt der Tabelle eine neue Spalte hinzu. Die ersten Befehle an letzer Stelle, die letzen Befehle als erste Spalte
  
  J_input <- rbind(J1_input, J2_input) 
  
  
  J_csv <- write.csv(J_input, file = "results/J_STR_3.csv", na = "", quote = FALSE, row.names = FALSE)
  return("results/J_STR_3.csv")
}


compare_call_to_Allele <- function(Run_01_file, Profiles_file, global_samples_file){
#einlesen der Run_01_short Tabelle
Run_01_short <- read.csv(Run_01_file, na.strings=c("","NA"), sep =",", header = TRUE, stringsAsFactors = FALSE)
colnames(Run_01_short) <- c("Run", "Sample", "Marker", "call_Allele", "Pattern", "Quality", "Reads", "Variants")

#einlesen der profiles Tabelle
Profiles <- read.csv(Profiles_file, na.strings=c("","NA"), sep =",", header = TRUE, stringsAsFactors = FALSE)
Profiles_melt <- melt(Profiles, id = "Marker")
colnames(Profiles_melt) = c("Marker", "Patient", "Allel")

#einlesen der global_sample Tabelle
global_samples <- read.csv(global_samples_file, na.strings=c("","NA"), sep =",", header = TRUE, stringsAsFactors = FALSE)

#zu der Run_01_short Tabelle wird die Spalte Patienten zugefügt, nach der global_samples Tabelle
Run_01_short$Patient <- global_samples$Patient[match(Run_01_short$Sample, global_samples$Sample)]
#Run_01_short$wAllel <- Profiles_melt$Allel[match(Run_01_short$Marker, Profiles_melt$Marker, Run_01_short$Patient, Profiles_melt$Patient)]

#Profiles_cast <- cast(Profiles_melt, formula = c(Marker, Patient) ~ Allel)
#Profiles_cast <- recast(Profiles_melt,formula =  Marker ~ variable, id.var = "Allel")
#Profiles_cast <- aggregate(cbind(Marker, Patient) ~ Allel, data = Profiles_melt)
#Profiles_cast <- with(Profiles_melt, tapply(Allel, list(Marker, Patient)))
#Profiles_cast <- reshape(Profiles_melt, idvar = "Allel")
#Profiles_cast <- reshape(Profiles_melt, idvar = c("Marker", "Patient"), timevar = "Allel", direction = "wide")
#Profiles_cast <- reshape(Profiles_melt, direction="wide", idvar=c("Patient"), timevar="Marker")
#Profiles_cast <- separate(Profiles_melt, col = "Allel", into = c("Allel_1", "Allel_2"))
#Profiles_cast <- aggregate(Allel ~ Marker + Patient, data= Profiles_melt, I)

#Die Tabelle Profiles_melt hatte die Allele in einer Spalte untereinander, die neue Tabelle Profiles_cast hat für jedes Alell eine Spalte
#Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt) Jede zweite Zeile der Tabelle beginnend mit der ersten Zeile wird in eine Tabelle umgeschrieben
#Das gleiche passiert auch beginnend mit der zweiten Zeile. Beides wird mit cbind in eine Tabelle  zusammengefügt. Dabei unterscheiden sich nur die Spalten mit den Allelen
#[,c(4,5)*-1] die 4. und die 5. Spalte werden gelöscht. Diese beiden Spalten enthielten nochmal Marker und Patient, genau wie die ersten beiden Spalten 
Profiles_cast <- cbind(Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt)), by=2),])[,c(4,5)*-1]
colnames(Profiles_cast) = c("Marker", "Patient", "Allel_1", "Allel_2")

#Run_01_short$trueAllel1 <- Profiles_cast$Allel_1[match(Run_01_short$Marker, Profiles_cast$Marker), match(Run_01_short$Patient, Profiles_cast$Patient)]
#Run_01_short$trueAllel1 <- Profiles_cast$Allel_1[match(match(Run_01_short$Marker, Profiles_cast$Marker), match(Run_01_short$Patient, Profiles_cast$Patient))]

# Vereint die Tabellen Run_01 und Profiles_cast miteinamder, sodass alle Spalten der Run_01 Tabelle erhalten bleiben und Allel1 und Allel2 dazu kommen
combined_tabel <- inner_join(Run_01_short, Profiles_cast)


#merge(Run_01_short, Profiles_cast, by.x = c("Patient", "Marker"), by.y = c("Patient", "Marker"), all.x = TRUE, all.y = TRUE)

#which( outer(combined_tabel$call_Allele, combined_tabel$Allel_1, "==") | outer(combined_tabel$call_Allele, combined_tabel$Allel_2, "=="), arr.ind=TRUE)

#neue Spalte für combined_tabel, diese Spalte heißt test. Wenn das call-Allel mit einem der korrekten Allele übereinstimmt, steht in der Spalte ein true, sonst ein wrong
for (i in 1:nrow(combined_tabel)){
if (combined_tabel$call_Allele[i] == combined_tabel$Allel_1[i] | combined_tabel$call_Allele[i] == combined_tabel$Allel_2[i]) {combined_tabel$test[i] = "true"} else {combined_tabel$test[i] = "wrong"}
}
return(combined_tabel)
}

