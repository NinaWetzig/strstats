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


compare_call_to_Allele <- function(Run_01_file, Profiles_file, global_samples_file) {
  
  Run_01_file <- "/home/nina/Downloads/Run_01_short.out"
  Profiles_file <- "/home/nina/projects/strstats/resources/profiles.csv" 
  global_samples_file <- "/home/nina/projects/strstats/resources/global_sample_overview_forensic.csv"
  
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
  
  #Die Tabelle Profiles_melt hatte die Allele in einer Spalte untereinander, die neue Tabelle Profiles_cast hat für jedes Alell eine Spalte
  #Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt) Jede zweite Zeile der Tabelle beginnend mit der ersten Zeile wird in eine Tabelle umgeschrieben
  #Das gleiche passiert auch beginnend mit der zweiten Zeile. Beides wird mit cbind in eine Tabelle  zusammengefügt. Dabei unterscheiden sich nur die Spalten mit den Allelen
  #[,c(4,5)*-1] die 4. und die 5. Spalte werden gelöscht. Diese beiden Spalten enthielten nochmal Marker und Patient, genau wie die ersten beiden Spalten 
  Profiles_cast <- cbind(Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt)), by=2),])[,c(4,5)*-1]
  colnames(Profiles_cast) = c("Marker", "Patient", "Allel_1", "Allel_2")
  
  
  # Vereint die Tabellen Run_01 und Profiles_cast miteinamder, sodass alle Spalten der Run_01 Tabelle erhalten bleiben und Allel1 und Allel2 dazu kommen
  combined_table <- inner_join(Run_01_short, Profiles_cast)
  # Allel1 und Allel2 Spalten von character zu numeric
  combined_table$Allel_1 <- as.numeric(as.character(combined_table$Allel_1))
  combined_table$Allel_2 <- as.numeric(as.character(combined_table$Allel_2))

  #neue Spalte für combined_table, diese Spalte heißt Result. Wenn das call-Allel mit einem der korrekten Allele übereinstimmt, steht in der Spalte ein true, sonst ein wrong
  for (i in 1:nrow(combined_table)) {
    if (combined_table$call_Allele[i] == combined_table$Allel_1[i] | combined_table$call_Allele[i] == combined_table$Allel_2[i]) {combined_table$Result[i] = "SN"} 
    else if (combined_table$call_Allele[i] == combined_table$Allel_1[i] - 1 | combined_table$call_Allele[i] == combined_table$Allel_2[i] -1) {combined_table$Result[i] = "stutter"}
    else {combined_table$Result[i] = "LN"}
  }
  return(combined_table)
}


#Depth of Coverage
DoC_function <- function(combined_table) {
  DoC <- aggregate(combined_table$Reads, by = list(Run = combined_table$Run, Sample = combined_table$Sample, Marker = combined_table$Marker), FUN = function(x) {sum = sum(x)})
  colnames(DoC) = c("Run", "Sample", "Marker", "DoC")
  return(DoC)
}


#erstellt Tabelle, die call-Allele mit korrekten Allelen vergleicht und mit SN, LN, stutter oder true bennent
Allele_function <- function(combined_table) {
  #erstellt eine Tabelle, nur aus den Allelen mit SN. Ordnet die Tabelle nach Reads (die höchsten und zweithöchsten Reads)
  #diese Tabelle enthält nun alle richtigen true Allele 
  combined_table_true <- combined_table[combined_table$Result == "SN", ]
  combined_table_true <-  group_by(combined_table_true, Marker, Sample, call_Allele) %>%  mutate(rank = rank(desc(Reads))) %>% arrange(rank) %>% filter(rank <= 1)
  combined_table_true$Result <- gsub(pattern = "SN", replacement="true", x = combined_table_true$Result)
  combined_table_true <- select(combined_table_true, -rank)
  
  #ersetzt die SN in combined_table, die laut combined_table_true, eigentlich "true" sein müssten.
  Allel_table <- left_join(combined_table,combined_table_true, by=c("Run", "Sample", "Marker","call_Allele","Pattern", "Quality", "Reads", "Variants", "Patient","Allel_1","Allel_2")) %>% mutate(Result=coalesce(Result.y,Result.x)) %>% select(-Result.x,-Result.y)

  #löscht ] und Zahlen in der Tabelle, damit die Pattern verglichen werden können. 
  Allel_table$Pattern2 <- gsub(pattern = "[[:punct:]]", replacement="", x = Allel_table$Pattern) 
  Allel_table$Pattern2 <- gsub(pattern = "[0-9]", replacement="", x = Allel_table$Pattern2)

  #Sortiert die combined_table aufsteigend nach Run, Sample, Marker und absteigend nach Reads
  Allel_table <- Allel_table[order(Allel_table$Run, Allel_table$Sample, Allel_table$Marker, -Allel_table$Reads),]

  # Diese Tabelle enthält jetzt die vollständige Spalte Result mit Sn, LN, stutter und true
  SN_LN_stutter_true_table <- as.data.frame(do.call("rbind",(by(Allel_table, INDICES = list(Sample = Allel_table$Sample, Marker = Allel_table$Marker), function(x){
    # find stutter patterns
    spatterns = subset(x, Result == "stutter")$Pattern2
    for (spattern in spatterns){
    
    has_true_pattern = (subset(x, Pattern2 == spattern & Result == "true") %>% nrow) > 0
      if (!has_true_pattern){
        x$Result[x$Result == "stutter" & x$Pattern2 == spattern] = "SN"
      }
   }
    return(x)
   # find stutter pattern
   # check stutter and true pattern are equal
   # set stutter to SN
  }))))
  return(SN_LN_stutter_true_table)
}
