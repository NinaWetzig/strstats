library(reshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(NMF)
library(gridExtra)
library(cowplot)

h_treshold <- 0.5

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
  raw_names <- c("Run", "Sample", "Marker", "Allele", "Pattern", "Reads")
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

combine_Allele <- function(Run_01_file, Profiles_file, global_samples_file, out_dir) {
  #Run_01_file <- "/home/nina/Downloads/Run_01_short.out"
  #Run_01_file <- "data/Run_07.out"
  #Profiles_file <- "resources/profiles.csv" 
  #global_samples_file <- "resources/global_sample_overview_forensic.csv"

  #einlesen der Run_01_short Tabelle
  Run_01_short <- read.csv(Run_01_file, na.strings=c("","NA"), sep =",", header = FALSE, stringsAsFactors = FALSE)
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
  
  ## NA mit "NA" ersetzen in Amel
  combined_table$Pattern <- apply(combined_table, 1 , FUN = function(x) {
    if (x[["Marker"]] == "Amelogenin") {
      return(ret = "NA") 
    } else {
        return(ret = x[["Pattern"]])
    }
  }) 
  
  combined_table <- aggregate(combined_table$Reads, by = list(Run = combined_table$Run, Sample = combined_table$Sample,Marker = combined_table$Marker,call_Allele = combined_table$call_Allele, Pattern = combined_table$Pattern, Patient = combined_table$Patient, Allel_1 = combined_table$Allel_1, Allel_2 = combined_table$Allel_2 ), FUN = function(x) sum = sum(x))
  colnames(combined_table)[9] <- "Reads"
  combined_table <- combined_table[,c(1, 2, 3, 4, 5, 9, 6, 7, 8)]
  
  #Variable, die angibt wo die Grafiken gespeichert werden
  write.csv(combined_table, file=paste(out_dir, "combined_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  return(combined_table)
}

compare_Allele <- function(combined_table, out_dir)  {
  combined_table <- group_by(combined_table, Marker, Sample, call_Allele) %>%  mutate(rank = rank(desc(Reads), ties.method = "last")) %>% arrange(rank)
  #im Fall von Amelogenin wird SN gesetzt
  combined_table$Result <- apply(combined_table, 1 , FUN = function(x) {
    if (!is.na(as.numeric(x[["call_Allele"]]))) {
      if (as.numeric(x[["call_Allele"]]) == as.numeric(x[["Allel_1"]]) || as.numeric(x[["call_Allele"]]) == as.numeric(x[["Allel_2"]])) {
        if (as.numeric(x[["Allel_1"]]) == as.numeric(x[["Allel_2"]])) {
          if (as.numeric(x[["rank"]]) <= 2) {
            return(ret = "true")
          } else {
            return(ret = "SN")
          }
        } else {
          if (as.numeric(x[["rank"]]) <= 1) {
            return(ret = "true")
          } else {
            return(ret = "SN")
          }
        }
      } else if (as.numeric(x[["call_Allele"]]) == as.numeric(x[["Allel_1"]])-1 || as.numeric(x[["call_Allele"]]) == as.numeric(x[["Allel_2"]])-1) {
        return(ret = "stutter")
      } else {
        return(ret = "LN")
      } 
    } else if (x[["call_Allele"]] == x[["Allel_1"]] || x[["call_Allele"]] == x[["Allel_2"]]) {
      if (x[["Allel_1"]] == x[["Allel_2"]]) {
        if (as.numeric(x[["rank"]]) <= 2) {
          return(ret = "true")
        } else {
          return(ret = "SN")
        }
      } else {
        if (as.numeric(x[["rank"]]) <= 1) {
          return(ret = "true")
        } else {
          return(ret = "SN")
        }
      }
    } else {
      return(ret = "SN")
    }
  })
  
  # In case of homocytoe samples, the first two alleles were set true. All true alleles with rank 2 are now compared to the corresponding rank 1 alleles and only remain true if there reads are at least - h_treshold - % of the reads of the rank 1 allele
  combined_table_true <- combined_table[combined_table$Result == "true",]
  
  combined_table_true$Result <- apply(combined_table_true, 1, FUN = function(x) {
    if (as.numeric(x[["rank"]]) == 2) {
      r1_reads <- filter(combined_table_true, Run == x[["Run"]], Sample == x[["Sample"]], Marker == x[["Marker"]], rank == as.numeric(1))$Reads
      if (!is.na(as.numeric(r1_reads))) {
        if (as.numeric(x[["Reads"]]) > (as.numeric(r1_reads)*as.numeric(h_treshold))) {
          return("true")
        } else {
          return("SN")
        }
      }
    } else {
      return("true")
    }
  })
  
  Allel_table <- left_join(combined_table,combined_table_true, by=c("Run", "Sample", "Marker","call_Allele","Pattern", "Reads", "Patient","Allel_1","Allel_2")) %>% mutate(Result=coalesce(Result.y,Result.x)) %>% select(-Result.x,-Result.y)
  Allel_table <- subset(Allel_table, select=-c(rank.x,rank.y))
  #löscht ] und Zahlen in der Tabelle, damit die Pattern verglichen werden können. 
  Allel_table$Pattern_short <- gsub(pattern = "[[:punct:]]", replacement="", x = Allel_table$Pattern) 
  Allel_table$Pattern_short <- gsub(pattern = "[0-9]", replacement="", x = Allel_table$Pattern_short)

  #Sortiert die combined_table aufsteigend nach Run, Sample, Marker und absteigend nach Reads
  Allel_table <- Allel_table[order(Allel_table$Run, Allel_table$Sample, Allel_table$Marker, -Allel_table$Reads),]

  # Diese Tabelle enthält jetzt die vollständige Spalte Result mit Sn, LN, stutter und true
  SN_LN_stutter_true_table <- as.data.frame(do.call("rbind",(by(Allel_table, INDICES = list(Sample = Allel_table$Sample, Marker = Allel_table$Marker), function(x){
  # find stutter patterns
  spatterns = subset(x, Result == "stutter")$Pattern_short
  for (spattern in spatterns){
    
    has_true_pattern = (subset(x, Pattern_short == spattern & Result == "true") %>% nrow) > 0
      if (!has_true_pattern){
        x$Result[x$Result == "stutter" & x$Pattern_short == spattern] = "SN"
      }
   }
    return(x)
   # find stutter pattern
   # check stutter and true pattern are equal
   # set stutter to SN
  }))))
  
  write.csv(SN_LN_stutter_true_table, file=paste(out_dir, "Allele_Result_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(SN_LN_stutter_true_table)
}

DoC_SCR_StR <- function(combined_table, SN_LN_stutter_true_table, out_dir) {
  
  #Depth of Coverage
  DoC <- aggregate(combined_table$Reads, by = list(Run = combined_table$Run, Marker = combined_table$Marker), FUN = function(x) {sum = sum(x)})
  colnames(DoC) = c("Run",  "Marker", "DoC")

  #diese Tabelle soll nur die Spalten enthalten, die für die SCR und DoC Berechnung benötigt werden
  SCR_Doc_table <- aggregate(SN_LN_stutter_true_table$Reads, by = list(Run = SN_LN_stutter_true_table$Run, Marker = SN_LN_stutter_true_table$Marker, Result = SN_LN_stutter_true_table$Result), FUN = function(x) sum = sum(x))
  SCR_Doc_table <- do.call(data.frame, SCR_Doc_table)
  colnames(SCR_Doc_table) <- c("Run", "Marker", "Result", "Reads")
  #Fügt Doc als Spalte hinzu
  SCR_Doc_table <- left_join(SCR_Doc_table,DoC)
  #neue Spalte mit SCR
  SCR_Doc_table$SCR <- SCR_Doc_table$Reads / SCR_Doc_table$DoC
  
  Stutter_Ratio <- filter(SCR_Doc_table, Result == "true" | Result == "stutter")
  runs <- levels(as.factor(SCR_Doc_table$Run))
  markers <- levels(as.factor(SCR_Doc_table$Marker))
  
  stutter_run_list <- data.frame()
  #Stutter_Ratio$St_Ratio <- ()
  
  for (run in runs) {
    temp_run <- filter(Stutter_Ratio, Run == run)
    for (marker in markers) {
      temp_marker <- filter(temp_run, Marker == marker)
      s <- temp_marker[temp_marker$Result == "stutter","Reads"]
      t <- temp_marker[temp_marker$Result == "true","Reads"]
      if (length(s) == 0 & length(t) == 0) {
        stutter_line <- NA
      }
      else if (length(s) == 0 & length(t) != 0) {
        stutter_line <- 0
      }
      else if (length(s) != 0 & length(t) == 0) {
        stutter_line <- "inf"
      }
      else if (length(s) != 0 & length(t) != 0) {
        stutter_line <- s/t
      }
      temp_line <- as.data.frame(cbind(run, marker, stutter_line))  
      stutter_run_list <- rbind(stutter_run_list, temp_line)
    }
    
  }
  colnames(stutter_run_list) <- c("Run", "Marker", "St_ratio")
  #Berechnung der Stutter Ratio
  #SCR_Doc_table$St_ratio <- SCR_Doc_table[SCR_Doc_table$Result == "stutter","Reads" ]/SCR_Doc_table[SCR_Doc_table$Result == "true","Reads" ]
  
  SN_Ratio <- filter(SCR_Doc_table, Result == "true" | Result == "stutter" | Result == "SN")
  
  SN_run_list <- data.frame()
  
  for (run in runs) {
    temp_run <- filter(SN_Ratio, Run == run)
    for (marker in markers) {
      temp_marker <- filter(temp_run, Marker == marker)
      s <- temp_marker[temp_marker$Result == "stutter","Reads"]
      t <- temp_marker[temp_marker$Result == "true","Reads"]
      sn <- temp_marker[temp_marker$Result == "SN","Reads"] 
      if (length(s) == 0 & length(t) == 0 & length(sn) == 0) {
        sn_line <- NA
      }
      else if (length(sn) != 0) {
        if (length(s) == 0 & length(t) == 0) {
          sn_line <- 1
        }
        else if (length(s) != 0 & length(t) == 0) {
          sn_line <- sn/(s+sn)
          }
        else if (length(s) == 0 & length(t) != 0) {
          sn_line <- sn/(t+sn)
        }
        else if (length(s) != 0 & length(t) != 0) {
          sn_line <- sn/(s+sn+t)
        }
      }
      else {
        sn_line <- 0
      }
      temp_line <- as.data.frame(cbind(run, marker, sn_line))  
      SN_run_list <- rbind(SN_run_list, temp_line)
    }
  }
  colnames(SN_run_list) <- c("Run", "Marker", "SN_ratio")
  
  
  #Berechnung SN Ratio
  #SCR_Doc_table$SN_ratio <- SCR_Doc_table[SCR_Doc_table$Result == "SN","Reads" ]/(SCR_Doc_table[SCR_Doc_table$Result == "stutter","Reads" ] + SCR_Doc_table[SCR_Doc_table$Result == "true","Reads" ] + SCR_Doc_table[SCR_Doc_table$Result == "SN","Reads" ])
  
  write.csv(SCR_Doc_table, file=paste(out_dir, "SCR_DoC_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  write.csv(stutter_run_list, file=paste(out_dir, "Stutter_ratio.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  write.csv(SN_run_list, file=paste(out_dir, "SN_ratio.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(SCR_Doc_table)
}  

ACR_function <- function(SN_LN_stutter_true_table, out_dir) {
  
  #Berechnung der Allele Coverage Ratio (ACR)
  Allele_Cov_Ratio <- aggregate(SN_LN_stutter_true_table$Reads, by = list(Run = SN_LN_stutter_true_table$Run, Marker = SN_LN_stutter_true_table$Marker, Sample = SN_LN_stutter_true_table$Sample, Result = SN_LN_stutter_true_table$Result, call_Allel = SN_LN_stutter_true_table$call_Allele), FUN = function(x) sum = sum(x))
  Allele_Cov_Ratio <- filter(Allele_Cov_Ratio, Result == "true")
  colnames(Allele_Cov_Ratio) <- c("Run", "Marker","Sample", "Result","call_Allele", "Reads")
  #Allele_Cov_Ratio <- Allele_Cov_Ratio %>% group_by(Run, Sample, Marker, Result, call_Allele) %>% summarise(Reads = sum(Reads))

  ACR_table <-do.call(rbind, 
                by(Allele_Cov_Ratio, 
                    INDICES = list(Run = Allele_Cov_Ratio$Run, Marker = Allele_Cov_Ratio$Marker, Sample = Allele_Cov_Ratio$Sample), 
                   FUN = function(x) {
                    # max_call_allel = max(as.numeric(x$call_Allele))
                     # max_call_allel_reads = x[x$call_Allele == max_call_allel, "Reads"]
                    
                     # min_call_allel = min(as.numeric(x$call_Allele))
                      #min_call_allel_reads = x[x$call_Allele == min_call_allel, "Reads"]
                     max_call_allel_reads = max(as.numeric(x$Reads))
                     min_call_allel_reads = min(as.numeric(x$Reads))
                     acr = (min_call_allel_reads/max_call_allel_reads) %>% unlist
                      x$ACR = acr
                      return(x)
                    })
               ) %>% as.data.frame
  result = list()
  
  write.csv(ACR_table, file=paste(out_dir, "ACR_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(ACR_table)
} 
  
create_plots <- function(ACR_table, combined_table, out_dir) {
  
  heatmap_table <- data.frame(combined_table$Run, combined_table$Sample, combined_table$Marker, combined_table$Reads)
  colnames(heatmap_table) <- c("Run", "Sample", "Marker", "Reads")   
  
  #Variable, die angibt wo die Grafiken gespeichert werden
  
  runs <- levels(heatmap_table$Run)
  
  #DoC Heatmap
  for (run in runs){
    # Für alle Run von 01 bis 10, also die, die in runs gespeichert sind
    data <- heatmap_table[heatmap_table$Run == run, ]
    # data wird alles genannt das in der Datei "input" unter Run mit "run" abgelegt war
    
    doc_x <- data[data$Marker == "Amelogenin", ]
    # hier wird gespeichert was aus "data" in der Marker-Spalte "Amelogenin" hat
    
    doc_x <- cast(doc_x, Sample ~ Marker, fun.aggregate = sum)
    # hier wird die doc_x Tabelle neu organisiert? (Marker gegen Sample)
    
    doc_x <- doc_x[,-1]
    # doc_x ohne die erste Spalte (ohne Sample)
    
    doc_a <- filter(data, Marker != "Amelogenin")
    # filtert in data nach Zeilen, die in der Marker-Spalte nicht Effective, Filtered oder Amelogenin haben.
    
    doc_a <- cast(doc_a, Sample ~ Marker, fun.aggregate = sum)
    # hier wird die doc_a Tabelle neu organisiert? (Marker gegen Sample)
    
    doc_m <- doc_a[,-1]
    # doc_a ohne die erste Spalte (ohne Sample)
    
    rownames(doc_m) <- doc_a[,1]
    # erste Spalte von doc_a wird als "rownames(doc_m)" gespeichert (Sample ist jetzt "rownames(doc_m)"?)
    
    doc_m <- data.matrix(doc_m)
    # erstellt Matrix aus doc_m
    
    doc_m <- doc_m[,order(colSums(doc_m), decreasing=TRUE)]
    # ordnet die Spalten nach Spaltensumme in absteigender Reihenfolge
    
    doc_m <- cbind(doc_m, doc_x)
    # verbindet doc_m und doc_x zu einer Matrix/Tabelle
    
    colnames(doc_m)[length(colnames(doc_m))] <- "Amelogenin" 
    # ?
    
    # Pfad anpassen
    write.table(doc_m, file = paste(out_dir, run, "_DoC.csv", sep = ""), append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", col.names = NA, qmethod = c("escape", "double"),fileEncoding = "")
    # Tabelle aus doc_m (also der doc_m die jetzt auch doc_x enthält)
    # eine Tabelle für jedes der run (Schleife)
    
    doc_m <- doc_m[order(rowSums(doc_m), decreasing=FALSE),]
    # ordner die Zeilen nach Zeilensummen in absteigender Reihenfolge
    
    # Pfad anpassen
    pdf(paste(out_dir, run, "_DoC.pdf", sep = ""), onefile = FALSE)
    aheatmap(doc_m, Rowv=NA, Colv=NA, color = "-RdYlBu2:100", main = run)
    # eine heatmap aus doc_m wird erstellt (Farben je nach Readzahl?)
    dev.off()
  }
  
  runs <- levels(ACR_table$Run)
  
  for (run in runs){
  
    data <- ACR_table[ACR_table$Run == run, ]  
    data$Marker <- as.factor(data$Marker)
    #ACR_table <- filter(ACR_table, Marker != "Amelogenin")
    
    #ACR Boxplots
    g <- ggplot(data, aes(x=reorder(Marker, ACR, FUN = median), y=ACR)) + geom_boxplot(color = "black", fill = "dodgerblue2")  
    g <- g + xlab("Marker") + geom_hline(yintercept=0.6, color =  "red")
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Allele Coverage Ratio per Marker") 
    g <- g + stat_summary(fun.y=mean, geom="point", color="white", fill="white") 
    g <- g + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.1))
    
    ggsave(file = paste(out_dir, run, "_ACR.pdf", sep = ""), plot = g, path = out_dir)
  }
}  

create_pies <- function(SN_LN_stutter_true_table, out_dir) {
  
  #Marker
  marker_table <- aggregate(SN_LN_stutter_true_table$Reads, by = list(Run = SN_LN_stutter_true_table$Run, Marker = SN_LN_stutter_true_table$Marker, Result = SN_LN_stutter_true_table$Result), FUN = function(x) sum = sum(x))
  marker_table <- filter(marker_table, Marker != "Amelogenin")
  markers <- levels(as.factor(marker_table$Marker))
  colnames(marker_table) <-  c("Run", "Marker", "Result", "Reads")

  #Average
  Average_table <- aggregate(SN_LN_stutter_true_table$Reads, by = list(Run = SN_LN_stutter_true_table$Run, Result = SN_LN_stutter_true_table$Result), FUN = function(x) sum = sum(x))
  Average_table$Marker <- "Average"
  Average_table <- Average_table[,c(1, 4, 2, 3)]
  colnames(Average_table) <-  c("Run", "Marker", "Result", "Reads")
  
  #SN:LN
  SNLN_table <- filter(marker_table, Result == "SN" | Result == "LN")
  SNLN_table$Marker <- "SN:LN"
  colnames(SNLN_table) <-  c("Run", "Marker", "Result", "Reads")
  
  pie_chart_table <- rbind(marker_table, Average_table, SNLN_table)
 
  #Pie charts zu allen Markern, Average und SN:LA 
  #pie_chart_table$Run <- as.factor(pie_chart_table$Run)
  #pie_chart_table$Marker <- as.factor(pie_chart_table$Marker)
  #pie_chart_table$Marker <- factor(pie_chart_table$Marker, levels=c("D10S1248","D12S391","D16S539","D18S51","D19S433", "D1S1656", "D21S11", "D22S1045","D2S1338","D2S441","D8S1179","FGA","SE33","TH01","vWA","Amelogenin","SN:LN","Average"))
  
  
  runs <- levels(as.factor(pie_chart_table$Run))
  markers <- c(markers, "SN:LN", "Average")
  
  for (run in runs) {
    
    data <- pie_chart_table[pie_chart_table$Run == run, ] 
    
    #pie_charts zur SCR
    #cols <- c("LN" = "red", "SN" = "brown4", "stutter" = "yellow", "true" = "blue")
    #pie  <- ggplot(data, aes(x="", y=Reads, fill=Result)) + 
     # geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
      #scale_fill_manual(values=cols)+
      #theme_minimal() +
      #theme(
       # axis.title.x = element_blank(),
       # axis.title.y = element_blank(),
       # panel.border = element_blank(),
        #panel.grid=element_blank(),
       # axis.ticks = element_blank(),
       # axis.text = element_blank()
      #) + facet_grid(~ Marker)
   # pie
    
    plots <- list()
    x <- 1
    for (marker in markers) {
      
      piedata <- data[data$Marker == marker, ]
      
      #pie_charts zur SCR
      cols <- c("LN" = "red", "SN" = "brown4", "stutter" = "yellow", "true" = "blue")
      pie  <- ggplot(piedata, aes(x="", y=Reads, fill=Result)) + 
        geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
        ggtitle(marker) + 
        scale_fill_manual(values=cols)+
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        ) +
        theme(legend.position="none")  
      
      plots[[x]] <- pie
      x <- x+1
      
    }
    pies <- do.call(grid.arrange,plots)
    
    legend <- get_legend(pie + theme(legend.position="bottom"))
    pie_chart <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
    
    ggsave("SCR.pdf", plot = pie_chart, path = out_dir) 
    write.csv(pie_chart_table, file=paste(out_dir, "plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
    
  }
  
  
  #StR plots
  StR_plot_table <- filter(SN_LN_stutter_true_table, Result == "true" | Result == "stutter")
  StR_plot_table <- aggregate(StR_plot_table$Reads, by = list(Run = StR_plot_table$Run, Marker = StR_plot_table$Marker, Result = StR_plot_table$Result), FUN = function(x) sum = sum(x))
  colnames(StR_plot_table) <-  c("Run", "Marker", "Result", "Reads")
  StR_plot_table <- filter(StR_plot_table, Marker != "Amelogenin")
  
  StR_plot_table$Run <- as.factor(StR_plot_table$Run)
  StR_plot_table$Marker <- as.factor(StR_plot_table$Marker)
  
  runs <- levels(StR_plot_table$Run)
  markers <- levels(StR_plot_table$Marker)
  
  for (run in runs) {
    
    data <- StR_plot_table[StR_plot_table$Run == run, ] 
    
    plots <- list()
    x <- 1
    for (marker in markers) {
      
      piedata <- data[data$Marker == marker, ]
      
      cols <- c("stutter" = "yellow", "true" = "blue")
      pie  <- ggplot(piedata, aes(x="", y=Reads, fill=Result)) + 
        geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
        ggtitle(marker) + 
        scale_fill_manual(values=cols)+
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        ) +
        theme(legend.position="none")  
      
      plots[[x]] <- pie
      x <- x+1
      print(pie)
      
    }
    pies <- do.call(grid.arrange,plots)
    
    legend <- get_legend(pie + theme(legend.position="bottom"))
    pies <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
    
    ggsave("StR.pdf", plot = pies, path = out_dir) 
    write.csv(StR_plot_table, file=paste(out_dir, "plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  }
 


#SNR plots
SNR_plot_table <- filter(SN_LN_stutter_true_table, Result == "true" | Result == "stutter" | Result == "SN")
SNR_plot_table <- aggregate(SNR_plot_table$Reads, by = list(Run = SNR_plot_table$Run, Marker = SNR_plot_table$Marker, Result = SNR_plot_table$Result), FUN = function(x) sum = sum(x))
colnames(SNR_plot_table) <-  c("Run", "Marker", "Result", "Reads")
SNR_plot_table <- filter(SNR_plot_table, Marker != "Amelogenin")

SNR_plot_table$Run <- as.factor(SNR_plot_table$Run)
SNR_plot_table$Marker <- as.factor(SNR_plot_table$Marker)

runs <- levels(SNR_plot_table$Run)
markers <- levels(SNR_plot_table$Marker)

  for (run in runs) {
    
    data <- SNR_plot_table[SNR_plot_table$Run == run, ] 
    
    plots <- list()
    x <- 1
    for (marker in markers) {
      
      piedata <- data[data$Marker == marker, ]
      
      cols <- c("stutter" = "blue", "true" = "blue", "SN" = "red")
      pie  <- ggplot(piedata, aes(x="", y=Reads, fill=Result)) + 
        geom_bar(width = 1, stat = "identity") + coord_polar("y", start=0) + 
        ggtitle(marker) + 
        scale_fill_manual(values=cols)+
        theme_minimal() +
        theme(
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_blank(),
          panel.grid=element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
        ) +
        theme(legend.position="none")  
      
      plots[[x]] <- pie
      x <- x+1
      print(pie)
      
    }
    pies <- do.call(grid.arrange,plots)
    
    legend <- get_legend(pie + theme(legend.position="bottom"))
    pies <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
    
    ggsave("SNR.pdf", plot = pies, path = out_dir) 
    write.csv(SNR_plot_table, file=paste(out_dir, "plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  }
  
}


