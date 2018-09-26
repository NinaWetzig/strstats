library(reshape)
library(tidyr)
library(dplyr)
library(ggplot2)
library(NMF)
library(gridExtra)
library(grid)
library(cowplot)

h_treshold <- 0.5

combine_Allele <- function(Data_File, Profiles_file, global_samples_file, out_dir) {
  #Data_File <- "/home/jbierm2m/strstats/data/Run_comb.out"
  #Profiles_file <- "resources/profiles.csv" 
  #global_samples_file <- "resources/global_sample_overview_forensic.csv"

  #einlesen der Run_01_short Tabelle
  input_data <- read.csv(Data_File, na.strings=c("","NA"), sep =",", header = FALSE, stringsAsFactors = FALSE)
  colnames(input_data) <- c("Run", "Sample", "Marker", "call_Allele", "Pattern", "Quality", "Reads", "Variants")
  
  #einlesen der profiles Tabelle
  Profiles <- read.csv(Profiles_file, na.strings=c("","NA"), sep =",", header = TRUE, stringsAsFactors = FALSE)
  Profiles_melt <- melt(Profiles, id = "Marker")
  colnames(Profiles_melt) = c("Marker", "Patient", "Allel")

  #einlesen der global_sample Tabelle
  global_samples <- read.csv(global_samples_file, na.strings=c("","NA"), sep =",", header = TRUE, stringsAsFactors = FALSE)
  
  #zu der Run_01_short Tabelle wird die Spalte Patienten zugefügt, nach der global_samples Tabelle
  input_data$Patient <- global_samples$Patient[match(input_data$Sample, global_samples$Sample)]
  
  #Die Tabelle Profiles_melt hatte die Allele in einer Spalte untereinander, die neue Tabelle Profiles_cast hat für jedes Alell eine Spalte
  #Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt) Jede zweite Zeile der Tabelle beginnend mit der ersten Zeile wird in eine Tabelle umgeschrieben
  #Das gleiche passiert auch beginnend mit der zweiten Zeile. Beides wird mit cbind in eine Tabelle  zusammengefügt. Dabei unterscheiden sich nur die Spalten mit den Allelen
  #[,c(4,5)*-1] die 4. und die 5. Spalte werden gelöscht. Diese beiden Spalten enthielten nochmal Marker und Patient, genau wie die ersten beiden Spalten 
  Profiles_cast <- cbind(Profiles_melt[seq(1, nrow(Profiles_melt), by=2),], Profiles_melt[seq(2, (nrow(Profiles_melt)), by=2),])[,c(4,5)*-1]
  colnames(Profiles_cast) = c("Marker", "Patient", "Allel_1", "Allel_2")
  
  # Vereint die Tabellen Run_01 und Profiles_cast miteinamder, sodass alle Spalten der Run_01 Tabelle erhalten bleiben und Allel1 und Allel2 dazu kommen
  combined_table <- inner_join(input_data, Profiles_cast)
  
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
  
  runs <- levels(as.factor(combined_table$Run))
  for (run in runs) {
    run_dir <- paste(out_dir, run, "/", sep = "")
    dir.create(run_dir, recursive = TRUE)
    
    run_table <- combined_table[combined_table$Run == run, ]
    run_big_table <- compare_Allele(run_table, run_dir)
    
    run_plots <- list()
    
    # DoC
    doc_heatmap <- doc(run_table, run_dir, run)
    run_plots <- list(run_plots, list(doc_heatmap))
    
    # ACR
    acr_boxplot <- acr(run_big_table, run_dir, run)
    run_plots <- list(run_plots, list(acr_boxplot))
    
  }
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

### to be remodeled

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
    
    
    
    data <- ACR_table[ACR_table$Run == run, ]  
    data$Marker <- as.factor(data$Marker)
    #ACR_table <- filter(ACR_table, Marker != "Amelogenin")
    
    #ACR Boxplots
    g <- ggplot(data, aes(x=reorder(Marker, ACR, FUN = median), y=ACR)) + geom_boxplot(color = "black", fill = "dodgerblue2")  
    g <- g + xlab("Marker") + geom_hline(yintercept=0.6, color =  "red")
    g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Allele Coverage Ratio per Marker") 
    g <- g + stat_summary(fun.y=mean, geom="point", color="white", fill="white") 
    g <- g + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.1))
    
    ggsave(file = paste(run, "ACR.pdf", sep = "_"), plot = g, path = out_dir)
  }
}

create_pies <- function(SN_LN_stutter_true_table, out_dir) {
  
  runs <- levels(as.factor(SN_LN_stutter_true_table$Run))
  for (run in runs) {
    data <- SN_LN_stutter_true_table[SN_LN_stutter_true_table$Run == run, ]
    
    #Marker w/o Amel
    marker_table <- aggregate(data$Reads, by = list(Run = data$Run, Marker = data$Marker, Result = data$Result), FUN = function(x) sum = sum(x))
    marker_table <- filter(marker_table, Marker != "Amelogenin")
    colnames(marker_table) <-  c("Run", "Marker", "Result", "Reads")
    
    #Average
    average_table <- aggregate(data$Reads, by = list(Run = data$Run, Result = data$Result), FUN = function(x) sum = sum(x))
    average_table$Marker <- "Average"
    average_table <- average_table[,c(1, 4, 2, 3)]
    colnames(average_table) <-  c("Run", "Marker", "Result", "Reads")
    
    #SN:LN
    snln_table <- filter(marker_table, Result == "SN" | Result == "LN")
    snln_table$Marker <- "SN:LN"
    colnames(snln_table) <-  c("Run", "Marker", "Result", "Reads")
    
    #SCR
    scr_table <- rbind(marker_table, average_table, snln_table)
    
    #StR
    str_table <- filter(data, Result == "true" | Result == "stutter")
    str_table <- aggregate(str_table$Reads, by = list(Run = str_table$Run, Marker = str_table$Marker, Result = str_table$Result), FUN = function(x) sum = sum(x))
    colnames(str_table) <-  c("Run", "Marker", "Result", "Reads")
    str_table <- filter(str_table, Marker != "Amelogenin")
    
    str_table$Run <- as.factor(str_table$Run)
    str_table$Marker <- as.factor(str_table$Marker)
    
    #SNR
    snr_table <- filter(data, Result == "true" | Result == "stutter" | Result == "SN")
    snr_table <- aggregate(snr_table$Reads, by = list(Run = snr_table$Run, Marker = snr_table$Marker, Result = snr_table$Result), FUN = function(x) sum = sum(x))
    colnames(snr_table) <-  c("Run", "Marker", "Result", "Reads")
    snr_table <- filter(snr_table, Marker != "Amelogenin")
    
    snr_table$Run <- as.factor(snr_table$Run)
    snr_table$Marker <- as.factor(snr_table$Marker)
    
    #marker levels
    markers <- levels(as.factor(marker_table$Marker))
    markers <- c(markers, "SN:LN", "Average")
    
    # create plots
    # scr
    plots <- list()
    x <- 1
    for (marker in markers) {
      piedata <- scr_table[scr_table$Marker == marker, ]
      
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
    
    ggsave(paste(run, "SCR.pdf", sep = "_"), plot = pies, path = out_dir) 
    write.csv(scr_table, file=paste(out_dir, run, "_SCR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
    # str
    plots <- list()
    x <- 1
    markers <- levels(as.factor(str_table$Marker))
    for (marker in markers) {
      
      piedata <- str_table[str_table$Marker == marker, ]
      
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
      
    }
    pies <- do.call(grid.arrange,plots)
    
    legend <- get_legend(pie + theme(legend.position="bottom"))
    pies <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
    
    
    ggsave(paste(run, "StR.pdf", sep = "_"), plot = pies, path = out_dir) 
    write.csv(str_table, file=paste(out_dir, run, "_StR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
    
    # sn
    plots <- list()
    x <- 1
    markers <- levels(as.factor(snr_table$Marker))
    for (marker in markers) {
      
      piedata <- snr_table[snr_table$Marker == marker, ]
      
      cols <- c("stutter" = "darkblue", "true" = "blue", "SN" = "red")
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
    pies <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
    
    ggsave(paste(run, "SNR.pdf", sep = "_"), plot = pies, path = out_dir) 
    write.csv(snr_table, file=paste(out_dir, run, "_SNR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  }
}


### streamlined functions

doc <- function(combined_table, out_dir, run) {
  data <- data.frame(combined_table$Run, combined_table$Sample, combined_table$Marker, combined_table$Reads)
  colnames(data) <- c("Run", "Sample", "Marker", "Reads")   
  
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
  
  pdf(paste(out_dir, run, "_DoC.pdf", sep = ""), onefile = FALSE)
  aheatmap(doc_m, Rowv=NA, Colv=NA, color = "-RdYlBu2:100", main = run)
  doc_heatmap <- recordPlot()
  # eine heatmap aus doc_m wird erstellt (Farben je nach Readzahl?)
  dev.off()
  return(doc_heatmap)
}

acr <- function(sn_ln_stutter_true_table, out_dir, run) {
  #Berechnung der Allele Coverage Ratio (ACR)
  Allele_Cov_Ratio <- aggregate(sn_ln_stutter_true_table$Reads, by = list(Run = sn_ln_stutter_true_table$Run, Marker = sn_ln_stutter_true_table$Marker, Sample = sn_ln_stutter_true_table$Sample, Result = sn_ln_stutter_true_table$Result, call_Allel = sn_ln_stutter_true_table$call_Allele), FUN = function(x) sum = sum(x))
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
  
  write.csv(ACR_table, file=paste(out_dir, run, "_ACR_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  ACR_table$Marker <- as.factor(ACR_table$Marker)
  
  #ACR Boxplots
  g <- ggplot(ACR_table, aes(x=reorder(Marker, ACR, FUN = median), y=ACR)) + geom_boxplot(color = "black", fill = "dodgerblue2")  
  g <- g + xlab("Marker") + geom_hline(yintercept=0.6, color =  "red")
  g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle("Allele Coverage Ratio per Marker") 
  g <- g + stat_summary(fun.y=mean, geom="point", color="white", fill="white") 
  g <- g + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.1))
  
  ggsave(file = paste(run, "ACR.pdf", sep = "_"), plot = g, path = out_dir)
  return(g)
}

scr <- function() {
  
}

snr <- function() {
  
}

str <- function() {
  
}
