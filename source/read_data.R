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
  #Data_File <- "/home/jbierm2m/strstats/data/Run_11.out"
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
  
  runs <- levels(as.factor(combined_table$Run))
  for (run in runs) {
    run_dir <- paste(out_dir, run, "/", sep = "")
    dir.create(run_dir, recursive = TRUE)
    
    run_table <- combined_table[combined_table$Run == run, ]
    run_big_table <- compare_Allele(run_table, run_dir, run)
    
    write.csv(run_table, file=paste(run_dir, run, "_combined_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
    
    run_plots <- list()
    
    # DoC
    doc_returns <- doc(run_table, run_dir, run)
    doc_m <- doc_returns$doc_m
    doc_m_expanded <- doc_returns$doc_m_exp
    ### recording the plot in doc and returning the recorded plot does not work, the object type always changes from recordedplot to a simple list and then can't be replayed

    # ACR
    acr_boxplot <- acr(run_big_table, run_dir, run)
    run_plots <- list(run_plots, list(acr_boxplot))
    
    scr_pie_chart <- scr(run_table, run_big_table, run_dir, run)
    run_plots <- list(run_plots, list(scr_pie_chart))
    
    str_pie_chart <- str(run_big_table, run_dir, run)
    run_plots <- list(run_plots, list(str_pie_chart))
    
    snr_pie_chart <- snr(run_big_table, run_dir, run)
    run_plots <- list(run_plots, list(snr_pie_chart))
    
    pdf(paste(out_dir, run, "_Report.pdf", sep = ""), onefile = TRUE)
    aheatmap(doc_m, Rowv=NA, Colv=NA, color = "-RdYlBu2:100", main = "Depth of Coverage", sub = "Depth of Coverage is the amount of reads that have been assigned to a marker, per sample")
    for(my.plot in run_plots) {
      print(my.plot)
    }
    dev.off()
    
    filename <- paste(out_dir, run, "_Report.pdf", sep = "")
    try(system(paste("pdftk ", filename, " cat 2-end output ", filename, ".tmp && mv ", filename, ".tmp ", filename, sep = "")))
  }
  return(combined_table)
}

compare_Allele <- function(combined_table, out_dir, run)  {
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
  
  write.csv(SN_LN_stutter_true_table, file=paste(out_dir, run, "_Allele_Result_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(SN_LN_stutter_true_table)
}

### streamlined functions

doc <- function(combined_table, out_dir, run) {
  data <- data.frame(combined_table$Run, combined_table$Sample, combined_table$Marker, combined_table$Reads)
  colnames(data) <- c("Run", "Sample", "Marker", "Reads")
  
  doc_x <- data[data$Marker == "Amelogenin", ]
  doc_x <- cast(doc_x, Sample ~ Marker, fun.aggregate = sum)
  doc_x <- doc_x[,-1]
  
  doc_a <- filter(data, Marker != "Amelogenin")
  doc_a <- cast(doc_a, Sample ~ Marker, fun.aggregate = sum)
  
  doc_m <- doc_a[,-1]
  rownames(doc_m) <- doc_a[,1]
  doc_m <- data.matrix(doc_m)
  doc_m <- doc_m[,order(colSums(doc_m), decreasing=FALSE)]
  doc_m <- cbind(doc_m, doc_x)
  colnames(doc_m)[length(colnames(doc_m))] <- "Amelogenin"
  
  doc_m_expanded <- rbind(doc_m, colSums(doc_m))
  doc_m_expanded <- cbind(doc_m_expanded, rowSums(doc_m_expanded))
  colnames(doc_m_expanded)[length(colnames(doc_m_expanded))] <- "Global"
  rownames(doc_m_expanded)[length(rownames(doc_m_expanded))] <- "Global"
  write.table(doc_m_expanded, file = paste(out_dir, run, "_DoC_expanded.csv", sep = ""), append = FALSE, quote = TRUE, sep = ",", eol = "\n", na = "NA", dec = ".", col.names = NA, qmethod = c("escape", "double"),fileEncoding = "")
  
  doc_m <- doc_m[order(rowSums(doc_m), decreasing=TRUE),]
  
  pdf(paste(out_dir, run, "_DoC.pdf", sep = ""), onefile = FALSE)
  aheatmap(doc_m, Rowv=NA, Colv=NA, color = "-RdYlBu2:100", main = "Depth of Coverage", sub = "Depth of Coverage is the amount of reads that have been assigned to a marker, per sample")
  dev.off()
  returns <- list("doc_m" = doc_m, "doc_m_exp" = doc_m_expanded)
  return(returns)
}

acr <- function(sn_ln_stutter_true_table, out_dir, run) {
  #Berechnung der Allele Coverage Ratio (ACR)
  allele_cov_ratio <- aggregate(sn_ln_stutter_true_table$Reads, by = list(Run = sn_ln_stutter_true_table$Run, Marker = sn_ln_stutter_true_table$Marker, Sample = sn_ln_stutter_true_table$Sample, Result = sn_ln_stutter_true_table$Result, call_Allel = sn_ln_stutter_true_table$call_Allele), FUN = function(x) sum = sum(x))
  allele_cov_ratio <- filter(allele_cov_ratio, Result == "true")
  colnames(allele_cov_ratio) <- c("Run", "Marker","Sample", "Result","call_Allele", "Reads")
  #Allele_Cov_Ratio <- Allele_Cov_Ratio %>% group_by(Run, Sample, Marker, Result, call_Allele) %>% summarise(Reads = sum(Reads))
  
  ACR_table <-do.call(rbind, 
                      by(allele_cov_ratio, 
                         INDICES = list(Run = allele_cov_ratio$Run, Marker = allele_cov_ratio$Marker, Sample = allele_cov_ratio$Sample), 
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
  caption <- "Allele Coverage Ratio shows the ratio of lower allele coverage to higher allele coverage. The red line marks the STR threshold."
  g <- ggplot(ACR_table, aes(x=reorder(Marker, ACR, FUN = median), y=ACR)) + geom_boxplot(color = "black", fill = "dodgerblue2") + labs(caption = caption)
  g <- g + xlab("Marker") + geom_hline(yintercept=h_treshold, color =  "red")
  g <- g + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.caption = element_text(size = 12, hjust = 0.5)) + ggtitle("Allele Coverage Ratio per Marker") 
  g <- g + stat_summary(fun.y=mean, geom="point", color="white", fill="white") 
  g <- g + coord_cartesian(ylim=c(0, 1)) + scale_y_continuous(breaks=seq(0, 1, 0.1))
  
  ggsave(file = paste(run, "ACR.pdf", sep = "_"), plot = g, path = out_dir)
  return(g)
}

scr <- function(combined_table, sn_ln_stutter_true_table, out_dir, run) {
  data <- sn_ln_stutter_true_table
  #Depth of Coverage
  DoC <- aggregate(combined_table$Reads, by = list(Run = combined_table$Run, Marker = combined_table$Marker), FUN = function(x) {sum = sum(x)})
  colnames(DoC) = c("Run",  "Marker", "DoC")
  
  #diese Tabelle soll nur die Spalten enthalten, die für die SCR und DoC Berechnung benötigt werden
  scr_doc_table <- aggregate(data$Reads, by = list(Run = data$Run, Marker = data$Marker, Result = data$Result), FUN = function(x) sum = sum(x))
  scr_doc_table <- do.call(data.frame, scr_doc_table)
  colnames(scr_doc_table) <- c("Run", "Marker", "Result", "Reads")
  #Fügt Doc als Spalte hinzu
  scr_doc_table <- left_join(scr_doc_table,DoC)
  #neue Spalte mit SCR
  scr_doc_table$SCR <- scr_doc_table$Reads / scr_doc_table$DoC
  
  write.csv(scr_doc_table, file=paste(out_dir, run, "_SCR_DoC_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  ### create plot
  #Marker w/o Amel
  marker_table <- aggregate(data$Reads, by = list(Run = data$Run, Marker = data$Marker, Result = data$Result), FUN = function(x) sum = sum(x))
  marker_table <- filter(marker_table, Marker != "Amelogenin")
  colnames(marker_table) <-  c("Run", "Marker", "Result", "Reads")
  
  #marker levels
  markers <- levels(as.factor(marker_table$Marker))
  markers <- c(markers, "SN:LN", "Average")
  
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
  scr_plot_table <- rbind(marker_table, average_table, snln_table)
  
  plots <- list()
  x <- 1
  for (marker in markers) {
    piedata <- scr_plot_table[scr_plot_table$Marker == marker, ]
    
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
  
  caption = "Sequence Coverage Ratio is the ratio of true allele, N-1 stutter, SN, and LN to DoC."
  legend <- get_legend(pie + theme(legend.position="bottom"))
  pie_chart <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
  pie_chart <- ggdraw(add_sub(pie_chart, caption))
  title <- ggdraw() + draw_label("Sequence Coverage Ratio", fontface='bold')
  pie_chart <- plot_grid(title, pie_chart, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(paste(run, "SCR.pdf", sep = "_"), plot = pie_chart, path = out_dir) 
  write.csv(scr_plot_table, file=paste(out_dir, run, "_SCR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  return(pie_chart)
}

snr <- function(sn_ln_stutter_true_table, out_dir, run) {
  data <- sn_ln_stutter_true_table
  snr_ratio <- filter(data, Result == "true" | Result == "stutter" | Result == "SN")
  
  markers <- levels(as.factor(snr_ratio$Marker))
  snr_run_list <- data.frame()
  
  for (marker in markers) {
    temp_marker <- filter(snr_ratio, Marker == marker)
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
    snr_run_list <- rbind(snr_run_list, temp_line)
  }
  colnames(snr_run_list) <- c("Run", "Marker", "SN_ratio")
  
  write.csv(snr_run_list, file=paste(out_dir, run, "_SN_ratio.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  
  ### prepare plot
  snr_plot_table <- filter(data, Result == "true" | Result == "stutter" | Result == "SN")
  snr_plot_table <- aggregate(snr_plot_table$Reads, by = list(Run = snr_plot_table$Run, Marker = snr_plot_table$Marker, Result = snr_plot_table$Result), FUN = function(x) sum = sum(x))
  colnames(snr_plot_table) <-  c("Run", "Marker", "Result", "Reads")
  snr_plot_table <- filter(snr_plot_table, Marker != "Amelogenin")
  
  snr_plot_table$Run <- as.factor(snr_plot_table$Run)
  snr_plot_table$Marker <- as.factor(snr_plot_table$Marker)
  
  ### create plot
  plots <- list()
  x <- 1
  markers <- levels(as.factor(snr_plot_table$Marker))
  for (marker in markers) {
    
    piedata <- snr_plot_table[snr_plot_table$Marker == marker, ]
    
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
  
  caption = "Sequence Noise Ratio is the ratio of SN to (SN + Stutter + True)."
  legend <- get_legend(pie + theme(legend.position="bottom"))
  pie_chart <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
  pie_chart <- ggdraw(add_sub(pie_chart, caption))
  title <- ggdraw() + draw_label("Sequence Noise Ratio", fontface='bold')
  pie_chart <- plot_grid(title, pie_chart, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(paste(run, "SNR.pdf", sep = "_"), plot = pie_chart, path = out_dir) 
  write.csv(snr_plot_table, file=paste(out_dir, run, "_SNR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(pie_chart)
}

str <- function(sn_ln_stutter_true_table, out_dir, run) {
  data <- sn_ln_stutter_true_table
  stutter_ratio <- filter(data, Result == "true" | Result == "stutter")
  
  markers <- levels(as.factor(stutter_ratio$Marker))
  stutter_run_list <- data.frame()
  
  for (marker in markers) {
    temp_marker <- filter(stutter_ratio, Marker == marker)
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
  colnames(stutter_run_list) <- c("Run", "Marker", "St_ratio")
  
  write.csv(stutter_run_list, file=paste(out_dir, run, "_Stutter_Ratio.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  ### create plot
  str_plot_table <- filter(data, Result == "true" | Result == "stutter")
  str_plot_table <- aggregate(str_plot_table$Reads, by = list(Run = str_plot_table$Run, Marker = str_plot_table$Marker, Result = str_plot_table$Result), FUN = function(x) sum = sum(x))
  colnames(str_plot_table) <-  c("Run", "Marker", "Result", "Reads")
  str_plot_table <- filter(str_plot_table, Marker != "Amelogenin")
  
  str_plot_table$Run <- as.factor(str_plot_table$Run)
  str_plot_table$Marker <- as.factor(str_plot_table$Marker)
  
  plots <- list()
  x <- 1
  markers <- levels(as.factor(str_plot_table$Marker))
  for (marker in markers) {
    
    piedata <- str_plot_table[str_plot_table$Marker == marker, ]
    
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
  
  caption = "Stutter Ratio is the ratio of N-1 stutter to true allele."
  legend <- get_legend(pie + theme(legend.position="bottom"))
  pie_chart <- plot_grid(pies, legend, ncol = 1, rel_heights = c(1, .2))
  pie_chart <- ggdraw(add_sub(pie_chart, caption))
  title <- ggdraw() + draw_label("Stutter Ratio", fontface='bold')
  pie_chart <- plot_grid(title, pie_chart, ncol=1, rel_heights=c(0.1, 1))
  
  ggsave(paste(run, "StR.pdf", sep = "_"), plot = pie_chart, path = out_dir) 
  write.csv(str_plot_table, file=paste(out_dir, run, "_StR_plot_table.csv", sep=""),row.names=F,col.names=F, quote = FALSE)
  
  return(pie_chart)
}
