names <- c("Marker", "HeLa", "Glenn", "Sarah" , "Trine" , "Nicole")

profiles <- read.csv("resources/profiles.csv", na.strings=c("","NA"), sep =",", col.names = names, header = TRUE, stringsAsFactors = FALSE)

names2 <- c("Marker", "Allele" , "Pattern" , "Pattern_Quality" , "Reads", "Variants")

J1_input <- read.csv("data/Run_01/variant_calling/J1_STR (copy).csv", na.strings=c("","NA"), sep = ",", col.names = names2, header = TRUE, stringsAsFactors = FALSE)
J2_input <- read.csv("data/Run_01/variant_calling/J2_STR.out", na.strings=c("","NA"), sep = ",", col.names = names2, header = TRUE, stringsAsFactors = FALSE)

# munging
profiles$Marker <- gsub(pattern = "\n", replacement="", x = profiles$Marker)
profiles$Marker <- gsub(pattern = " ", replacement="", x = profiles$Marker)


vWA <- J1_input[J1_input$Marker == "vWA", ]
TH01 <- J1_input[J1_input$Marker == "TH01", ]
c <- rbind(vWA, TH01)
head(subset(J1_input, Marker %in% c("vWA", "TH01")))
# mit subset kann die Tabelle leichter erstellt werden


#J1_input$Sample <- "J1"
#J2_input$Sample <- "J2"
J1_input <- cbind(Sample = "J1" ,J1_input )
J2_input <- cbind(Sample = "J2" ,J2_input )
#Beides fügt der Tabelle eine neue Spalte hinzu. Die ersten Befehle an letzer Stelle, die letzen Befehle als erste Spalte

J_input <- rbind(J1_input, J2_input) 


write.csv(J_input, file = "results/J_STR.csv", na = "", sep =",", quote = FALSE, row.names = FALSE)
