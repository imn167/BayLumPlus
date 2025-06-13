library(ArchaeoPhases)


#reinitialize memory
rm(list=ls())

setwd("D:/Users/Pc/Desktop")

Path=c("R/DataManipulation/")

########### sans micocouliers et sans UCIAMS-109995 (reconnu par BayLum comme un outlier)


C14_SampleNames = c("OxA-21261","OxA-21262","UCIAMS-98210","UCIAMS-103134","P-782","UCIAMS-103138","OxA-23247",
                    "OxA-9774","OxA-9946","UCIAMS-98208","UCIAMS-98209","OxA-9947","OxA-9775","OxA-9948","OxA-23523",
                    "UCIAMS-103135","UCIAMS_103136","UCIAMS-103137",'OxA-23248',"UCIAMS-103141","UCIAMS-103140",
                    "OxA-27087","OxA-9949","UCIAMS-103139","OxA-23249","OxA-9950","UCIAMS-109991",
                    "OxA-23250","UCIAMS-109992","OxA-9776","OxA-23251","UCIAMS-109993","OxA-9892",
                    "OxA-9777","UCIAMS-109994","OxA-9778","OxA-23252",
                    "OxA-9893","PL-980525A","AA-27982")

C14_Nb_sample = length(C14_SampleNames)

C14ages = c(8033,7955,7940,7955,8092,7920,8027,7935,7980,7965,7990,7985,8090,8090,7931,7970,7940,7965,8082,7980,
            8025,8000,8050,7970,8024,8030,8035,8085,8030,7985,8137,8160,8150,8160,8210,8240,
            8199,8155,8390,8195)

C14agesEr = c(39,40,30,25,98,25,37,50,55,25,25,50,55,50,38,25,25,30,37,25,25,50,40,25,35,50,30,36,30,55,36,
              30,50,50,30,55,36,50,90,80)

AC14_WithStratiWithout109995= AgeC14_Computation(Data_C14Cal=C14ages, Data_SigmaC14Cal=C14agesEr,
                                                SampleNames=C14_SampleNames, Nb_sample = C14_Nb_sample,
                                                PriorAge = rep(c(7, 13), C14_Nb_sample), SavePdf
                                                = TRUE,
                                                OutputFileName = c("MCMCplot",
                                                                   "HPD_CalC-14Curve", "summary"),
                                                OutputFilePath = Path, SaveEstimates = TRUE,
                                                OutputTableName = c("AllC14"), OutputTablePath
                                                = c(""),
                                                StratiConstraints = c(), sepSC = c(","), Model       #contraintes stratigraphiques
                                                = c("full"),
                                                CalibrationCurve = c("IntCal20"), Iter = 5000,
                                                t = 5,
                                                n.chains = 3, quiet = FALSE)

edges <- cbind(1:(C14_Nb_sample-1), 2:C14_Nb_sample)
G = igraph::graph_from_edgelist(edges)

plot(G,
     vertex.label = V(G)$name,
     vertex.color = "lightblue",
     vertex.size = 5,
     edge.arrow.size = 0.5,
     layout = layout_with_sugiyama(G))

AC14_WithStratiWithout109995$Sampling[[1]]

PlotIsotonicCurve(G, AC14_WithStratiWithout109995, level = .68)



