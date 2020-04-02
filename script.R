# charger les librairies
require(DESeq2)
# se placer dans le bon répertoire de travail
setwd("~/data")

# importer les données du plan expérimental
planExperimental <- read.table(file = "planExperimental.tsv", sep = "\t", header = TRUE, row.names = 1)
# mettre en forme les données du plan expérimental

# importer les données de mapping
resultatsMapping <- read.table(file = "StatisticsReportFINAL.resultsCounts.MappingResults.csv", header = TRUE)
# mettre en forme les données de mapping
  # supprimer la dernière ligne
resultatsMapping <- resultatsMapping[1:73013,]
  # ne conserver que la première colonne et les colonnes qui sont nommées *.map
resultatsMapping <- resultatsMapping[,c(1,grep("*.map",colnames(resultatsMapping)))]
  # nommer les lignes d'après le contenu de la première colonne
rownames(resultatsMapping) <- resultatsMapping$TranscriptID
  # supprimer la première colonne
resultatsMapping$TranscriptID <- NULL
  # renommer les colonnes en enlevant le "NG.7655_" en début et le ".map" en fin de libellé
colnames(resultatsMapping) <- gsub("NG.7655_","",colnames(resultatsMapping))
colnames(resultatsMapping) <- gsub("*.map","",colnames(resultatsMapping))
  # ranger les colonnes par ordre alphabétique de nom de colonne
resultatsMapping <- resultatsMapping[,order(colnames(resultatsMapping))]

# définition des variables utiles pour les fonctions de DESeq2
  # modèle d'analyse statistique
modele <- "~arrosage"
  # définition de ce que l'on veut comparer dans le modèle
contraste <- c("arrosage", "Normal","Stress")
  # définition du niveau de contrôle du taux de faux-positifs (pas plus de cinq pour mille faux-positifs)
fdrLevel <- 0.005  #####attention pas un p-value mais un probabilité pour un faux positive (regarder proportion des +ve et les -ve)

# structurer les données de mapping pour utilisation par DESeq2
countsData <- DESeqDataSetFromMatrix(countData = resultatsMapping, colData = planExperimental, design = formula(modele) )
# faire l'analyse statistique des données de mapping par DESes2 (modèle de distribution binomiale négative)
ddsFiltered <- DESeq(object = countsData, fitType = "mean", parallel = TRUE, quiet = TRUE)
# récupérer les résultats chiffrés de l'analyse statistique
resultsDifferentials <- results(object = ddsFiltered, alpha = fdrLevel, parallel = TRUE, contrast = contraste)  #contraste pour comparer control et stresser

# sous-tirer les transcrits différentiels sur-exprimés
up <- resultsDifferentials[resultsDifferentials$padj <= fdrLevel & resultsDifferentials$log2FoldChange>0 & ! is.na(resultsDifferentials$padj), ]
  # ordonner ces transcrits en commençant par le plus significatif (pvalue ajustée la plus faible)
up <- up[order(up$padj), ]
  # en faire un dataframe
up <- as.data.frame(up)

# idem pour les transcrits différentiellement sous-exprimés
down <- resultsDifferentials[resultsDifferentials$padj <= fdrLevel & resultsDifferentials$log2FoldChange<0 & ! is.na(resultsDifferentials$padj), ]
down <- down[order(down$padj), ]
down <- as.data.frame(down)

# idem pour les transcrits non différentiels
nonDiff <- resultsDifferentials[resultsDifferentials$padj > fdrLevel & ! is.na(resultsDifferentials$padj), ]
nonDiff <- nonDiff[order(nonDiff$padj),]
nonDiff <- as.data.frame(nonDiff)

                    baseMean log2FoldChange      lfcSE     stat        pvalue
Potri.012G036000.1 3401.0468       6.649614 0.19147868 34.72770 3.009561e-264

                            padj
Potri.012G036000.1 2.871637e-260

#base de donnée fitozomes


#basemean=moyenne entre stressé et normal
stats= stats de deSEQ
p-value= brutes
padj=p-value ajusté avec faux positive

p-value: la prob observer les donnée soumis au test (donnée experimantal) sous hyp null (hasard sous hyp null)

# réponse à la question : les 30 transcrits les plus sur-exprimés
head(x = up, n = 30)

# complément avec les représentations graphiques
# charger les librairies nécessaires
require(FactoMineR)
require(gplots)
require(RColorBrewer)
require(ggplot2)
require(corrplot)
require(wordcloud)

# calcul de l'ACP
  # transformer les données de base pour les rendre gaussiennes
rldResults <- rlogTransformation(object = ddsFiltered, blind = FALSE)

  # associer les données transformées et le plan expérimental
dataAcpRld <- cbind.data.frame(planExperimental, t( assay(rldResults) ) )

  # calculer l'ACP
resultsAcpRld <- PCA(dataAcpRld, ncp = nrow(dataAcpRld), scale.unit = TRUE, quali.sup = 1, graph = FALSE)

# classification des échantillons
  # calcul de la matrice de distance
distRld <- as.matrix( dist( t( assay(rldResults) ) ) )
  # représentation graphique par 'heatmap' du regroupement des échantillons
heatmap.2(x = distRld, trace = "none", dendrogram = "both", col = rev(brewer.pal(11, "PRGn")), key=FALSE, 
          labRow = gsub("NG-7655_", "", rownames(distRld)), labCol = planExperimental[rownames(distRld),"arrosage"], 
          cexRow = 1, cexCol = 1, main="", margins = c(10,10), 
          RowSideColors = brewer.pal(3,"Set3")[as.numeric(planExperimental$arrosage)], 
          ColSideColors = c("lightblue","moccasin")[as.numeric(planExperimental$arrosage)])

# placement des échantillons sur le plan principal de la PCA
  # tracé du cadre de fond
plot(x = resultsAcpRld$ind$coord[,1], y = resultsAcpRld$ind$coord[,2], type = "n", xlab = paste0("Dim 1 (",round(resultsAcpRld$eig[1,2],2)," %)"), ylab = paste0("Dim 2 (",round(resultsAcpRld$eig[2,2],2)," %)"), asp = 1)
  # réservation de l'espace de traçage
rect(xleft = par("usr")[1], ybottom = par("usr")[3], xright = par("usr")[2], ytop = par("usr")[4], col = "gray90")
  # grille de repérage sur le fond du graphique
grid(col = "white", lty = 1)
  # ligne sur l'axe des abscisses
abline(h = 0, lty = 1, col = "darkgrey")
  # ligne sur l'axe des ordonnées
abline(v = 0, lty = 1, col = "darkgrey")
  # tracé des échantillons sur le graphique
points(x = resultsAcpRld$ind$coord[,1], y = resultsAcpRld$ind$coord[,2],pch = 16, cex = 3, col = brewer.pal(3,"Set3")[as.numeric(planExperimental$arrosage)] )
  # ajout des noms d'échantillons
wordcloud::textplot(x = resultsAcpRld$ind$coord[,1], y = resultsAcpRld$ind$coord[,2], words = row.names(resultsAcpRld$ind$coord), new = FALSE, show.lines = FALSE, cex = 1)

# contribution des échantillons à la variance globale
corrplot(corr = resultsAcpRld$ind$cos2, is.corr = FALSE, tl.cex = .8, method = "square")

# placement des ellipes de confiance sur le plan principal de l'ACP
plotellipses(model = resultsAcpRld, level = 0.99, magnify = 1, cex = 0.5, pch.means = NA, pch = 20, keepvar = "quali.sup")
# tracé des transcrits sur le plan principal de l'ACP
  # cadre du graphique de la PCA
plot.PCA(x = resultsAcpRld, title = "PCA", xlim = c(-1,1), ylim = c(-1,1) )
  # transcrits non différentiels
points(x = resultsAcpRld$var$coord[rownames(nonDiff),1:2], cex = 0.4, pch = 16, col = "grey")
  # transcrits sur-exprimés
points(x = resultsAcpRld$var$coord[rownames(up),1:2], cex = 0.4, pch = 16, col = "green")
  # transcrits sous-exprimés
points(x = resultsAcpRld$var$coord[rownames(down),1:2], cex = 0.4, pch = 16, col = "magenta")
  # légende du graphique
legend(x = -2, y = -0.5, legend = c(paste0(nrow(up)," up-regulated"), paste0(nrow(down)," down-regulated"), paste0(nrow(nonDiff)," non-differential")), 
       text.col = c("green", "magenta", "grey"), cex = 0.8, pch = 16, col = c("green", "magenta", "grey") )

