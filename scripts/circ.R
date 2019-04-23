library(circlize)
library(randomcoloR)
library(ggplot2)

df <- read.table("../meta/ligand_receptor_villi.tsv", sep="\t", header=T, stringsAsFactors = F)
allnames <- c(df$from, df$to)
shortnames <- sapply(strsplit(allnames, "_"), "[[", 2)

pdf("../plots/circ_receptor_ligand.pdf", width=10, height=10) 

chordDiagram(df[,1:3], transparency = 1, annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = 0.1),
             directional = 1, direction.type = "arrows",
             link.arr.col = "black", link.arr.length = 0.2, link.arr.lwd = df$width,
             order=allnames[order(allnames)])

ec <- grep("EC", allnames, value=T)
highlight.sector(ec, track.index = 1, col = "red", 
                 text = "EC", cex = 0.8, text.col = "white", niceFacing = TRUE)
fb <- grep("FB", allnames, value=T)
highlight.sector(fb, track.index = 1, col = "blue", 
                 text = "FB", cex = 0.8, text.col = "white", niceFacing = TRUE)
hb <- grep("HB", allnames, value=T)
highlight.sector(hb, track.index = 1, col = "yellow", 
                 text = "HB", cex = 0.8, text.col = "black", niceFacing = TRUE)
mf <- grep("MF", allnames, value=T)
highlight.sector(mf, track.index = 1, col = "green", 
                 text = "MF", cex = 0.8, text.col = "black", niceFacing = TRUE)
vct <- grep("VCT", allnames, value=T)
highlight.sector(vct, track.index = 1, col = "pink", 
                 text = "VCT", cex = 0.8, text.col = "black", niceFacing = TRUE)
evt <- grep("EVT", allnames, value=T)
highlight.sector(evt, track.index = 1, col = "orange", 
                 text = "EVT", cex = 0.8, text.col = "white", niceFacing = TRUE)
sct <- grep("SCT", allnames, value=T)
highlight.sector(sct, track.index = 1, col = "cyan", 
                 text = "SCT", cex = 0.8, text.col = "black", niceFacing = TRUE)


shortnames <- unique(shortnames)
palette <- distinctColorPalette(length(shortnames))
for (i in 1:length(shortnames)) {
  name <- grep(shortnames[i], allnames, value=T)
  highlight.sector(name, track.index = 2, col = palette[i], 
                 text = shortnames[i], cex = 0.6, text.col = "black", niceFacing = TRUE)
}

dev.off()
