library(plyr)
library(pheatmap)


brain_protein = subset(brain_protein, select = -c(1,2))
#save.image(file = "Env.RData", version = NULL, ascii = FALSE, compress = FALSE, safe = TRUE)
a1 = cbind.data.frame(adrenal_gland_protein, appendix_protein, bone_marrow_protein, brain_protein, colon_protein, duodenum_protein, endometrium_uterine_protein, esophagus_protein, fallopian_tube_oviduct_protein, fat_adipose_tissue_protein, gallbladder_protein, heart_protein, kidney_protein, liver_protein, lung_protein, lymph_node_protein, ovary_protein, pancreas_protein, placenta_protein, prostate_protein, rectum_protein, salivary_gland_protein, small_intestine_protein, smooth_muscle_protein, spleen_protein, stomach_protein, testis_protein, thyroid_protein, tonsil_protein, urinary_bladder_protein)
a1 = a1[complete.cases(a1), ]
#save(a1, file = "matrix.RData")

#a2 = a1
a2 = merge(a2,a1,by="gene")

#b = c("DCAF1")
#b = as.data.frame(b)
#b2 = which(a1 == "DCAF1", arr.ind = TRUE)

a_DCAF7 = a1[13620,-c(1)]
row.names(a_DCAF7) = c("DCAF7")
a_DCAF13 = a1[4781,-c(1)]
row.names(a_DCAF13) = c("DCAF13")
a_DCAF8L1_L2 = a1[15762,-c(1)]
row.names(a_DCAF8L1_L2) = c("DCAF8L1/2")

a2 = rbind.data.frame(a_DCAF7, a_DCAF13, a_DCAF8L1_L2)
#a3 = a2
#a4 = rowMeans(a3)
#for (i in seq(1:length(rownames(a3)))) {
#  for (j in seq(1:length(colnames(a3)))) {
#    a3[i,j] = (a3[i,j] - a4[i])/a4[i]
#  }
#}
#a3[3,27] = 5

#4x10
pheatmap(a2, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), scale = "row" , cellwidth = 14, cellheight = 14, treeheight_row = 0, treeheight_col = 0, cluster_rows = TRUE, cluster_cols = TRUE)
#pheatmap(a3, treeheight_row = 0, treeheight_col = 0)


a4 = a1
a4$gene = as.character(a4$gene)
#save(a4, file = "matrix2.RData")
a4[11277,1] = "DCAF10"
a4[4,1]

a4[which(a4$gene == "DCAF7"), ]
#load("matrix2.RData")


colnames(protein_list)
b = protein_list$DCAFS
b = as.character(b)

d = NULL
for (i in seq(1:length(b))) {
  c = a4[which(a4$gene == b[i]), ]
  d = rbind.data.frame(d, c)
}

rownames(d) = NULL
d = ddply(d, "gene", numcolwise(max))
e = d$gene
f = d[,-1]
rownames(f) = e

g = rowSums(f)
g = as.data.frame(g)
colnames(g) = c("Intensity")
g = log10(g)

j = re$DCAFS
j = as.character(j)
f = f[j]
k = re2$DCAFS
k = as.character(k)
k = k[k != ""]
f = f[match(k, rownames(f)), ]


f1 = rowMeans(f)
for (i in seq(1:length(rownames(f)))) {
  for (j in seq(1:length(colnames(f)))) {
    if (f[i,j] == 0) {
      f[i,j] = NA
    } else {
      f[i,j] = (f[i,j] - f1[i])/f1[i]
    }
  }
}

#f[f == 0] <- NA

pheatmap(f, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), breaks = seq(-5, 5, length.out = 101), legend_breaks = c(-4, -2, 0, 2, 4), annotation_row = g, cellwidth = 14, cellheight = 14, treeheight_col = 0, treeheight_row = 0, cluster_cols = FALSE, cluster_rows = FALSE, na_col = "grey")
#pheatmap(f, color = colorRampPalette(c("navy", "white", "firebrick3"))(100), scale = "row", annotation_row = g, cellwidth = 14, cellheight = 14, treeheight_col = 0, treeheight_row = 0, cluster_cols = FALSE, cluster_rows = TRUE)

#DCAF
#CUL2
#CUL3_BBK
#CUL5
#FBOX
#NR
#DCAF_OR
#DCAF_B


#manual

#CUL2
#f = f[c("duodenum","small.intestine","brain","rectum","lung","adrenal.gland","liver","placenta","gallbladder","endometrium.uterine","kidney","fallopian.tube.oviduct","prostate","lymph.node","stomach","colon","urinary.bladder","testis","smooth.muscle","salivary.gland","tonsil","thyroid","pancreas","pancreas","heart","ovary","esophagus","spleen","appendix","fat.adipose.tissue","bone.marrow")]
