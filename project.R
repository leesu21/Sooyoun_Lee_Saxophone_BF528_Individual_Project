## PART 4
library(tidyverse)
library(data.table)

# data loading
dt <- fread("combat_adj.csv",drop = 1) %>% as.data.frame()

# filter
dt1 <- dt %>% mutate(row=fread("combat_adj.csv")$V1) %>%
  mutate(filter1=ifelse(rowSums(.>log(15,2))>.2*134,1,0)) %>%
  rowwise() %>% mutate(var=var(c_across(1:134)),
                       mean=mean(c_across(1:134))) %>% as.data.frame() %>%
  mutate(filter2=ifelse(var/median(var)>qchisq(0.99,133)/133,1,0)) %>%
  mutate(cv=sqrt(var)/mean,
         filter3=ifelse(cv>0.186,1,0))

# save the data for the biologist
dt2 <- dt1 %>% filter(filter1==1,filter2==1)
rownames(dt2) <- dt2$row
cat(c("Number of probes left:",dim(dt2)[1]))
write.csv(dt2[1:134], file = "part4_biologist.csv")

# save the data that pass all thresholds
dt3 <- dt1 %>% filter(filter1==1,filter2==1,filter3==1)
rownames(dt3) <- dt3$row
cat(c("Number of probes left:",dim(dt3)[1]))
write.csv(dt3[1:134], file = "part4.csv")


## PART 5

# Hierarchical Clustering
hc <- hclust(dist(t(dt3[1:134])))
plot(hc, labels = FALSE, main = "", sub="")
rect.hclust(hc,k=2)
table(cutree(hc,k=2))

# Heatmap
metadata <- read.csv("proj_metadata.csv")
colorbar <- ifelse(metadata$cit.coloncancermolecularsubtype=="C3","red","blue")
heatmap(as.matrix(dt3[1:134]), ColSideColors = colorbar, labRow = FALSE, labCol = FALSE)

# Welch t-test
grp <- as.numeric(cutree(hc,k=2))
welch_t <- data.frame()
for (i in 1:nrow(dt3)) {
  x <- dt3[i,which(grp==1)]
  y <- dt3[i,which(grp==2)]
  t <- t.test(x,y)
  welch_t <- rbind(welch_t,data.frame(id=rownames(dt3)[i],
                                      t=t$statistic,
                                      p=t$p.value))
}
welch_t$p_adj <- p.adjust(welch_t$p,method = "fdr")

# save the data
write.csv(welch_t, file = "part5.csv", row.names = FALSE)

# p_adj less than 0.05
cat(c("Number of probes left with p_adj<0.05:",nrow(welch_t[welch_t$p_adj<0.05,])))


#Same test for the biologist
hc1 <- hclust(dist(t(dt2[1:134])))

# Welch t-test
grp1 <- as.numeric(cutree(hc1,k=2))
welch_t1 <- data.frame()
for (i in 1:nrow(dt2)) {
  x <- dt2[i,which(grp1==1)]
  y <- dt2[i,which(grp1==2)]
  t <- t.test(x,y)
  welch_t1 <- rbind(welch_t1,data.frame(id=rownames(dt2)[i],
                                        t=t$statistic,
                                        p=t$p.value))
}
welch_t1$p_adj <- p.adjust(welch_t1$p,method = "fdr")

# save the data
write.csv(welch_t1, file = "part5_biologist.csv", row.names = FALSE)

# p_adj less than 0.05
cat(c("Number of probes left with p_adj<0.05(biologist):",nrow(welch_t1[welch_t1$p_adj<0.05,])))


## PART 6
if (!require("BiocManager")) install.packages("BiocManager")
if (!require("GSEABase")) install.packages("GSEABase")
if (!require("affy")) install.packages("affy")
if (!require("hgu133plus2.db")) install.packages("hgu133plus2.db")

# calling libraries
library(BiocManager)
library(hgu133plus2.db)
library(affy)
library(GSEABase)

biological_results <- read.csv("part5_biologist.csv", col.names = c('PROBEID', 't', 'p', 'p_adjust'))
genesymbol <- AnnotationDbi::select(hgu133plus2.db, biological_results$PROBEID, c('SYMBOL'))

#paste
genesymbol <- genesymbol %>%
  arrange(PROBEID,SYMBOL) %>%
  group_by(PROBEID) %>% summarise(SYMBOL=paste(SYMBOL,collapse="|"))

#merging the result of biological with genesymbol
merge_symbol <- merge(biological_results, genesymbol, on = 'PROBEID') %>%
  filter(!(is.na(SYMBOL) | SYMBOL == "" | SYMBOL == "NA"))

#to same genesymbol with its ids
countsignificant <- merge_symbol %>% group_by(SYMBOL) %>% count() %>% filter(n>=2) %>% arrange(-n)

#finding the most significant with common genesymbols
idiff <- data.frame(PROBEID = character(), t = numeric(), p = numeric(), p_adjust = numeric(), SYMBOL = character() )
for (i in countsignificant$SYMBOL) {
  x <- merge_symbol[merge_symbol$SYMBOL==i, ]
  x <- x[x$p_adjust == min(x$p_adjust), ]
  merge_symbol <- merge_symbol[!merge_symbol$SYMBOL == i,]
  idiff <- rbind(idiff, x)
}
merge_symbol <- rbind(merge_symbol, idiff)

#loading gene sets
hallmarks <- getGmt('h.all.v7.3.symbols.gmt')
GO <- getGmt('c5.go.v7.3.symbols.gmt')
KEGG <- getGmt('c2.cp.kegg.v7.3.symbols.gmt')

#geneset length
length(hallmarks)      #50
length(GO)             #10183
length(KEGG)           #186

# values in decreasing values of t-stats
merge_symbol <- merge_symbol %>% arrange(-t)
head(merge_symbol)
tail(merge_symbol)

#top 1000 up regulated & down regulated
up_1000 <- merge_symbol %>% head(1000)
down_1000 <- merge_symbol %>% tail(1000)

#selecting top 10
up_10 <- merge_symbol %>% head(10)
down_10 <- merge_symbol %>% tail(10)

#store the top10 up & down regulated
write.csv(down_10, "10_downregulated_genes.csv")
write.csv(up_10, "10_upregulated_genes.csv")

#genes that did not expressed
not_diffexp_up <- merge_symbol %>% filter(!merge_symbol$SYMBOL %in% up_1000$SYMBOL)
not_diffexp_down <- merge_symbol %>% filter(!merge_symbol$SYMBOL %in% down_1000$SYMBOL)

fishertest <- function(gl, gs, nde) {  #gl = genelist, gs= geneset, nde= not differentially expressed
  diffexp_ings <- length(intersect(gl,gs))
  diffexp_notgs <- length(gl) - diffexp_ings
  notde_ings <- length(intersect(nde,gs))
  notde_notgs <- length(nde) - notde_ings
  return(c(diffexp_ings,diffexp_notgs,notde_ings,notde_notgs))
}

#fisher test result of the hallmark
hallmarks_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(hallmarks)) {
  geneid <- geneIds(hallmarks[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  hallmarks_results[nrow(hallmarks_results)+1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  hallmarks_results[nrow(hallmarks_results)+1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')
}
hallmarks_results <- hallmarks_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

#fisher test result of the KEGG
kegg_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(KEGG)) {
  geneid <- geneIds(KEGG[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  kegg_results[nrow(kegg_results)+1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  kegg_results[nrow(kegg_results)+1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')
}
kegg_results <- kegg_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

#fisher test result of the kegg
go_results <- data.frame(setname = character(), pvalue = numeric(), estimate = numeric(), exp = character(), stringsAsFactors = FALSE)

for (i in 1:length(GO)) {
  geneid <- geneIds(GO[i])
  fisher_up <- fishertest(up_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_up$SYMBOL)
  fisher_down <- fishertest(down_1000$SYMBOL, geneid[[names(geneid)]], not_diffexp_down$SYMBOL)
  up <- fisher.test(matrix(fisher_up,nrow=2))
  down <- fisher.test(matrix(fisher_down, nrow=2))
  go_results[nrow(go_results)+1, ] <- c(names(geneid), up$p.value, up$estimate, 'UP')
  go_results[nrow(go_results)+1, ] <- c(names(geneid), down$p.value, down$estimate, 'Down')
}
go_results <- go_results %>% mutate(pvalue = as.numeric(pvalue), estimate = as.numeric(estimate))

go_results$BH <- p.adjust(go_results$pvalue, method = "BH", n = length(go_results$pvalue))
write.csv(go_results, "final_go.csv")

kegg_results$BH <- p.adjust(kegg_results$pvalue, method = "BH", n = length(kegg_results$pvalue))
write.csv(kegg_results, "final_kegg.csv")

hallmarks_results$BH <- p.adjust(hallmarks_results$pvalue, method = "BH", n = length(hallmarks_results$pvalue))
write.csv(hallmarks_results, "final_hallmarks.csv")
