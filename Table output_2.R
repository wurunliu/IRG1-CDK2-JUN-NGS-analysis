mergedcounts<- read.table("mergedcounts M0 M3 C3 J3.txt", header = TRUE, quote = '\t',skip =0)
counts <- mergedcounts[,c(1,5,6,7,8,9,10)]
write.table(counts,"M3 C3 counts.txt", row.names = FALSE, quote = FALSE)



