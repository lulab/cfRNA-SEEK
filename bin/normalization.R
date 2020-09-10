library("RUVSeq")
library("optparse")
option_list <- list( 
		    make_option("--counts", help="path of raw counts"),
		    make_option("--classes", help = "path of batch info"),
		    make_option("--tmm", help="TMM normalized expression"),
		    make_option("--ruv",help="RUVg normalized expression"),
		    make_option("--anova",help="output anova table")
	           )
opt <- parse_args(OptionParser(option_list=option_list))

cat("Load matrix and sample classes\n")
mat <- read.table(file=opt$counts,header = TRUE,sep="\t",check.names = F,row.names=1,stringsAsFactors = F)
labels <- read.csv(file=opt$classes,header=TRUE,sep="\t",stringsAsFactors = F,check.names = F)
sample_ids <- colnames(mat)
labels <- subset(labels,labels$sample_id %in% sample_ids)
sample_ids <- labels$sample_id
group <- as.factor(labels$label)
mat <- mat[,sample_ids]

cat("Identify stable genes\n")
design <- model.matrix(~group)
y <- DGEList(counts=mat, group=group)
y <- calcNormFactors(y, method="TMM")
y <- estimateDisp(y, design)
mat.tmm <- cpm(y)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2:6)
top <- topTags(lrt, n=nrow(mat))$table
n.mask <- as.integer(0.8*nrow(as.matrix(mat)))
empirical <- rownames(mat)[which(!(rownames(mat) %in% rownames(top)[1:n.mask]))]

cat("Remove batch effect\n")
res <- RUVg(x=log(mat.tmm+0.1), cIdx=empirical, k=2,isLog = T)
mat.tmm.ruv <- exp(res$normalizedCounts)

cat("Write to output\n")
write.table(mat.tmm,opt$tmm,quote=F,sep="\t")
write.table(mat.tmm.ruv,opt$ruv,quote=F,sep="\t")
write.table(top,opt$anova,quote=F,sep="\t")




