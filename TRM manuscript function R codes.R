#Fig.1C 1D 1F & Fig.S3 
REP_total <- function(data, pre_donor, pre_reci){
  l <- length(data)
  Rep <- list()
  Rep[1] <- as.data.frame(pre_donor)
  Rep[2] <- as.data.frame(pre_reci)
  add <- c(Rep[[1]], Rep[[2]])
  for (i in 1:l) {
    Rep[i+2] <- as.data.frame(setdiff(rownames(data[data[,i]>0,]), add))
    add <- c(add, Rep[[i+2]])
  }
  matrix <- matrix(nrow = l, ncol = (l+2))  
  rownames(matrix) <- colnames(data)
  colnames(matrix) <- c(rep(c("pre donor", "pre reci", colnames(data)), 1))
  for (m in 1:(l+2)){
    for (n in 1:l) {
      matrix[n,m] <- sum(data[intersect(rownames(data[data[,n]>0,]), Rep[[m]]),n])/sum(data[data[,n]>0,n])*100
    }
  }
  return(matrix)
}

#Fig.1E 1G
JSD_total <- function(data){
  l <- length(data)
  matrix <- matrix(nrow = 1, ncol = l-1)  
  colnames(matrix) <- c(colnames(data)[-1])
  for (i in 1:(l-1)) {
    matrix[1,i] <- jensen_shannon(data[,i],data[,(i+1)])
  }
  return(matrix)
}

jensen_shannon <- function(p, q){
  p = p / sum(p)
  q = q / sum(q)
  Hj = shannon.entropy(0.5 *(p+q)) 
  Hp = shannon.entropy(p) 
  Hq = shannon.entropy(q)
  jsd = Hj - 0.5*(Hp+Hq)
  return(jsd)
}

shannon.entropy <- function(p){
  if (min(p) < 0 || sum(p) <= 0)
    return(NA)
  p.norm <- p[p>0]/sum(p)
  -sum(log2(p.norm)*p.norm)
}

#Fig.2A Fig.5A 5B Fig.S9A
Share <- function(data, gut, MLN, donor){
  matrix <- as.data.frame(matrix(nrow = 4, ncol = ncol(data)))
  rownames(matrix) <- c("gutonlyFre%", "MLNonlyFre%", "sharedFre%", "unmappableFre%")
  colnames(matrix) <- colnames(data)
  
  error <- intersect(donor, union(gut, MLN))
  gutonly <- setdiff(gut, union(MLN, donor))
  MLNonly <- setdiff(MLN, union(gut, donor))
  share <- setdiff(intersect(gut, MLN),donor)
  donor <- setdiff(donor, union(gut, MLN))
  data <- data[!(rownames(data) %in% c(error, donor)),]
  unmappable <- setdiff(rownames(data), c(gut, MLN, donor))
  data1 <- normalize(data)
  
  for (i in 1:ncol(data)){
    matrix[1,i] <- sum(data1[gutonly,i])*100
    matrix[2,i] <- sum(data1[MLNonly,i])*100
    matrix[3,i] <- sum(data1[share,i])*100
    matrix[4,i] <- sum(data1[unmappable,i])*100
  }
  return(matrix)
}

#Fig.2B
Share_v1 <- function(data, gut, MLN_spleen, donor){
  l <- ncol(data)
  matrix <- matrix(nrow = l, ncol = 12)  
  colnames(matrix) <- c( "gutonly Pre-transplant clones detected", "gutonly post-Tx clones undetected", 
                         "MLN_spleenonly Pre-transplant clones detected", "MLN_spleenonly post-Tx clones undetected",
                         "shared Pre-transplant clones detected", "shared post-Tx clones undetected", 
                         "gutonly MLN_spleenonly Odds Ratio", "gutonly MLN_spleenonly P-value", 
                         "gutonly shared Odds Ratio", "gutonly shared P-value",
                         "MLN_spleenonly shared Odds Ratio", "MLN_spleenonly shared P-value")
  rownames(matrix) <- colnames(data)
  
  error <- intersect(donor, union(gut, MLN_spleen))
  gutonly <- setdiff(gut, union(MLN_spleen, donor))
  MLN_spleenonly <- setdiff(MLN_spleen, union(gut, donor))
  share <- setdiff(intersect(gut, MLN_spleen),donor)
  donor <- setdiff(donor, union(gut, MLN_spleen))
  data <- data[!(rownames(data) %in% c(error, donor)),]
  unmappable <- setdiff(rownames(data), c(gut, MLN_spleen, donor))
  
  for (i in 1:l) {
    matrix[i,1] <- length(gutonly)
    matrix[i,2] <- length1(data[gutonly,i])
    matrix[i,3] <- length(MLN_spleenonly)
    matrix[i,4] <- length1(data[MLN_spleenonly,i])
    matrix[i,5] <- length(share)
    matrix[i,6] <- length1(data[share,i])
    a <- fisher.test(matrix(c(matrix[i,2], matrix[i,1], matrix[i,4], matrix[i,3]),nrow=2))
    matrix[i,7] <- a[[3]]
    matrix[i,8] <- a[[1]]
    b <- fisher.test(matrix(c(matrix[i,2], matrix[i,1], matrix[i,6], matrix[i,5]),nrow=2))
    matrix[i,9] <- b[[3]]
    matrix[i,10] <- b[[1]]
    c <- fisher.test(matrix(c(matrix[i,4], matrix[i,3], matrix[i,6], matrix[i,5]),nrow=2))
    matrix[i,11] <- c[[3]]
    matrix[i,12] <- c[[1]]
  }
  return(matrix)
}

#Fig.2C 2D Fig.3B
matrix <- as.data.frame(matrix(nrow=length(data) , ncol=length(data) ))
colnames(matrix) <- colnames(data)
rownames(matrix) <- colnames(data)
for (m in 1:ncol(matrix)){
  for (n in 1:nrow(matrix)) {
    matrix[n,m] <- sum(data[intersect(rownames(data)[data[,n]>0], share),m])/sum(data[share,m])*100
  }
}

#Fig.3A
library(scales)
library(ggplot2)

data <- normalize(data)
data <- data[which(rowSums(data) > 0),]
data$Type <- NA
data[(data$Blood>0) & (data$Tissue == 0), c("Type")] <- c("Blood_Only")
data[(data$Blood == 0) & (data$Tissue>0), c("Type")] <- c("ileum_Only")
data[(data$Blood>0) & (data$Tissue>0), c("Type")] <- c("Shared")
SharedinGut <- sum(data[data$Type == c("Shared"),c("Tissue")])/sum(data[,c("Tissue")])*100
SharedinBlood <- sum(data[data$Type == c("Shared"),c("Blood")])/sum(data[,c("Blood")])*100
data[data$Blood == 0, c("Blood")] <- 0.000001
data[data$Tissue == 0, c("Tissue")] <- 0.000001

color <- c("#F08080","#90EE90","#ADD8E6")

dfPlot = ggplot(data, aes(x = Tissue, y = Blood, fill = Type)) +
  geom_point(size=6.5, shape=21, show.legend = F, color="gray20") + 
  scale_colour_manual(values = color, aesthetics = "fill")+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = 
                       element_blank(),panel.grid.minor = element_blank()) +theme(axis.line = 
                                                                                    element_line(colour = "black"), axis.title = element_text(size=20, family="Arial",
                                                                                                                                              color="black", face = "bold", vjust=0.5, hjust=0.5), axis.text = 
                                                                                    element_text(size = 20, family="Arial", color = "black", face="bold", vjust=0.5,hjust=0.5)) + 
  scale_y_log10(limits =  c(0.000001,1),breaks=trans_breaks("log10", function(x)10^x),labels =   trans_format("log10", math_format(10^.x))) + 
  scale_x_log10(limits =  c(0.000001,1),breaks=trans_breaks("log10", function(x)10^x),labels =   trans_format("log10", math_format(10^.x))) +
  xlab("gut") + ylab("blood") 

print(dfPlot)

#Fig.3C 3D Fig.4D 4E
cosine <- function(p, q){
  l <- length(p)
  a <- 0
  b <- 0
  c <- 0
  for (i in 1:l) {
    a <- a + p[i]*q[i]
    b <- b + p[i]*p[i]
    c <- c + q[i]*q[i]
  }
  d <- a / (sqrt(b)*sqrt(c))
  return(d)
}

#Fig.4B 4C Fig.S6
data <- data[,c(ileum, colon,native_colon)]

ileum_only <- rownames(data[data$ileum >0 & (data$colon == 0 & data$native_colon == 0), ])
colon_only <- rownames(data[data$colon >0 & (data$ileum == 0 & data$native_colon == 0), ])
native_only <- rownames(data[data$native_colon >0 & (data$ileum == 0 & data$colon == 0), ])
ileum_colon_share <- rownames(data[(data$ileum >0 & data$colon >0) & data$native_colon == 0,])
ileum_native_colon_share <- rownames(data[(data$ileum >0 & data$native_colon >0) & data$colon == 0,])
colon_native_colon_share <- rownames(data[(data$colon >0 & data$native_colon >0) & data$ileum == 0,])
triple_share <- rownames(data[(data$colon >0 & data$native_colon >0) & data$ileum >0,])

matrix <- as.data.frame(matrix(nrow=1, ncol=9))
colnames(matrix) <- c("ileum_Non_shared","ileum_Double_shared","ileum_Triple_shared", 
                      "colon_Non_shared","colon_Double_shared","colon_Triple_shared",
                      "native_colon_Non_shared","native_colon_Double_shared","native_colon_Triple_shared")
matrix[1,1] <- sum(data[ileum_only, 1])
matrix[1,2] <- sum(data[c(ileum_colon_share, ileum_native_colon_share), 1])
matrix[1,3] <- sum(data[triple_share, 1])
matrix[1,4] <- sum(data[colon_only, 2])
matrix[1,5] <- sum(data[c(ileum_colon_share, colon_native_colon_share), 2])
matrix[1,6] <- sum(data[triple_share, 2])
matrix[1,7] <- sum(data[native_only, 3])
matrix[1,8] <- sum(data[c(colon_native_colon_share, ileum_native_colon_share), 3])
matrix[1,9] <- sum(data[triple_share, 3])

#Fig.S7 Fig.S8 & HvG/NonHvG define
rcd4 <-  data[,c("R4U","R4L")]
rcd4 <-  normalize(rcd4)
rCD4HVG <-  rownames(rcd4[rcd4[,2]>0.00002 & rcd4[,2] > rcd4[,1]*2,])
rCD4NonHVG <- setdiff(rownames(data[data[,7]>0|data[,9]>0,]), rCD4HVG)
rcd8 <-  data[,c("R8U","R8L")]
rcd8 <-  normalize(rcd8)
rCD8HVG <-  rownames(rcd8[rcd8[,2]>0.00002 & rcd8[,2] > rcd8[,1]*2,])
rCD8NonHVG <- setdiff(rownames(data[data[,8]>0|data[,10]>0,]), rCD8HVG)
dcd4 <-  data[,c("D4U","D4L")]
dcd4 <-  normalize(dcd4)
rCD4GVH <-  rownames(dcd4[dcd4[,2]>0.00002 & dcd4[,2] > dcd4[,1]*2,])
rCD4NonGVH <- setdiff(rownames(data[data[,3]>0|data[,5]>0,]), rCD4GVH)
dcd8 <-  data[,c("D8U","D8L")]
dcd8 <-  normalize(dcd8)
rCD8GVH <-  rownames(dcd8[dcd8[,2]>0.00002 & dcd8[,2] > dcd8[,1]*2,])
rCD8NonGVH <- setdiff(rownames(data[data[,4]>0|data[,6]>0,]), rCD8GVH)
CD4HVG <- setdiff(rCD4HVG,c(rCD4NonHVG,rCD4GVH,rCD4NonGVH,rCD8HVG,rCD8NonHVG,rCD8GVH,rCD8NonGVH))
CD4NonHVG <- setdiff(rCD4NonHVG,c(rCD4HVG,rCD4GVH,rCD4NonGVH,rCD8HVG,rCD8NonHVG,rCD8GVH,rCD8NonGVH))
CD8HVG <- setdiff(rCD8HVG,c(rCD4NonHVG,rCD4GVH,rCD4NonGVH,rCD4HVG,rCD8NonHVG,rCD8GVH,rCD8NonGVH))
CD8NonHVG <- setdiff(rCD8NonHVG,c(rCD4HVG,rCD4GVH,rCD4NonGVH,rCD8HVG,rCD4NonHVG,rCD8GVH,rCD8NonGVH))

##Thank you!
