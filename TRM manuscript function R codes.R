#Fig.1C 1D 1F & Fig.S4 
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

#example for function REP_total
#D4U pre-Tx donor unstimulated CD4
#D8U pre-Tx donor unstimulated CD8
#D4L pre-Tx donor CD4 CDSE low from donor vs recipient MLR
#D8L pre-Tx donor CD8 CDSE low from donor vs recipient MLR
#R4U pre-Tx recipient unstimulated CD4
#R8U pre-Tx recipient unstimulated CD8
#R4L pre-Tx recipient CD4 CDSE low from recipient vs donor MLR
#R8L pre-Tx recipient CD8 CDSE low from recipient vs donor MLR
colnames(Pt14_reorder)
[1] "Sum..Templates."                     "Present.In"                          "Pt14 D4U --- 168380"                 "Pt14 D8U --- 172991"                
[5] "Pt14 D4L --- 107684"                 "Pt14 D8L --- 83220"                  "Pt14 R4U --- 218304"                 "Pt14 R8U --- 169987"                
[9] "Pt14 R4L --- 73163"                  "Pt14 R8L --- 111507"                 "Pt14 PBMC POD30 --- 33681"           "Pt14 PBMC POD456 --- 50353"         
[13] "Pt14 PBMC POD1764 --- 23755"         "Pt14 ileum POD9 --- 9392"            "Pt14 ileum POD16 --- 50586"          "Pt14 ileum POD156 --- 9970"         
[17] "Pt14 ileum POD226 --- 27639"         "Pt14 ileum POD527 --- 22714"         "Pt14 ileum POD717 --- 46427"         "Pt14 ileum POD1764 --- 47833"       
[21] "Pt14 colon POD527 --- 830"           "Pt14 colon POD717 --- 13482"         "Pt14 colon POD1764 --- 855"          "Pt14 native colon POD527 --- 32172" 
[25] "Pt14 native colon POD717 --- 7311"   "Pt14 native colon  POD1764 --- 3075" "Pt14 PBMC POD456 R4U --- 18655"      "Pt14 PBMC POD456 R4L --- 1810"      
[29] "Pt14 PBMC POD456 R8U --- 14693"      "Pt14 PBMC POD456 R8L --- 791"

data <- Pt14_reorder[,c(14,15,16,17,18,19)]
pre_donor <- rownames(Pt14_reorder[Pt14_reorder[,3]>0|Pt14_reorder[,4]>0|Pt14_reorder[,5]>0|Pt14_reorder[,6]>0, ])
pre_reci <- rownames(Pt14_reorder[Pt14_reorder[,7]>0|Pt14_reorder[,8]>0|Pt14_reorder[,9]>0|Pt14_reorder[,10]>0, ])
matrix <- REP_total(data, pre_donor, pre_reci)

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

#Fig.2A Fig.S6A
Share <- function(data, A, B, donor){
  matrix <- as.data.frame(matrix(nrow = 5, ncol = ncol(data)))
  rownames(matrix) <- c("A", "B", "ABshared", "unmappable", "donor")
  colnames(matrix) <- colnames(data)
  
  error <- intersect(donor, union(A, B))
  Aonly <- setdiff(A, union(B, donor))
  Bonly <- setdiff(B, union(A, donor))
  share <- setdiff(intersect(A, B),donor)
  donor <- setdiff(donor, union(A, B))
  data <- data[!(rownames(data) %in% c(error, donor)),]
  unmappable <- setdiff(rownames(data), c(A, B, donor))
  data1 <- normalize(data)
  
  for (i in 1:ncol(data)){
    matrix[1,i] <- sum(data1[Aonly,i])*100
    matrix[2,i] <- sum(data1[Bonly,i])*100
    matrix[3,i] <- sum(data1[share,i])*100
    matrix[4,i] <- sum(data1[unmappable,i])*100
    matrix[5,i] <- sum(data1[donor,i])*100
  }
  return(matrix)
}

#example for function Share
colnames(Pt20_new)
[1] "Sum..Templates."                    "Present.In"                         "Pt20 SP D4U --- 86451"              "Pt20 SP D8U --- 56329"             
[5] "Pt20 SP D4L --- 67406"              "Pt20 SP D8L --- 57582"              "Pt20 MLN R4U --- 186655"            "Pt20 MLN R8U --- 157716"           
[9] "Pt20 MLN R4L --- 76200"             "Pt20 MLN R8L --- 60132"             "Pt20 native colon POD0 --- 212129"  "Pt20 PBMC POD16 --- 46030"         
[13] "Pt20 PBMC POD26 --- 906"            "Pt20 PBMC POD39 --- 33287"          "Pt20 PBMC POD68 --- 14467"          "Pt20 PBMC POD110 --- 25611"        
[17] "Pt20 PBMC POD242 --- 16612"         "Pt20 PBMC POD306 --- 20737"         "Pt20 PBMC POD524 T --- 68925"       "Pt20 PBMC POD635 --- 79855"        
[21] "Pt20 PBMC POD796 --- 20218"         "Pt20 ileum POD16 --- 479"           "Pt20 ileum POD23 --- 30225"         "Pt20 ileum POD37 --- 108"          
[25] "Pt20 ileum POD68 --- 3510"          "Pt20 ileum POD104 --- 41422"        "Pt20 ileum POD250 T --- 3497"       "Pt20 ileum POD306 --- 25379"       
[29] "Pt20 ileum POD523 T --- 4419"       "Pt20 ileum POD796 --- 371119"       "Pt20 BM POD125 T --- 27452"         "Pt20 BM POD521 T --- 15616"        
[33] "Pt20 colon POD250 T --- 1291"       "Pt20 native colon POD250 T --- 116"

data <- Pt20_new[,c(25,26,27,28,29,30)]
lymphoid <- rownames(Pt20_new[Pt20_new[,7]>0|Pt20_new[,8]>0|Pt20_new[,9]>0|Pt20_new[,10]>0, ])
gut <- rownames(Pt20_new[Pt20_new[,11]>0, ])
donor <- rownames(Pt20_new[Pt20_new[,3]>0|Pt20_new[,4]>0|Pt20_new[,5]>0|Pt20_new[,6]>0, ])
matrix <- Share(data, lymphoid, gut, donor)

#Fig.2B Fig. 5A, Fig.5B, Table S4, Table S10
Share_v1 <- function(data, A, B, donor){
  l <- ncol(data)
  matrix <- matrix(nrow = l, ncol = 12)  
  colnames(matrix) <- c( "Aonly Pre-transplant clones detected", "Aonly post-Tx clones undetected", 
                         "Bonly Pre-transplant clones detected", "Bonly post-Tx clones undetected",
                         "shared Pre-transplant clones detected", "shared post-Tx clones undetected", 
                         "Aonly Bonly Odds Ratio", "Aonly Bonly P-value", 
                         "Aonly shared Odds Ratio", "Aonly shared P-value",
                         "Bonly shared Odds Ratio", "Bonly shared P-value")
  rownames(matrix) <- colnames(data)
  
  error <- intersect(donor, union(A, B))
  Aonly <- setdiff(A, union(B, donor))
  Bonly <- setdiff(B, union(A, donor))
  share <- setdiff(intersect(A, B),donor)
  donor <- setdiff(donor, union(A, B))
  data <- data[!(rownames(data) %in% c(error, donor)),]
  unmappable <- setdiff(rownames(data), c(A, B, donor))
  
  for (i in 1:l) {
    matrix[i,1] <- length(Aonly)
    matrix[i,2] <- length1(data[Aonly,i])
    matrix[i,3] <- length(Bonly)
    matrix[i,4] <- length1(data[Bonly,i])
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

#example for function Share_v1
colnames(Pt20_new)
[1] "Sum..Templates."                    "Present.In"                         "Pt20 SP D4U --- 86451"              "Pt20 SP D8U --- 56329"             
[5] "Pt20 SP D4L --- 67406"              "Pt20 SP D8L --- 57582"              "Pt20 MLN R4U --- 186655"            "Pt20 MLN R8U --- 157716"           
[9] "Pt20 MLN R4L --- 76200"             "Pt20 MLN R8L --- 60132"             "Pt20 native colon POD0 --- 212129"  "Pt20 PBMC POD16 --- 46030"         
[13] "Pt20 PBMC POD26 --- 906"            "Pt20 PBMC POD39 --- 33287"          "Pt20 PBMC POD68 --- 14467"          "Pt20 PBMC POD110 --- 25611"        
[17] "Pt20 PBMC POD242 --- 16612"         "Pt20 PBMC POD306 --- 20737"         "Pt20 PBMC POD524 T --- 68925"       "Pt20 PBMC POD635 --- 79855"        
[21] "Pt20 PBMC POD796 --- 20218"         "Pt20 ileum POD16 --- 479"           "Pt20 ileum POD23 --- 30225"         "Pt20 ileum POD37 --- 108"          
[25] "Pt20 ileum POD68 --- 3510"          "Pt20 ileum POD104 --- 41422"        "Pt20 ileum POD250 T --- 3497"       "Pt20 ileum POD306 --- 25379"       
[29] "Pt20 ileum POD523 T --- 4419"       "Pt20 ileum POD796 --- 371119"       "Pt20 BM POD125 T --- 27452"         "Pt20 BM POD521 T --- 15616"        
[33] "Pt20 colon POD250 T --- 1291"       "Pt20 native colon POD250 T --- 116"

data <- Pt20_new[,c(25,26,27,28,29,30)]
lymphoid <- rownames(Pt20_new[Pt20_new[,7]>0|Pt20_new[,8]>0|Pt20_new[,9]>0|Pt20_new[,10]>0, ])
gut <- rownames(Pt20_new[Pt20_new[,11]>0, ])
donor <- rownames(Pt20_new[Pt20_new[,3]>0|Pt20_new[,4]>0|Pt20_new[,5]>0|Pt20_new[,6]>0, ])
matrix <- Share_v1(data, lymphoid, gut, donor)

#Fig.S7A
library(tidyverse)

data_share <- data[share,]
A <- as.data.frame(table(data_share[,1]))
colnames(aaa) <- c("copy_number", "A") 

B <- as.data.frame(table(data_share[,2]))
colnames(bbb) <- c("copy_number", "B") 

C <- as.data.frame(table(data_share[,3]))
colnames(ccc) <- c("copy_number", "C") 

matrix <- full_join(A, full_join(B, C))


data_gut_only <- data[gut_only,]
A <- as.data.frame(table(data_gut_only[,1]))
colnames(aaa) <- c("copy_number", "A") 

B <- as.data.frame(table(data_gut_only[,2]))
colnames(bbb) <- c("copy_number", "B") 

C <- as.data.frame(table(data_gut_only[,3]))
colnames(ccc) <- c("copy_number", "C") 

matrix <- full_join(A, full_join(B, C))


data_lymphoid_only <- data[gut_only,]
A <- as.data.frame(table(data_lymphoidt_only[,1]))
colnames(aaa) <- c("copy_number", "A") 

B <- as.data.frame(table(data_lymphoid_only[,2]))
colnames(bbb) <- c("copy_number", "B") 

C <- as.data.frame(table(data_lymphoid_only[,3]))
colnames(ccc) <- c("copy_number", "C") 

matrix <- full_join(A, full_join(B, C))

#example 
names(Pt18_new)
data <- Pt18_new[,c(23,24,25)]
colnames(data)
[1] "Pt18 ileum POD105 --- 6742"  "Pt18 ileum POD307 --- 11692" "Pt18 ileum POD357 --- 96857"
donor <- rownames(Pt18_new[Pt18_new[,3]>0|Pt18_new[,4]>0|Pt18_new[,5]>0|Pt18_new[,6]>0, ])
gut <- rownames(Pt18_new[Pt18_new[,11]>0, ])
lymphoid <- rownames(Pt18_new[Pt18_new[,7]>0|Pt18_new[,8]>0|Pt18_new[,9]>0|Pt18_new[,10]>0, ])

error <- intersect(donor, union(gut, lymphoid))
gutonly <- setdiff(gut, union(lymphoid, donor))
lymphoidonly <- setdiff(lymphoid, union(gut, donor))
share <- setdiff(intersect(gut, lymphoid),donor)
donor <- setdiff(donor, union(gut, lymphoid))
data <- data[!(rownames(data) %in% c(error, donor)),]

sub <- rownames(Pt18_new[Pt18_new[,7]>0|Pt18_new[,8]>0, ])
sub2 <- intersect(sub, lymphoidonly)
sub2_share <- intersect(sub, share)

data_sub2_lymphoidonly <- data[sub2,]
A <- as.data.frame(table(data_sub2_lymphoidonly[,1]))
colnames(A) <- c("copy_number", "Pt18 ileum POD105 --- 6742") 

B <- as.data.frame(table(data_sub2_lymphoidonly[,2]))
colnames(B) <- c("copy_number", "Pt18 ileum POD307 --- 11692") 

C <- as.data.frame(table(data_sub2_lymphoidonly[,3]))
colnames(C) <- c("copy_number", "Pt18 ileum POD357 --- 96857") 

matrix <- full_join(A, full_join(B, C))
matrix$copy_number <- as.numeric(as.character(matrix$copy_number))


data_sub2_gutonly <- data[gutonly,]
A <- as.data.frame(table(data_sub2_gutonly[,1]))
colnames(A) <- c("copy_number", "Pt18 ileum POD105 --- 6742") 

B <- as.data.frame(table(data_sub2_gutonly[,2]))
colnames(B) <- c("copy_number", "Pt18 ileum POD307 --- 11692") 

C <- as.data.frame(table(data_sub2_gutonly[,3]))
colnames(C) <- c("copy_number", "Pt18 ileum POD357 --- 96857") 

matrix <- full_join(A, full_join(B, C))
matrix$copy_number <- as.numeric(as.character(matrix$copy_number))


data_sub2_share <- data[sub2_share,]
A <- as.data.frame(table(data_sub2_share[,1]))
colnames(A) <- c("copy_number", "Pt18 ileum POD105 --- 6742") 

B <- as.data.frame(table(data_sub2_share[,2]))
colnames(B) <- c("copy_number", "Pt18 ileum POD307 --- 11692") 

C <- as.data.frame(table(data_sub2_share[,3]))
colnames(C) <- c("copy_number", "Pt18 ileum POD357 --- 96857") 

matrix <- full_join(A, full_join(B, C))
matrix$copy_number <- as.numeric(as.character(matrix$copy_number))

#Fig.S7B Table S5
Share_v1_R4U_R8U <- function(data, gut, lymphoid, donor, R4U_R8U){
  l <- ncol(data)
  matrix <- matrix(nrow = l, ncol = 12)  
  colnames(matrix) <- c( "gutonly Pre-transplant clones detected", "gutonly post-Tx clones undetected", 
                         "lymphoidonly Pre-transplant clones detected", "lymphoidonly post-Tx clones undetected",
                         "shared Pre-transplant clones detected", "shared post-Tx clones undetected", 
                         "gutonly lymphoidonly Odds Ratio", "gutonly lymphoidonly P-value", 
                         "gutonly shared Odds Ratio", "gutonly shared P-value",
                         "lymphoidonly shared Odds Ratio", "lymphoidonly shared P-value")
  rownames(matrix) <- colnames(data)
  
  error <- intersect(donor, union(gut, lymphoid))
  gutonly <- setdiff(gut, union(lymphoid, donor))
  lymphoidonly <- setdiff(lymphoid, union(gut, donor))
  share <- setdiff(intersect(gut, lymphoid),donor)
  donor <- setdiff(donor, union(gut, lymphoid))
  data <- data[!(rownames(data) %in% c(error, donor)),]
  unmappable <- setdiff(rownames(data), c(gut, lymphoid, donor))
  R4U_R8U2 <- intersect(R4U_R8U, lymphoidonly)
  share <- intersect(R4U_R8U, share)
  
  share <- intersect(share, rownames(data[(data[,1] + data[,2])==1 & data[,5]==1, ]))
  R4U_R8U2 <- intersect(R4U_R8U2, rownames(data[(data[,1] + data[,2]) ==1, ]))
  gutonly <- intersect(gutonly, rownames(data[data[,5]==1, ]))
  
  matrix[,1] <- length1(data[gutonly,5])
  matrix[,3] <- length1(data[R4U_R8U2,1])
  matrix[,5] <- length1(data[share,5])
  
  for (i in 1:l) {
    matrix[i,2] <- length1(data[gutonly,i])
    matrix[i,4] <- length1(data[R4U_R8U2,i])
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

#example for function Share_v1_R4U_R8U
names(Pt18_new)
[1] "Sum..Templates."                    "Present.In"                         "Pt18 SP D4U --- 80509"              "Pt18 SP D8U --- 37862"             
[5] "Pt18 SP D4L --- 43547"              "Pt18 SP D8L --- 80174"              "Pt18 SP R4U --- 41891"              "Pt18 SP R8U --- 145868"            
[9] "Pt18 SP R4L --- 46336"              "Pt18 SP R8L --- 82850"              "Pt18 intestine POD0 --- 4628"       "Pt18 PBMC POD5 --- 1159"           
[13] "Pt18 PBMC POD7 --- 1817"            "Pt18 PBMC POD14 --- 28858"          "Pt18 PBMC POD28 --- 40747"          "Pt18 PBMC POD98 --- 108893"        
[17] "Pt18 PBMC POD314 --- 57969"         "Pt18 PBMC POD357 --- 26279"         "Pt18 PBMC POD357 D --- 2578"        "Pt18 ileum POD5.7 --- 6326"        
[21] "Pt18 ileum POD18 --- 3934"          "Pt18 ileum POD35 --- 456"           "Pt18 ileum POD105 --- 6742"         "Pt18 ileum POD307 --- 11692"       
[25] "Pt18 ileum POD357 --- 96857"        "Pt18 stomach POD357 --- 17204"      "Pt18 duodenum POD357 --- 7284"      "Pt18 colon POD357 --- 14076"       
[29] "Pt18 native colon POD357 --- 17731" "Pt18 BM POD357 --- 4298"            "Pt18 BM POD357 D --- 1430"   
data <- Pt18_new[,c(7,8,9,10, 11, 23,24,25)]
donor <- rownames(Pt18_new[Pt18_new[,3]>0|Pt18_new[,4]>0|Pt18_new[,5]>0|Pt18_new[,6]>0, ])
gut <- rownames(Pt18_new[Pt18_new[,11]>0, ])
lymphoid <- rownames(Pt18_new[Pt18_new[,7]>0|Pt18_new[,8]>0|Pt18_new[,9]>0|Pt18_new[,10]>0, ])
sub <- rownames(Pt18_new[Pt18_new[,7]>0|Pt18_new[,8]>0, ])
matrix <- Share_v1_R4U_R8U(data, gut, lymphoid, donor, sub)

#Fig.2C 2D 
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
data[(data$blood>0) & (data$tissue == 0), c("Type")] <- c("Blood_Only")
data[(data$blood == 0) & (data$tissue>0), c("Type")] <- c("ileum_Only")
data[(data$blood>0) & (data$tissue>0), c("Type")] <- c("Shared")
SharedinGut <- sum(data[data$Type == c("Shared"),c("tissue")])/sum(data[,c("tissue")])*100
data[, c("Blood")] <- data[, c("blood")]
data[, c("Tissue")] <- data[, c("tissue")]
data[data$blood == 0, c("Blood")] <- 0.000001
data[data$tissue == 0, c("Tissue")] <- 0.000001


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

#example
colnames(HC_UC)
[1] "Sum..Templates."             "Present.In"                  "C01_Blood_control"           "C01_Tissue_control"          "C02_Blood_control"           "C02_Tissue_control"         
[7] "C03_Blood_control"           "C03_Tissue_control"          "C04_Blood_control"           "C04_Tissue_control"          "P01_Blood_mild_disease"      "P01_Tissue_mild_disease"    
[13] "P02_Blood_mild_disease"      "P02_Tissue_mild_disease"     "P03_Blood_moderate_disease"  "P03_Tissue_moderate_disease" "P04_Blood_moderate_disease"  "P04_Tissue_moderate_disease"
[19] "P05_Blood_moderate_disease"  "P05_Tissue_moderate_disease" "P06_Blood_moderate_disease"  "P06_Tissue_moderate_disease" "P07_Blood_severe_disease"    "P07_Tissue_severe_disease"  
[25] "P08_Blood_severe_disease"    "P08_Tissue_severe_disease"   "P09_Blood_severe_disease"    "P09_Tissue_severe_disease" 
data <- HC_UC[,c(15,16)]
data <- normalize(data)
data <- data[which(rowSums(data) > 0),]
data$Type <- NA
data[(data$P03_Blood_moderate_disease>0) & (data$P03_Tissue_moderate_disease == 0), c("Type")] <- c("Blood_Only")
data[(data$P03_Blood_moderate_disease == 0) & (data$P03_Tissue_moderate_disease>0), c("Type")] <- c("ileum_Only")
data[(data$P03_Blood_moderate_disease>0) & (data$P03_Tissue_moderate_disease>0), c("Type")] <- c("Shared")
SharedinGut <- sum(data[data$Type == c("Shared"),c("P03_Tissue_moderate_disease")])/sum(data[,c("P03_Tissue_moderate_disease")])*100
data[, c("Blood")] <- data[, c("P03_Blood_moderate_disease")]
data[, c("Tissue")] <- data[, c("P03_Tissue_moderate_disease")]
data[data$P03_Blood_moderate_disease == 0, c("Blood")] <- 0.000001
data[data$P03_Tissue_moderate_disease == 0, c("Tissue")] <- 0.000001


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

#Fig.3B
matrix <- as.data.frame(matrix(nrow=length(data) , ncol=length(data) ))
colnames(matrix) <- colnames(data)
rownames(matrix) <- colnames(data)
for (m in 1:ncol(matrix)){
  for (n in 1:nrow(matrix)) {
    matrix[n,m] <- sum(data[intersect(rownames(data)[data[,n]>0], rownames(data)[data[,m]>0]),n])/sum(data[,n])*100
  }
}

#Fig.3C Fig.4D 4E Fig. S10B
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

#Fig.4B 4C Fig.S8 Fig.S10A
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

#HvG/NonHvG define in Fig.5 & Fig.S11
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

#ITx pts bulk sequencing data match path_data in Fig.5 & Fig.S11
library(dplyr)
library(tidyverse)
data <- unique(left_join(data, path_data))

#Fig. S9
library(Seurat)

#read data
mj026_data <- Read10X(data.dir = "/Users/MJ26")
mj027_data <- Read10X(data.dir = "/Users/MJ27")
mj028_data <- Read10X(data.dir = "/Users/MJ28")

#wash data
mj026_me <- CreateSeuratObject(counts = mj026_data, project = "Pt21_Pre_Spl", min.cells = 5)
mj026_me$pos <- "Pt21_Pre_Spl"
mj026_me$mt <- PercentageFeatureSet(mj026_me, pattern = "^MT-")
mj026_me_Q1 = quantile(mj026_me$nFeature_RNA)[2]
mj026_me_Q3 = quantile(mj026_me$nFeature_RNA)[4]
mj026_me_interQ = (mj026_me_Q3-mj026_me_Q1)*1.5
mj026_me_upperb = mj026_me_Q3 + mj026_me_interQ
mj026_me_lowerb = mj026_me_Q1 - mj026_me_interQ
mj026_me <- subset(mj026_me, subset = nFeature_RNA > mj026_me_lowerb & nFeature_RNA < mj026_me_upperb  & mt < 15)

mj027_me <- CreateSeuratObject(counts = mj027_data, project = "Pt21_Pre_IEL", min.cells = 5)
mj027_me$pos <- "Pt21_Pre_IEL"
mj027_me$mt <- PercentageFeatureSet(mj027_me, pattern = "^MT-")
mj027_me_Q1 = quantile(mj027_me$nFeature_RNA)[2]
mj027_me_Q3 = quantile(mj027_me$nFeature_RNA)[4]
mj027_me_interQ = (mj027_me_Q3-mj027_me_Q1)*1.5
mj027_me_upperb = mj027_me_Q3 + mj027_me_interQ
mj027_me_lowerb = mj027_me_Q1 - mj027_me_interQ
mj027_me <- subset(mj027_me, subset = nFeature_RNA > mj027_me_lowerb & nFeature_RNA < mj027_me_upperb  & mt < 15)

mj028_me <- CreateSeuratObject(counts = mj028_data, project = "Pt21_Pre_LPL", min.cells = 5)
mj028_me$pos <- "Pt21_Pre_LPL"
mj028_me$mt <- PercentageFeatureSet(mj028_me, pattern = "^MT-")
mj028_me_Q1 = quantile(mj028_me$nFeature_RNA)[2]
mj028_me_Q3 = quantile(mj028_me$nFeature_RNA)[4]
mj028_me_interQ = (mj028_me_Q3-mj028_me_Q1)*1.5
mj028_me_upperb = mj028_me_Q3 + mj028_me_interQ
mj028_me_lowerb = mj028_me_Q1 - mj028_me_interQ
mj028_me <- subset(mj028_me, subset = nFeature_RNA > mj028_me_lowerb & nFeature_RNA < mj028_me_upperb  & mt < 15)

#combine data
ifnb.list_pre <- list(mj026_me, mj027_me, mj028_me)

ifnb.list_pre <- lapply(X = ifnb.list_pre, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 20000)
})

features_pre <- SelectIntegrationFeatures(object.list = ifnb.list_pre, nfeatures = 20000)

immune.anchors_pre <- FindIntegrationAnchors(object.list = ifnb.list_pre, anchor.features = features_pre)
immune.combined_pre <- IntegrateData(anchorset = immune.anchors_pre)

DefaultAssay(immune.combined_pre) <- "integrated"

#Fig. S8a
DimPlot(immune.combined_pre, reduction = "umap")

#Fig. S8b
DoHeatmap(immune.combined_pre, features = c("CD3E", "CD4", "CD8A", "CD8B", "TRDC",
                                            "CD69", "ITGAE", "ITGA1", "CXCR6", "RUNX3", "PRDM1", 
                                            "CCR7", "KLF2", "S1PR1", "SELL"),
          size = 5, angle = 0, hjust = 0.5) +
  NoLegend()

#Fig. S8c & Fig. S8d
mj26_contig <- read.csv("~/Desktop/MJ26/all_contig_annotations.csv")
mj27_contig <- read.csv("~/Desktop/all_contig_annotations.csv")
mj28_contig <- read.csv("~/Desktop/all_contig_annotations.csv")
Pt21_new_truncted <- read.csv("~/Desktop/bulk sequencing/Pt21_new_truncted.csv")
ITx_Pt21_aamatch<- read.csv("~/Desktop/bulk sequencing/ITx_Pt21_aamatch.csv")

#shared defined from 10X TCR
Pt21_preSpl <- mj26_contig[mj26_contig$chain == c("TRB") & mj26_contig$productive == c("true") & mj26_contig$is_cell == c("true") & mj26_contig$high_confidence == c("true"),
                           c("barcode","v_gene","j_gene","cdr3")]
Pt21_preIEL <- mj27_contig[mj27_contig$chain == c("TRB") & mj27_contig$productive == c("true") & mj27_contig$is_cell == c("true") & mj27_contig$high_confidence == c("true"),
                           c("barcode","v_gene","j_gene","cdr3")]
Pt21_preLPL <- mj28_contig[mj28_contig$chain == c("TRB") & mj28_contig$productive == c("true") & mj28_contig$is_cell == c("true") & mj28_contig$high_confidence == c("true"),
                           c("barcode","v_gene","j_gene","cdr3")]

#shared defined from bulk TCR
library(tidyverse)

bulk_gut <- rownames(Pt21_new_truncted[Pt21_new_truncted[,11]>0,])
bulk_spl <- rownames(Pt21_new_truncted[Pt21_new_truncted[,7]>0|Pt21_new_truncted[,8]>0|Pt21_new_truncted[,9]>0|Pt21_new_truncted[,10]>0,] )

ITx_Pt21_aamatch_pre <- ITx_Pt21_aamatch[ITx_Pt21_aamatch$sample_name %in% c("Pt21_Pre_Recipient_CD4_CFSElo", "Pt21_Pre_Recipient_CD4_unstim",              
                                                                             "Pt21_Pre_Recipient_CD8_CFSElo", "Pt21_Pre_Recipient_CD8_unstim",
                                                                             "Pt21_MVTx_POD0_recip_ileum_Bx"),]
bulk_Gut <- ITx_Pt21_aamatch_pre[ITx_Pt21_aamatch_pre$rearrangement %in% bulk_gut, c("amino_acid", "v_gene", "j_gene")]
bulk_Gut <- unique(bulk_Gut)

bulk_Gut <- bulk_Gut %>%
  mutate(v_gene = str_remove_all(.$v_gene, "C")) %>%
  mutate(j_gene = str_remove_all(.$j_gene, "C"))

for (n in 1:nrow(bulk_Gut)) {
  if(bulk_Gut[n,2] %in% c("unresolved")){bulk_Gut[n,2] <- NA}
  else
  {
    bulk_Gut[n,2] <- gsub("V01-", "V1-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V02-", "V2-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V03-", "V3-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V04-", "V4-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V05-", "V5-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V06-", "V6-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V07-", "V7-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V08-", "V8-", bulk_Gut[n,2])
    bulk_Gut[n,2] <- gsub("V09-", "V9-", bulk_Gut[n,2])
  }
}

for (n in 1:nrow(bulk_Gut)) {
  bulk_Gut[n,2] <- gsub("-01", "-1", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-02", "-2", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-03", "-3", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-04", "-4", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-05", "-5", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-06", "-6", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-07", "-7", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-08", "-8", bulk_Gut[n,2])
  bulk_Gut[n,2] <- gsub("-09", "-9", bulk_Gut[n,2])
}

bulk_Gut[bulk_Gut$v_gene %in% c("TRBVA-or09_02"), 2] <- c("TRBVAor9-2")

for (n in 1:nrow(bulk_Gut)) {
  if(bulk_Gut[n,3] %in% c("unresolved")){bulk_Gut[n,3] <- NA}
  else
  {
    bulk_Gut[n,3] <- gsub("J01-", "J1-", bulk_Gut[n,3])
    bulk_Gut[n,3] <- gsub("J02-", "J2-", bulk_Gut[n,3])
  }
}

for (n in 1:nrow(bulk_Gut)) {
  bulk_Gut[n,3] <- gsub("-01", "-1", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-02", "-2", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-03", "-3", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-04", "-4", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-05", "-5", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-06", "-6", bulk_Gut[n,3])
  bulk_Gut[n,3] <- gsub("-07", "-7", bulk_Gut[n,3])
}

colnames(bulk_Gut) <- c("cdr3", "v_gene", "j_gene")

aaa <- mj26_contig[,c("barcode", "v_gene", "j_gene", "cdr3")]
bulk_Gut <- left_join(bulk_Gut, aaa)
bulk_Gut <- bulk_Gut[!is.na(bulk_Gut$barcode),]
bulk_Gut <- unique(paste0(bulk_Gut$barcode, "_1"))

bulk_Spl <- ITx_Pt21_aamatch_pre[ITx_Pt21_aamatch_pre$rearrangement %in% bulk_spl, c("amino_acid", "v_gene", "j_gene")]
bulk_Spl <- unique(bulk_Spl)

bulk_Spl <- bulk_Spl %>%
  mutate(v_gene = str_remove_all(.$v_gene, "C")) %>%
  mutate(j_gene = str_remove_all(.$j_gene, "C"))

for (n in 1:nrow(bulk_Spl)) {
  if(bulk_Spl[n,2] %in% c("unresolved")){bulk_Spl[n,2] <- NA}
  else
  {
    bulk_Spl[n,2] <- gsub("V01-", "V1-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V02-", "V2-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V03-", "V3-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V04-", "V4-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V05-", "V5-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V06-", "V6-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V07-", "V7-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V08-", "V8-", bulk_Spl[n,2])
    bulk_Spl[n,2] <- gsub("V09-", "V9-", bulk_Spl[n,2])
  }
}

for (n in 1:nrow(bulk_Spl)) {
  bulk_Spl[n,2] <- gsub("-01", "-1", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-02", "-2", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-03", "-3", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-04", "-4", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-05", "-5", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-06", "-6", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-07", "-7", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-08", "-8", bulk_Spl[n,2])
  bulk_Spl[n,2] <- gsub("-09", "-9", bulk_Spl[n,2])
}

bulk_Spl[bulk_Spl$v_gene %in% c("TRBVA-or09_02"), 2] <- c("TRBVAor9-2")
bulk_Spl[bulk_Spl$v_gene %in% c("TRBV20-or09_02"), 2] <- c("TRBV20or9-2")

for (n in 1:nrow(bulk_Spl)) {
  if(bulk_Spl[n,3] %in% c("unresolved")){bulk_Spl[n,3] <- NA}
  else
  {
    bulk_Spl[n,3] <- gsub("J01-", "J1-", bulk_Spl[n,3])
    bulk_Spl[n,3] <- gsub("J02-", "J2-", bulk_Spl[n,3])
  }
}

for (n in 1:nrow(bulk_Spl)) {
  bulk_Spl[n,3] <- gsub("-01", "-1", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-02", "-2", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-03", "-3", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-04", "-4", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-05", "-5", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-06", "-6", bulk_Spl[n,3])
  bulk_Spl[n,3] <- gsub("-07", "-7", bulk_Spl[n,3])
}

colnames(bulk_Spl) <- c("cdr3", "v_gene", "j_gene")

bbb <- mj27_contig[,c("barcode", "v_gene", "j_gene", "cdr3")]
ccc <- mj28_contig[,c("barcode", "v_gene", "j_gene", "cdr3")]
ddd <- unique(rbind(bbb, ccc))
bulk_Spl <- left_join(bulk_Spl, ddd)
bulk_Spl <- bulk_Spl[!is.na(bulk_Spl$barcode),]
bulk_Spl <- c(unique(paste0(bulk_Spl$barcode, "_2")), unique(paste0(bulk_Spl$barcode, "_3")))

#combine
shared_spl <- paste0(intersect(Pt21_preSpl$barcode, c(Pt21_preIEL$barcode, Pt21_preLPL$barcode)), "_1")
shared_gut <- c(paste0(intersect(Pt21_preSpl$barcode, Pt21_preIEL$barcode), "_2"), paste0(intersect(Pt21_preSpl$barcode, Pt21_preLPL$barcode), "_3"))
only_spl <- paste0(setdiff(Pt21_preSpl$barcode, c(Pt21_preIEL$barcode, Pt21_preLPL$barcode)), "_1")
only_gut <- c(paste0(setdiff(Pt21_preIEL$barcode, Pt21_preSpl$barcode), "_2"), paste0(setdiff(Pt21_preLPL$barcode, Pt21_preSpl$barcode), "_3"))

aaa <- as.data.frame(as.matrix(immune.combined_pre@active.ident))
colnames(aaa) <- c("cate")
aaa$cate <- c("others")
aaa[rownames(aaa) %in% c(shared_spl, bulk_Gut), c("cate")] <- c("shared_spleen") 
aaa[rownames(aaa) %in% c(shared_gut, bulk_Spl), c("cate")] <- c("shared_gut") 
aaa[rownames(aaa) %in% setdiff(only_spl, bulk_Gut), c("cate")] <- c("spleen_only") 
aaa[rownames(aaa) %in% setdiff(only_gut, bulk_Spl), c("cate")] <- c("gut_only") 

cate <- aaa$cate
immune.combined_pre$cate <- cate

Idents(immune.combined_pre) <- immune.combined_pre@meta.data$seurat_clusters
DimPlot(immune.combined_pre, reduction = "umap", pt.size = 0.3, split.by='cate')

table(Idents(immune.combined_pre), immune.combined_pre$cate)

#Fig. S8e
library(EnhancedVolcano)

deg1 = FindMarkers(immune.combined_pre, 
                   ident.1 = 'shared_gut', 
                   ident.2 = 'shared_spleen',
                   logfc.threshold = 0,
                   min.pct = 0)
write.csv(deg1,'~/Desktop/DEG shared_gut vs shared_spleen.csv')

deg2 = FindMarkers(immune.combined_pre,
                   ident.1 = 'shared_gut', 
                   ident.2 = 'gut_only',
                   logfc.threshold = 0,
                   min.pct = 0)
write.csv(deg2,'~/Desktop/DEG shared_gut vs gut_only.csv')

deg3 = FindMarkers(immune.combined_pre,
                   ident.1 = 'shared_spleen', 
                   ident.2 = 'spleen_only',
                   logfc.threshold = 0,
                   min.pct = 0)
write.csv(deg3,'~/Desktop/DEG shared_spleen vs spleen_only.csv')

EnhancedVolcano(deg1,
                lab = rownames(deg1),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 6.0,
                xlim = c(-2, 2),
                drawConnectors = T)

EnhancedVolcano(deg2,
                lab = rownames(deg2),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 6.0,
                xlim = c(-3.5, 3.5),
                drawConnectors = T,
                selectLab = c("GZMB", "GZMA", "IFNG", "PRF1", "NKG7", "GNLY", "KLRB1"))

EnhancedVolcano(deg3,
                lab = rownames(deg3),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.5,
                pointSize = 2.0,
                labSize = 6.0,
                xlim = c(-2.5, 2.5),
                drawConnectors = T,
                selectLab = c( "CXCR6", "CCR7", "KLF2", "SELL"))


##Thank you!
