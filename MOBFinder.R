# remove all existing variables
rm(list=ls())

# load R package
library(randomForest)

# create function: append column
column_append <- function(idata){
  data_append <- data.frame(matrix(0, nrow = nrow(idata), ncol = 13))
  colnames(data_append) <- c(paste0("class", 1:11), "max_scores", "max_class")
  idata <-cbind(idata, data_append)
  return(idata)
}

# input file: [column 1: id] [column 2: id_length] [column 3: segment_length] [column 4-103: vector]
# create function: calculate and return classification probability
model_predict <- function(idata, model_100_400, model_401_800, model_801_1200, model_1201_1600){
  colnames(idata) <- c('C1','C2','C3',paste0('V',1:100))
  cols <- ncol(idata)
  idata <- column_append(idata)
  if(any(idata$C3 %in% 1201:1600)){
    idata_1201_1600 <- subset(idata, idata$C3 <= 1600 & idata$C3 > 1200)
    idata_1201_1600[,(cols+1):(cols+11)] <- idata_1201_1600[,3] / idata_1201_1600[,2] * as.data.frame(predict(model_1201_1600, newdata = idata_1201_1600[,c(4:103)], type = "prob"))
    idata_sum <- idata_1201_1600
  }
  if(any(idata$C3 %in% 801:1200)){
    idata_801_1200 <-subset(idata, idata$C3 <= 1200 & idata$C3 > 800)
    idata_801_1200[,(cols+1):(cols+11)] <- idata_801_1200[,3] / idata_801_1200[,2] * as.data.frame(predict(model_801_1200, newdata = idata_801_1200[,c(4:103)], type = "prob"))
    if (!exists("idata_sum")){
      idata_sum <- idata_801_1200
    }
    else{
      idata_sum <- rbind(idata_sum, idata_801_1200)
    }
  }
  if(any(idata$C3 %in% 401:800)){
    idata_401_800 <- subset(idata, idata$C3 <= 800 & idata$C3 > 400)
    idata_401_800[,(cols+1):(cols+11)] <- idata_401_800[,3] / idata_401_800[,2] * as.data.frame(predict(model_401_800, newdata = idata_401_800[,c(4:103)], type = "prob"))
    if (!exists("idata_sum")){
      idata_sum <- idata_401_800
    }
    else{
      idata_sum <- rbind(idata_sum, idata_401_800)
    }
  }
  if(any(idata$C3 %in% 4:400)){
    idata_100_400 <- subset(idata, idata$C3 <= 400)
    idata_100_400[,(cols+1):(cols+11)] <- idata_100_400[,3] / idata_100_400[,2] * as.data.frame(predict(model_100_400, newdata = idata_100_400[,c(4:103)], type = "prob"))
    if (!exists("idata_sum")){
      idata_sum <- idata_100_400
    }
    else{
      idata_sum <- rbind(idata_sum, idata_100_400)
    }
  }
  idata_sum <- idata_sum[order(as.numeric(rownames(idata_sum))),]
  idata_prob <- aggregate(idata_sum[, (cols+1):(cols+13)], by = list(idata_sum[, 1]), FUN = sum)
  colnames(idata_prob)[1:12] <- c("id", "MOBB", "MOBC", "MOBF", "MOBH", "MOBL", "MOBM", "MOBP", "MOBQ", "MOBT", "MOBV", "non-mob")
  idata_prob <- idata_prob[match(unique(idata$C1), idata_prob$id),]
  idata_prob$max_scores <- apply(idata_prob[, 2:12], MARGIN = 1, max)
  idata_prob$max_class <- names(idata_prob[, 2:12])[apply(idata_prob[, 2:12], 1, which.max)]
  idata_prob <- idata_prob[,c(1,14,2:12)]
  colnames(idata_prob)[2:13] <- c("class", "mobb_score", "mobc_score", "mobf_score", "mobh_score", "mobl_score", "mobm_score", "mobp_score", "mobq_score", "mobt_score", "mobv_score", "non-mob_score")
  rownames(idata_prob) <- NULL
  return(idata_prob)
}

# load model
model_01 <- readRDS("./model/model_100_400.rds")
model_02 <- readRDS("./model/model_401_800.rds")
model_03 <- readRDS("./model/model_801_1200.rds")
model_04 <- readRDS("./model/model_1201_1600.rds")

# get input parameter
args <- commandArgs(trailingOnly = TRUE)

# get input file list
input_list <- read.table(args[1], sep = "\t", header = FALSE)

# predict with our model, output result
# condition 1: predict contigs
# condition 2: predcit binnings
if(length(args) == 1){
  for (rown in 1:nrow(input_list)){
    input_w2v <- read.table(input_list[rown,1], sep = "\t", header = FALSE)
    pred_result <- model_predict(input_w2v, model_01, model_02, model_03, model_04)
    write.table(pred_result, file = input_list[rown,2], sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
} else{
  for (rown in 1:nrow(input_list)) {
    input_w2v <- read.table(input_list[rown,1], sep = "\t", header = FALSE)
    input_w2v$V2 <- sum(input_w2v$V2)
    pred_result <- model_predict(input_w2v, model_01, model_02, model_03, model_04)
    if (!exists("result")){
      result <- pred_result
    } else{
      result <- rbind(result, pred_result)
    }
    write.table(result, file = args[2],sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  }
}


