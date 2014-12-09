library('ROCR')

# CPAT output file (lncRNA and mRNA) with labels in argument
args <- commandArgs(trailingOnly = TRUE)
cpatinfile <-args[1]


data=read.table(file=cpatinfile,header=T,sep="\t")
attach(data)
#names(data)

# Split in 10 fold cross validation
number_row = nrow(data)
d1 =seq(1,as.integer(number_row/10))
d2 =seq(as.integer(number_row/10)+1,as.integer(2*number_row/10))
d3 = seq(as.integer(2*number_row/10)+1,as.integer(3*number_row/10))
d4 = seq(as.integer(3*number_row/10)+1,as.integer(4*number_row/10))
d5 = seq(as.integer(4*number_row/10)+1,as.integer(5*number_row/10))
d6 = seq(as.integer(5*number_row/10)+1,as.integer(6*number_row/10))
d7 = seq(as.integer(6*number_row/10)+1,as.integer(7*number_row/10))
d8 = seq(as.integer(7*number_row/10)+1,as.integer(8*number_row/10))
d9 = seq(as.integer(8*number_row/10)+1,as.integer(9*number_row/10))
d10 = seq(as.integer(9*number_row/10)+1,number_row)


# Create list ROCR_data with the 10 datasets Response versus Labels
Response = list(data$coding_prob[d1],data$coding_prob[d2],data$coding_prob[d3],data$coding_prob[d4],data$coding_prob[d5],data$coding_prob[d6],data$coding_prob[d7],data$coding_prob[d8],data$coding_prob[d9],data$coding_prob[d10])
Labls = list(data$label[d1],data$label[d2],data$label[d3],data$label[d4],data$label[d5],data$label[d6],data$label[d7],data$label[d8],data$label[d9],data$label[d10])
ROCR_data = list(predictions=Response,Labels=Labls)

# Compute prediction with ROCR
pred <- prediction(ROCR_data$predictions, ROCR_data$Labels)

# Measure performance with tpr (sens) and tnr (spec) => see http://cran.r-project.org/web/packages/ROCR/ROCR.pdf
perf <- performance(pred,"sens","spec")

# Apply a function to get the 10 cutoffs that maximize the sens and spec : Thanks Oliver Sander
cutoffs= sapply(1:length(perf@alpha.values), function(i) { perf@alpha.values[[i]][which.max(perf@x.values[[i]]+perf@y.values[[i]])] } )

# Print 10 cutoofs
cutoffs

# Get the mean of the 10 cutoffs
mean(cutoffs)

