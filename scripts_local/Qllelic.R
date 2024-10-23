library(Qllelic)

########################

sample = 'GM19102'

intable <- paste(file = '/path/to/Qllelic_input_tables/', sample, '_count.txt', sep = '')

designMatrix <- BuildDesign(experimentNames = 'sample', techReps = c(3)) # change techReps depending on the number of replicates

condition <- which(designMatrix$experimentNames == "sample")

geneCountTab <- GetGatkPipelineTabs(inFiles = intable, nReps = designMatrix$techReps)

aiTable <- Reduce(function(x,y){merge(x, y, by="ID")},
                  lapply(unlist(designMatrix$replicateNums[condition]), function(i){
                    CountsToAI(df = geneCountTab, reps = i, thr=10)
                  })
)

colnames(aiTable) <- c("ID", paste("Gene AI - replicate", unlist(designMatrix$replicateNums[condition])))

covTable <- MeanCoverage(df = geneCountTab, reps = unlist(designMatrix$replicateNums[condition]), thr=10)
colnames(covTable) <- c("ID", "Mean Allelic Coverage")

aicovTable <- merge(aiTable, covTable, by="ID")


QCCs <- lapply(1:length(designMatrix$techReps), function(i){
  resCC <- ComputeCorrConstantsForAllPairsReps(inDF = geneCountTab,
                                               vectReps = unlist(designMatrix$replicateNums[i]),
                                               EPS = 1.05, thr = 10)
  sapply(resCC, function(x){x$fittedCC})
})

designMatrix <- BuildDesign(experimentNames = 'sample', techReps = c(3),
                            corrConst = sapply(QCCs, function(cc){mean(cc)})) #change techReps depending on the number of replicates

AICI_tables <- lapply(1:length(designMatrix$techReps), function(i){
  PerformBinTestAIAnalysisForConditionNPoint_knownCC(inDF = geneCountTab, vectReps = designMatrix$replicateNums[i], 
                                                     vectRepsCombsCC = designMatrix$corrConst[i], thr = 10)
})

write.table(AICI_tables[[1]], file = paste('/path/to/Qllelic_results/', sample, '_results.txt', sep = ''), sep = '\t')


