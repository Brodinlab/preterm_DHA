library(dplyr)
library(randomForest)
library(ROSE)
library(caret)
library(readr)
library(pROC)

df <- read_csv('data/normalized_rf_input.csv') %>%
  mutate(f_Breastmilk_group = if_else(f_Breastmilk_group=='donated', 'donated', 'MOM'),
         f_preterm_class = as.factor(f_preterm_class))
set.seed(42)
ctrl <- trainControl(method='repeatedcv', number=5, repeats=3, classProbs =T, summaryFunction = twoClassSummary, savePredictions = "all")
# ctrl <- trainControl(method="cv", 
#                      summaryFunction=twoClassSummary, 
#                      classProbs=T,
#                      savePredictions = T)
rose.rf <- function(df, ntree=1000){
  set.seed(42)
  train.rose <- ROSE(f_Breastmilk_group ~., data = df, N=400)$data # ROSE for unbalanced data set
  set.seed(42)
  rf.model <- train(
    f_Breastmilk_group~., 
    data=train.rose, 
    method='rf', 
    trControl=ctrl,
    metric='ROC',
    # params for RF
    # na.action = na.omit, 
    ntree = ntree,
    importance = TRUE,
  )
}

# model1 
m1 <- rose.rf(df)
pdf('figures/Figure 7/7b.pdf')
varImpPlot(m1$finalModel, 
           sort = T,
           n.var = 30,
           main = "Top 30 - Variable Importance")
dev.off()

# model2 remove NK clusters
df2 <- df %>% dplyr::select(-c("f_NK_67", "f_NK_75", "f_NK_84"))
m2 <- rose.rf(df2)
# varImpPlot(m2$finalModel, 
#            sort = T,
#            n.var = 30,
#            main = "Top 30 - Variable Importance")

plot.auc<- function(model, df){
  predictions <- as.data.frame(predict(model, df, type = "prob"))
  predictions$predict <- names(predictions)[1:2][apply(predictions[,1:2], 1, which.max)]
  predictions$observed <- df$f_Breastmilk_group
  predictions$ob_prob <- sapply(c(1:dim(predictions)[1]), function(x) {predictions[x, as.character(predictions[x, 'observed'])]})
  roc.res <- pROC::roc(predictions$observed, predictions$ob_prob, auc=TRUE)
  plot(roc.res, col = "gray60")
  title(paste0('AUC=', round(roc.res$auc, digits = 3)))
  return(roc.res)
}
auc1 <- plot.auc(m1,df)
auc2 <- plot.auc(m2,df2)

df <- rbind(data.frame(specificity = auc1$specificities,
                       sensitivity = auc1$sensitivities,
                       model = 'complete model'),
            data.frame(specificity = auc2$specificities,
                       sensitivity = auc2$sensitivities,
                       model = 'no NK subclusters'))
ggplot(df, aes(x=specificity, y=sensitivity, color=model)) + 
  geom_line() + theme_classic() +
  geom_abline(intercept = 1, slope = -1, color = 'grey')
ggsave('figures/Figure 7/7a.pdf')


