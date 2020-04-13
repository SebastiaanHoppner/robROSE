rm(list = ls())

# Load packages -----------------------------------------------------------------------------------
library(robROSE)

library(ROSE)
library(rpart)
library(ellipse)
library(ggplot2)
library(hmeasure)
library(partykit)
library(smotefamily)



# Load and preprocess data ------------------------------------------------------------------------
setwd("~/Desktop/PhD/PhD KUL/RobROSE/data/")
load("artificialdata.RData")
set.seed(2019)

X <- hacide.train
colnames(X)[1:2] <- c("X1", "X2")

X <- X[-c(1000, 999, 985),]
X[which(X$class == 1 & X$X1 < -3), ] <- cbind.data.frame(X1 = -2.6, X2 = -0.9, class = 1)
Xred <- X[-c(982, 996), ]

rm(chessdata_v1, circledata_v1, circledata_v2, circledata_v3, circledata_v4,
   hacide.train, hacide.test, rectangledata_v1, rectangledata_v2)



# Build plot function -----------------------------------------------------------------------------
printer <- function (X, title, syndata = NULL,
                     color.1.orig = "black", shape.1.orig = 18 , size.1.orig = 5,
                     color.1.syn  = "red"  , shape.1.syn  = 19,  size.1.syn  = 2) {
  Var.names.orig <- colnames(X)[1:3]
  colnames(X)[1:3] <- c("X1", "X2", "class")
  xlimit <- c(min(X[, 1:2]) - 1, max(X[, 1:2]) - 1) # c(floor(min(X[, 1])), ceiling(max(X[, 1])))
  ylimit <- c(min(X[, 1:2]) - 0.5,     max(X[, 1:2])) # c(floor(min(X[, 2])), ceiling(max(X[, 2])))
  plot.original.data <- ggplot(X, aes(x = X1, y = X2, group = class)) +
    geom_point(aes(shape = class, color = class, size = class)) +
    scale_color_manual(values = c("blue", color.1.orig)) +
    scale_shape_manual(values = c(1, shape.1.orig)) +
    scale_size_manual(values = c(1.8, size.1.orig)) +
    xlab("") +
    ylab("") +
    xlim(xlimit) +
    ylim(ylimit) +
    coord_fixed() +
    ggtitle(title) +
    theme(legend.position = "none") +
    theme(text = element_text(size = 27),
          axis.title.y = element_text(angle=0, vjust = 0.5))
  if (is.null(syndata)) {
    return(plot.original.data)
  } else {
    colnames(syndata)[1:3] = c("X1", "X2", "class")
    return(plot.original.data +
             geom_point(data = syndata, color = color.1.syn, shape = shape.1.syn, size = size.1.syn))
  }
}



# Original data -----------------------------------------------------------------------------------
printer(X, "Original data")
printer(Xred, "Reduced original data")

logit    <- glm(class ~ ., data = X,    family = "binomial")
logitRed <- glm(class ~ ., data = Xred, family = "binomial")
# scores_logit    <- logit$fitted.values
# scores_logitRed <- logitRed$fitted.values
# HMeasure(true.class = X$class,    scores_logit,    threshold = 0.5)$metrics[c("AUC", "MER", "ER", "Precision", "Recall", "F")]
# HMeasure(true.class = Xred$class, scores_logitRed, threshold = 0.5)$metrics[c("AUC", "MER", "ER", "Precision", "Recall", "F")]
slope        <- coef(logit)[2] / (-coef(logit)[3])
intercept    <- coef(logit)[1] / (-coef(logit)[3])
slopeRed     <- coef(logitRed)[2] / (-coef(logitRed)[3])
interceptRed <- coef(logitRed)[1] / (-coef(logitRed)[3])

printer(X, "Logistic regression") +
  geom_abline(intercept = intercept, slope = slope, linetype = 1, color = "black", size = 1)

printer(X, "Logistic regression") +
  geom_abline(intercept = intercept,    slope = slope,    linetype = 1, color = "black", size = 1) +
  geom_abline(intercept = interceptRed, slope = slopeRed, linetype = 2, color = "black", size = 1)


cart    <- as.party(rpart(class ~ ., data = X,    control = rpart.control(maxdepth = 3, minsplit = 4)))
cartRed <- as.party(rpart(class ~ ., data = Xred, control = rpart.control(maxdepth = 3, minsplit = 4)))
plot(cart)
plot(cartRed)

seg_df <- data.frame(X1    = c(-1.496),
                     X2    = c(-2.648),
                     xend = c(Inf),
                     yend = c(-2.648),
                     class = factor(1))
seg_df_red <- data.frame(X1    = c(-1.5),
                         X2    = c(-2.65),
                         xend = c(Inf),
                         yend = c(-2.65),
                         class = factor(1))

printer(X, "Classification tree (CART)") +
  geom_vline(xintercept = -1.496, linetype = 1, size = 1) +
  geom_segment(data = seg_df, aes(X1, X2, xend = xend, yend = yend), linetype = 1, size = 1)

printer(X, "Classification tree (CART)") +
  geom_vline(xintercept = -1.496, linetype = 1, size = 1) +
  geom_segment(data = seg_df, aes(X1, X2, xend = xend, yend = yend), linetype = 1, size = 1) +
  geom_vline(xintercept = -1.5, linetype = 2, size = 1) +
  geom_segment(data = seg_df_red, aes(X1, X2, xend = xend, yend = yend), linetype = 2, size = 1)




# SMOTE -------------------------------------------------------------------------------------------
set.seed(2019)
dataSMOTE <- SMOTE(X = X[, 1:2], target = X$class, K = 5, dup_size = 0)$syn_data
printer(X, syndata = dataSMOTE, title = "SMOTE")



# ROSE --------------------------------------------------------------------------------------------
nVar <- ncol(X)
n0 <- max(table(X$class))
n1 <- min(table(X$class))
i.orig.1 <- which(X$class == 1)

set.seed(2019)
dataROSE1 <- ROSE(class ~ ., data = X, N = nrow(X), seed = 2019,
                  p = 0.999, hmult.majo = 0.001, hmult.mino = 1)$data
dataROSE1 <- dataROSE1[dataROSE1$class == 1, ][1:(n0-n1), ]
printer(X, syndata = dataROSE1, title = "ROSE (h.minor = 1) (default)")


set.seed(2019)
h.minor <- 0.15
dataROSE2 <- ROSE(class ~ ., data = X, N = nrow(X), seed = 2019,
                  p = 0.999, hmult.majo = 0.001, hmult.mino = h.minor)$data
dataROSE2 <- dataROSE2[dataROSE2$class == 1, ][1:(n0-n1), ]
printer(X, syndata = dataROSE2, title = paste0("ROSE (h.minor = ", h.minor, ")"))


Sigma1 <- diag(diag(cov(X[i.orig.1, 1:2])))
level <- 0.03
plot.ROSE2 <- printer(X, syndata = dataROSE2, title = paste0("ROSE (h.minor = ", h.minor, ")"))
for (i in 1:17) {
  ellipsedata <- cbind.data.frame(ellipse(x = Sigma1, centre = c(X[i.orig.1, 1:2][i, 1], X[i.orig.1, 1:2][i, 2]),
                                          level = level, t = sqrt(qchisq(level, nVar-1)), npoints = 1000),
                                  factor(rep(1, 1000), levels = c(0,1)))
  colnames(ellipsedata) = c("X1", "X2", "class")
  plot.ROSE2 <- plot.ROSE2  + geom_path(data = ellipsedata, color = "black")
}
plot.ROSE2



# robROSE -----------------------------------------------------------------------------------------
set.seed(2019)
RobROSE <- robROSE(class ~ ., data = X, r = 0.5,
                   type = "equalsize", const = 1.5, seed = 2019)
printer(X, syndata = RobROSE$data, title = "robROSE")


robSigma1 <- covMcd(X[i.orig.1,1:2])$cov
level <- RobROSE$hmult * 0.6
level[which.max(level)] <- level[which.max(level)] / 0.6 * 1
plot.robROSE <- printer(X, syndata = RobROSE$data, title = "robROSE")
i_selected <- c(1,3,4,5,6,7,8,9,10,11,12,13,14,15,17)
for (i in i_selected) {
  ellipsedata <- cbind.data.frame(ellipse(x = robSigma1, centre = c(X[i.orig.1, 1:2][i, 1], X[i.orig.1, 1:2][i, 2]),
                                          level = level[which(i_selected == i)], t = sqrt(qchisq(level[which(i_selected == i)], nVar-1)), npoints = 1000),
                                  factor(rep(1, 1000), levels = c(0,1)))
  colnames(ellipsedata) = c("X1", "X2", "class")
  plot.robROSE <- plot.robROSE  + geom_path(data = ellipsedata, color = "black")
}
plot.robROSE


