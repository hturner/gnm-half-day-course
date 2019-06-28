## ----setup, include = FALSE----------------------------------------------
library(knitr)
opts_chunk$set(fig.path = 'figure/beamer-', fig.align = 'center',
               fig.show = 'hold', size = 'footnotesize')
## markup inline code http://stackoverflow.com/a/16406120/173755
knit_hooks$set(inline = function(x) {
  if (is.numeric(x)) return(knitr:::format_sci(x, 'latex'))
  highr:::hi_latex(x)
})
# make the printing fit on the page
options(width = 70, digits = 3, show.signif.stars=FALSE)
par(mar = c(4, 4, .1, .1)) # reduce space above/below plots
set.seed(1121)   # make the results repeatable
library(gnm)


## ----glm, eval = FALSE---------------------------------------------------
## glm(y ~ row + col, family = poisson)


## ----quasiIndep, eval = FALSE--------------------------------------------
## y ~ row + col + Diag(row, col)


## ----quasiSymm, eval = FALSE---------------------------------------------
## y ~ row + col + Symm(row, col)


## ----Symm, eval = FALSE--------------------------------------------------
## y ~ Symm(row, col)


## ----mentalHealth--------------------------------------------------------
xtabs(count ~ SES + MHS, mentalHealth)


## ----trtContr------------------------------------------------------------
mentalHealth$MHS <- C(mentalHealth$MHS, treatment)
mentalHealth$SES <- C(mentalHealth$SES, treatment)


## ----RC------------------------------------------------------------------
RC <- gnm(count ~ SES + MHS + Mult(SES, MHS), family = poisson,
          data = mentalHealth, verbose = FALSE, ofInterest = "Mult")
coef(RC)


## ----colScores-----------------------------------------------------------
colProbs <- with(mentalHealth, tapply(count, MHS, sum) / sum(count))
colScores <- getContrasts(RC, pickCoef(RC, "[.]MHS"), ref = colProbs,
                          scaleRef = colProbs, scaleWeights = colProbs)
colScores


## ----rowScores-----------------------------------------------------------
rowProbs <- with(mentalHealth, tapply(count, SES, sum) / sum(count))
rowScores <- getContrasts(RC, pickCoef(RC, "[.]SES"), ref = rowProbs,
                          scaleRef = rowProbs, scaleWeights = rowProbs)


## ----assoc---------------------------------------------------------------
phi <- pickCoef(RC, "[.]SES", value = TRUE)
psi <- pickCoef(RC, "[.]MHS", value = TRUE)
sqrt(sum(rowProbs*(phi - sum(rowProbs*phi))^2)) *
         sqrt(sum(colProbs*(psi - sum(colProbs*psi))^2))



## ----conformity, echo = FALSE--------------------------------------------
conformity <- read.table("/media/veracrypt1/Work/Repos/CRAN/gnm-svn/DataSets/Van_der_Slik/conformity.txt",
                         colClasses = c("character", "numeric", "numeric",
                         "factor", "factor", rep("numeric", 6)))

## ----A-------------------------------------------------------------------
A <- gnm(MCFM ~ -1 +
             AGEM + MRMM + FRMF + MWORK + MFCM + Dref(MOPLM, FOPLF),
           family = gaussian, data = conformity, verbose = FALSE)


## ----w, message = FALSE--------------------------------------------------
w <- DrefWeights(A)
w


## ----wCI-----------------------------------------------------------------
w$MOPLM["weight"] + qnorm(c(0.025, 0.975)) * w$MOPLM["se"]


## ----A2, echo = FALSE----------------------------------------------------
A2 <- update(A, . ~ -1 +  AGEM + MRMM + FRMF + MWORK + MFCM + FOPLF)
anova(A2, A, test = "Chisq")


## ----F-------------------------------------------------------------------
F <- gnm(MCFM ~ -1 + AGEM + MRMM + FRMF + MWORK + MFCM +
          Dref(MOPLM, FOPLF, delta = ~ 1 + MFCM),
          family = gaussian, data = conformity, verbose = FALSE)


## ----wF, message = FALSE-------------------------------------------------
DrefWeights(F)


## ----TypeII--------------------------------------------------------------
TypeII <- function(x){
  list(predictors = list(a = 1, h = 1),
       variables = list(substitute(x)))
}
class(TypeII) <- "nonlin"


## ----paste0--------------------------------------------------------------
term = function(predLabels, varLabels){
    paste0(predLabels[1], "*", varLabels[1], "/(1 + ",
           predLabels[1], "*", predLabels[2], "*", varLabels[1], ")")
}
term(c("a", "h"), "x")


## ----sprintf-------------------------------------------------------------
term = function(predLabels, varLabels){
    sprintf("%s * %s / (1 + %s * %s * %s)",
            predLabels[1], varLabels[1],
            predLabels[1], predLabels[2], varLabels[1])
}


## ----nonlin--------------------------------------------------------------
TypeII <- function(x){
  list(predictors = list(a = 1, h = 1),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels){
           sprintf("%s * %s / (1 + %s * %s * %s)",
                   predLabels[1], varLabels[1],
                   predLabels[1], predLabels[2], varLabels[1])
})
}
class(TypeII) <- "nonlin"


## ----prey----------------------------------------------------------------
Density <- rep(c(2,5,10,15,20,30), each = 4)
Eaten <- c(1,1,0,0,2,2,1,1,1,2,3,2,2,2,3,3,3,3,4,3,3,3,4,3)


## ----mod1----------------------------------------------------------------
mod1 <- gnm(Eaten ~ -1 + TypeII(Density), start = c(a = 0.1, h = 0.1),
            family = quasipoisson(link = "identity"))


## ----mod1Summary, echo = FALSE-------------------------------------------
summary(mod1)


## ----factor--------------------------------------------------------------
TypeII <- function(C, x){
  list(predictors = list(a = substitute(C), h = substitute(C)),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels){
           sprintf("%s * %s / (1 + %s * %s * %s)",
                   predLabels[1], varLabels[1],
                   predLabels[1], predLabels[2], varLabels[1])
})
}
class(TypeII) <- "nonlin"


## ----factorResult--------------------------------------------------------
Catchment <- factor(rep(1:2, 6, each = 2))
mod2 <- gnm(Eaten ~ -1 + TypeII(Catchment, Density),
            start = rep(0.2, 4),
            family = quasipoisson(link = "identity"))
coef(mod2)


## ----formula-------------------------------------------------------------
TypeII <- function(f, x){
  list(predictors = list(a = f, h = f),
       variables = list(substitute(x)),
       term = function(predLabels, varLabels){
           sprintf("(%s) * (%s)/ (1 + (%s) * (%s) * %s)",
                   predLabels[1], varLabels[1],
                   predLabels[1], predLabels[2], varLabels[1])
})
}
class(TypeII) <- "nonlin"


## ----formulaResult-------------------------------------------------------
mod2 <- gnm(Eaten ~ -1 + TypeII(~ 1 + Catchment, Density),
            start = c(0.2, -0.1, 0.2, -0.1),
            family = quasipoisson(link = "identity"))
coef(mod2)


## ----binomial, eval = FALSE----------------------------------------------
## count <- with(voting, percentage/100 * total)
## yvar <- cbind(count, voting$total - count)


## ----upward, eval = FALSE------------------------------------------------
## origin <- as.numeric(as.character(voting$origin))
## destination <- as.numeric(as.character(voting$destination))
## upward <- origin > destination


## ----inOut, eval = FALSE-------------------------------------------------
## in1 <- origin != 1 & destination == 1
## out1 <- origin == 1 & destination != 1

