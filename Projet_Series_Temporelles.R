# Romain Lefranc
# Maxime Denizan


######
# Q1 #
######


# Importation des données
setwd("/Users/Romain/Documents/Romain/ENSAE 2A 2020-2021/S2/Séries temporelles/Projet-Series-Temporelles")
getwd()
datafile <- "valeurs_mensuelles.csv"
data <- read.csv(datafile, sep=";")

data

# Packages séries temporelles
require(zoo)
require(tseries)
require(fUnitRoots)

typeof(data$Dates) # integer
dates_char <- as.character(data$Dates)
typeof(dates_char) # character
dates_char[1];dates_char[length(dates_char)]
dates <- as.yearmon(seq(from=2011+1/12,to=2019+2/12,by=1/12))

X_t <- zoo(data$Indice,order.by=dates)
typeof(data$X_t)
typeof(X_t)

# plot de la série
# chercher une tendance ou saisonalité
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = "Indice brut de la production industrielle") 
acf(X_t, main = "")

# plot du logarithme de la série
W_t = log(X_t)
plot(W_t, xlab = "Années", ylab = "Logarithme de l'indice") 
acf(W_t, main = "")



######
# Q2 #
######


# correction de la tendance

lt <- lm(W_t ~ dates) # régression du log de la série sur les dates
summary(lt)
dates

# résidus de la régression du logarithme de la série sur les dates
# série corigée
Z_t <- lt$residuals
acf(Z_t, main = ""); pacf(Z_t, main = "")

# Test racine unitaire sur le logarithme de la série
# Test de Dickey-Fuller
adf <- adfTest(W_t, lag=0, type="ct")
adf
# function test d'auto-corrélation
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients))

# L’absence d’autocorrélation des résidus est rejetéee au moins une fois (Q(4) à Q(9)), le test ADF avec aucun retard n’est donc pas valide. 
# Ajoutons des retards de ∆Xt jusqu’à ce que les résidus ne soient plus autocorrélée.
adfTest_valid <- function(series,kmax,type){
  k <- 0
  noautocorr <- 0
  while (noautocorr==0) {
    cat(paste0("ADF with ",k, " lags: residuals OK? "))
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1; cat("OK \n")}
    else cat("nope \n")
    k <- k + 1}
  return(adf)
  }
# Puis    
adf <- adfTest_valid(W_t,24,"ct")
# Il a fallu considérer 9 retards au test ADF pour supprimer l’autocorrélation des résidus.
adf

dspread <- diff(W_t,1)
plot(dspread)

adf <- adfTest_valid(dspread,24, type="nc")
# Il n’a pas été nécessaire d’inclure des retards dans le test ADF (test DF simple).
adf


x <- dspread
par(mfrow=c(1,2))
acf(x);pacf(x)



# Test de Phillips-Perron
pp.test(W_t)


######
# Q3 #
######


# plots avant et après la transformation de la série
plot(X_t, xlab = "Années", ylab = expression(X[t])) ; plot(Z_t, xlab = "Années", ylab = expression(Z[t]))

######
# Q4 #
######


# Découvrir les valeurs de p* et q*
acf(Z_t); pacf(Z_t)

# Test de validité

# fonction Qtests déjà fait

# Test de nullité des coefficients

signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}


######
# Q5 #
######


# p = 0 et q = 1
arima001 <- arima(Z_t,c(0,0,1))
Qtests(arima001$residuals, 24, 1)
signif(arima001)

# p = 0 et q = 2
arima002 <- arima(Z_t,c(0,0,2))
Qtests(arima002$residuals, 24, 2)
signif(arima002)
ma2 <- arima002

# p = 1 et q = 1
arima101 <- arima(Z_t,c(1,0,1))
Qtests(arima101$residuals, 24, 2)
signif(arima101)

# p = 1 et q = 0
arima100 <- arima(Z_t,c(1,0,0))
Qtests(arima100$residuals, 24, 1)
signif(arima100)
ar1 <- arima100

# p = 1 et q = 2
arima102 <- arima(Z_t,c(1,0,2))
Qtests(arima102$residuals, 24, 3)
signif(arima102)                    


# AIC et BIC
models <- c("ma2", "ar1"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))


######
# Q8 #
######

prediction <- predict(ma2,2)
prediction

#On r?cup?re le sigma? de notre bruit

var_prev_1 <- ma2$sigma2
var_prev_2 <- ma2$sigma2*(1+(ma2$coef[1] + ma2$coef[2])**2)

ma2$sigma2
  
#Formule th?orique test
Bound_sup <- rbind(prediction$pred[1] +1.96*sqrt(var_prev_1),prediction$pred[2] +1.96*sqrt(var_prev_2))
Bound_inf <- rbind(prediction$pred[1] -1.96*sqrt(var_prev_1),prediction$pred[2] -1.96*sqrt(var_prev_2))

#sigma2 nous donne : the MLE of the innovations variance.

date_pred <- c(2013+9/12, 2013+10/12) 
coef_reg <- lt$coefficients


dev.off()
plot(NULL,NULL,xlim=c(2011+7/12, 2013+10/12),ylim=c(400,800),
     xlab="Ann?es",ylab=expression(X[t]))
lines(X_t, type = "o", pch = 16)# les donn?es de base


polygon(x=c(date_pred, rev(date_pred)),
        y=c(c(exp(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1]), 
            exp(Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1])),
            rev(c(exp(Bound_sup[1] + coef_reg[2]*date_pred[1] + coef_reg[1]), 
                  exp(Bound_sup[2] + coef_reg[2]*date_pred[1] + coef_reg[1])))),col="grey",border=NA) #?a cr?e l'IC
lines(c(X_t,exp(prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1])) )
points(date_pred, exp(prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1]), col="red", pch=16) # predictions

# Anexe

dev.off()
par(mfrow = c(1,2))
plot(X_t, xlab = "Ann?es", ylab = expression(X[t])); acf(X_t, main = "")
plot(W_t, xlab = "Ann?es", ylab = expression(W[t])); acf(W_t, main = "")
plot(Z_t, xlab = "Ann?es", ylab = expression(Z[t])); acf(Z_t, main = "")

