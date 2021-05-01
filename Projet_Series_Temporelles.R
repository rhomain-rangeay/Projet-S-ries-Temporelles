# Romain Lefranc
# Maxime Denizan


######
# Q1 #
######


# Importation des données

path_romain <- "/Users/Romain/Documents/Romain/ENSAE 2A 2020-2021/S2/Series temporelles/Projet de series temporelles"
#path_maxime <- 

setwd(path_romain)
#setwd(path_maxime)

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
dates <- as.yearmon(seq(from=2010+11/12,to=2019+3/12,by=1/12))


X_t <- zoo(data$Indice,order.by=dates)
typeof(data$X_t)
typeof(X_t)

# plot de la série
# chercher une tendance ou saisonalité
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = "Indice brut de la production industrielle") 
acf(X_t, main = "")

######## QUESTION A ECLAIRCIR #######
# lien heteroscédasticité autocorrélogramme ? Notre série est-elle heteroscedastique ?

######
# Q2 #
######


# correction de la tendance

lt <- lm(X_t ~ dates) # régression du log de la série sur les dates
summary(lt)

# résidus de la régression du logarithme de la série sur les dates
# série corigée
Z_t <- lt$residuals
acf(Z_t, main = ""); pacf(Z_t, main = "")

##### A FAIRE ######
# Il faudrait l'ACF et le PACF de Diff
#acf(diff_X_t, main = ""); pacf(diff_X_t, main = "")
# Au regard des graphiques, on a l'intuition que la série est différenciée semble stationnaire. On va vérifier cette supposition avec des tests.


# Test racine unitaire sur le logarithme de la série
# Test de Dickey-Fuller
adf <- adfTest(X_t, lag=0, type="ct")
adf

# fonction test d'auto-corrélation
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals, 24,length(adf@test$lm$coefficients))

# L’absence d’autocorrélation des résidus est rejetéee au moins une fois (Q(4) à Q(9)), le test ADF avec aucun retard n’est donc pas valide. 
# Ajoutons des retards de ∆Xt jusqu’à ce que les résidus ne soient plus autocorrélés.
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
adf <- adfTest_valid(X_t,24,"ct")
# Il a fallu considérer 8 retards au test ADF pour supprimer l’autocorrélation des résidus.
adf

# La racine unitaire n’est pas rejetée à un seuil de 95% pour la série en niveau, la série est donc au moins I(1) (intégrée d'ordre 1). 

# Testons maintenant la racine unitaire pour la série différenciée dspread. 
diff_X_t <- diff(X_t,1)

# La représentation graphique précédente semble montrer l’absence de constante et de tendance non nulle. Vérifions avec une régression :
summary(lm(diff_X_t ~ dates[-1])) #sans la première date car on a différencié la série

# La p-valeur est largement supérieure à 0.05.
# donc on ne rejette pas le test de nullité du coefficient associé à dates dans la régression de la série différenciée sur les dates
# CONCLUSION : Il y a bien ni constante ni tendance significative.
data
adf <- adfTest_valid(diff_X_t,24, type="nc")
# On ajoute 7 retards dans le test ADF (test DF simple).
adf

# Le test rejette la racine unitaire (p-value<0.05), on dira donc que la série différenciée est ”stationnaire”. 
# diff_X_t est donc I(1).

# Test de Phillips-Perron
pp.test(diff_X_t)


######
# Q3 #
######


# plots avant et après la transformation de la série
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = expression(X[t])) ; plot(diff_X_t, xlab = "Années", ylab = expression(Delta*X[t]))

#####
# Q4 follasse #
#####
mat <- matrix(NA,nrow=19+1,ncol=3+1) #matrice vide à remplir
rownames(mat) <- paste0("p=",0:19) #renomme les lignes
colnames(mat) <- paste0("q=",0:3) #renomme les colonnes
AICs <- mat #matrice des AIC non remplie
BICs <- mat #matrice des BIC non remplie
pqs <- expand.grid(0:19,0:3) #toutes les combinaisons possibles de p et q
  for (row in 1:dim(pqs)[1]){ #boucle pour chaque (p,q)
    p <- pqs[row,1] # récupère p
    q <- pqs[row,2] # récupère q
    estim <- try(arima(diff_X_t,c(p,0,q),include.mean = F)) #tente d’estimer l’ARIMA 
    AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigne l’AIC 
    BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigne le BIC
}

AICs==min(AICs) #affiche le modèle minimisant l’AIC
AICs
#L'ARIMA (7,1,2) minimise l'AIC

BICs==min(BICs) #affiche le modèle minimisant le BIC
BICs #affiche les BICs
#L'ARIMA (0,1,1) minimise l'AIC


######
# Q4 #
######


# Découvrir les valeurs de p* et q*
acf(diff_X_t, main = ""); pacf(diff_X_t, main = "")

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


# p = 19 et q = 0
arima1910 <- arima(X_t,c(19,1,0))
Qtests(arima1910$residuals, 24, 19)
signif(arima1910)

# p = 7 et q = 2
arima712 <- arima(X_t,c(7,1,2))
Qtests(arima712$residuals, 24, 1)
signif(arima712)

# p = 0 et q = 2
#arima002 <- arima(diff_X_t,c(0,0,2))
#Qtests(arima002$residuals, 24, 2)
#signif(arima002)
#ma2 <- arima002

# p = 1 et q = 1
arima101 <- arima(diff_X_t,c(1,0,1))
Qtests(arima101$residuals, 24, 2)
signif(arima101)

# p = 1 et q = 0
arima100 <- arima(diff_X_t,c(1,0,0))
Qtests(arima100$residuals, 24, 1)
signif(arima100)
ar1 <- arima100

# p = 1 et q = 2
arima102 <- arima(diff_X_t,c(1,0,2))
Qtests(arima102$residuals, 24, 3)
signif(arima102)                    


# AIC et BIC
models <- c("ma2", "ar1"); names(models) <- models
apply(as.matrix(models),1, function(m) c("AIC"=AIC(get(m)), "BIC"=BIC(get(m))))


######
# Q8 #
######

prediction <- predict(arima1910,2)
prediction

#On récupère le sigma de notre bruit

var_prev_1 <- arima1910$sigma2
var_prev_1
var_prev_2 <- arima1910$sigma2*(1+(arima1910$coef[1] + arima1910$coef[2])**2)
var_prev_2

#Formule théorique test
Bound_sup <- rbind(prediction$pred[1] +1.96*sqrt(var_prev_1),prediction$pred[2] +1.96*sqrt(var_prev_2))
Bound_inf <- rbind(prediction$pred[1] -1.96*sqrt(var_prev_1),prediction$pred[2] -1.96*sqrt(var_prev_2))

#sigma2 nous donne : the MLE of the innovations variance.

date_pred <- c(2019+04/12, 2019+05/12) 
coef_reg <- lt$coefficients
date_pred
coef_reg[2]

dev.off()
plot(NULL,NULL,xlim=c(2017+11/12, 2019+05/12),ylim=c(100,240),
     xlab="Années",ylab=expression(X[t]))
lines(X_t, type = "o", pch = 16)# les données de base

Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1]
Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1]

c(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
            Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1])
rev(c(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
            Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1]))
c(X_t,prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1])
rev(date_pred)

polygon(x=c(date_pred, rev(date_pred)),
        y=c(c(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
            Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1]),
            rev(c(Bound_sup[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
                  Bound_sup[2] + coef_reg[2]*date_pred[1] + coef_reg[1]))),col="grey",border=NA) #ça crée l'IC
lines(c(X_t,prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1]))
points(date_pred, prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1], col="red", pch=16) # predictions


prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1]

# Anexe

dev.off()
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = expression(X[t])); acf(X_t, main = "")
plot(W_t, xlab = "Années", ylab = expression(W[t])); acf(W_t, main = "")
plot(Z_t, xlab = "Années", ylab = expression(Z[t])); acf(Z_t, main = "")

