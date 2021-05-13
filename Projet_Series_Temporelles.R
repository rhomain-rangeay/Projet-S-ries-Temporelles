# Projet séries temporelles ENSAE 2A
# Romain Lefranc
# Maxime Denizan


######
# Q1 #
######


# Importation des données et préparation de l'espace de travail
path_romain <- "/Users/Romain/Documents/Romain/ENSAE 2A 2020-2021/S2/Series temporelles/Projet de series temporelles"
#path_maxime <- " /Users/maximedenizan/Documents/GitHub/Projet de series temporelles "
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

# Formatage des données de la série X_t
typeof(data$Dates) # integer
dates_char <- as.character(data$Dates)
typeof(dates_char) # character
dates_char[1];dates_char[length(dates_char)]
dates <- as.yearmon(seq(from=2010+11/12,to=2019+3/12,by=1/12))
X_t <- zoo(data$Indice,order.by=dates)


# Plot de la série X_t et de son autocorrélogramme pour repérer une tendance ou une saisonnalité au premier regard
dev.off()
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = "Indice de la production industrielle") 
acf(X_t, main = "")

######
# Q2 #
######


# Correction de la tendance
lt <- lm(X_t ~ dates) # régression de la série sur les dates
summary(lt)

# Résidus de la régression de la série X_t sur les dates, puis affichage des corrélations totales et partielles
Z_t <- lt$residuals
acf(Z_t, main = ""); pacf(Z_t, main = "")

# Test racine unitaire sur la série : test de Dickey-Fuller augmenté avec constante et tendance
adf <- adfTest(X_t, lag=0, type="ct")
adf

# Fonction test d'auto-corrélation
Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

Qtests(adf@test$lm$residuals,24,length(adf@test$lm$coefficients))

# L’absence d’autocorrélation des résidus est rejetéee systématiquement (Q(4) à Q(24)), le test ADF avec aucun retard n’est donc pas valide. 
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


###### Étude de la série différenciée ######
diff_X_t <- diff(X_t,1) # Série différenciée

# Testons maintenant la racine unitaire pour la série différenciée diff_X_t dont la stationarité semble relativement plus probable. 
par(mfrow = c(1,2))
plot(diff_X_t, xlab = "Années", ylab = expression(Delta*X[t])) 
acf(diff_X_t, main = "")

# La représentation graphique de la série différenciée précédente semble d'abord montrer l’absence de constante et de tendance non nulle.
# Vérifions avec une régression :
summary(lm(diff_X_t ~ dates[-1])) #sans la première date car on a différencié la série
# p-valeur est largement supérieure à 0.05 (0.887). Donc on ne rejette pas le test de nullité du coefficient associé à dates dans la régression de la série différenciée sur les dates
# CONCLUSION : Il y a bien ni constante ni tendance significative.

# On peut donc tester la stationnarité de diff_X_t avec le test ADF dans le cas sans constante ni tendance, en vérifiant l’absence autocorréation des résidus
adf <- adfTest_valid(diff_X_t,24, type="nc")
# On ajoute 6 retards dans le test ADF (test DF simple) pour qu'il soit valide.
adf
# Le test rejette la racine unitaire (p-value<0.05), on dira donc que la série différenciée est ”stationnaire”. 
# diff_X_t est donc I(1).

######
# Q3 #
######

# plots avant et après la transformation de la série
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = expression(X[t])) ; plot(diff_X_t, xlab = "Années", ylab = expression(Delta*X[t]))


######
# Q4 #
######

# Découvrir les valeurs de p* et q*
acf(diff_X_t, main = ""); pacf(diff_X_t, main = "")

# Les sous-modèles possibles sont les ARIMA(p,d,q) tels que p <= p*, d <= d*, q <= q*.
# Dans l’absolu, on cherche un modèle :
# — bien ajusté : les coefficients estimés (notamment les coefficients des ordres AR et MA les plus élevés) sont statistiquement significatifs
# — valide : les résidus ne sont pas autocorrélés

# Fonction de test des significativités individuelles des coefficients
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))*2
  return(rbind(coef,se,pval))
}

# Fonction pour estimer un arima et en vérifier l’ajustement et la validité
modelchoice <- function(p,q,data=diff_X_t, k=24) {
  estim <- try(arima(data, c(p,0,q),optim.control=list(maxit=20000)))
  if (class(estim)=="try-error") return(c("p"=p,"q"=q,"arsignif"=NA,"masignif"=NA,"resnocorr"=NA, "ok"=NA))
  arsignif <- if (p==0) NA else signif(estim)[3,p]<=0.05
  masignif <- if (q==0) NA else signif(estim)[3,p+q]<=0.05
  resnocorr <- sum(Qtests(estim$residuals,24,length(estim$coef)-1)[,2]<=0.05,na.rm=T)==0 
  checks <- c(arsignif,masignif,resnocorr)
  ok <- as.numeric(sum(checks,na.rm=T)==(3-sum(is.na(checks))))
  return(c("p"=p,"q"=q,"arsignif"=arsignif,"masignif"=masignif,"resnocorr"=resnocorr,"ok"=ok))
  }

# Fonction pour estimer et vérifier tous les arima(p,q) avec p<=pmax et q<=max
armamodelchoice <- function(pmax,qmax) {
  pqs <- expand.grid(0:pmax,0:qmax) 
  t(apply(matrix(1:dim(pqs)[1]),1,function(row) {
    p <- pqs[row,1]; q <- pqs[row,2]
    cat(paste0("Computing ARMA(",p,",",q,") nn"))
    modelchoice(p,q)}))
  }

# On choisira p <= 3 et q <= 20 au regard des autocorrélations totales et partielles.
pmax = 3
qmax = 20

# On teste les différents modèles au regard des pmax et qmax identifiés
armamodels <- armamodelchoice(pmax,qmax) # estime tous les arima 
selec <- armamodels[armamodels[,"ok"]==1&!is.na(armamodels[,"ok"]),] #modèles bien ajustés et valides
selec
# Commentaire : On a 11 modèles bien ajustés et valides .

# On crée une liste des ordres p et q des modèles candidats
pqs <- apply(selec,1,function(row) list("p"=as.numeric(row[1]),"q"=as.numeric(row[2]))) 

# On renomme les éléments de la liste
names(pqs) <- paste0("arma(",selec[,1],",",selec[,2],")") 

# On crée une liste des modèles
models <- lapply(pqs, function(pq) arima(diff_X_t,c(pq[["p"]],0,pq[["q"]]))) 
vapply(models, FUN.VALUE=numeric(2), function(m) c("AIC"=AIC(m),"BIC"=BIC(m))) #calcule les AIC et BIC des modèles candidats

# L'ARIMA (0,1,18) minimise l'AIC 
arima0118 <- arima(X_t,c(0,1,18))

# L'ARIMA (0,1,1) minimise le BIC
arima011 <- arima(X_t,c(0,1,1))

#####################################
# Q8 : prédiction avec ARIMA(0,1,1) #
#####################################

prediction <- predict(arima011,2)
prediction

#On récupère le sigma de notre bruit
var_prev_1 <- arima011$sigma2 #sigma2 nous donne : la variance du bruit d'innovation.
var_prev_1
var_prev_2 <- arima011$sigma2*(1+(1-arima011$coef[1])**2)
var_prev_2

#Formule théorique test
Bound_sup <- rbind(prediction$pred[1] +1.96*sqrt(var_prev_1),prediction$pred[2] +1.96*sqrt(var_prev_2))
Bound_inf <- rbind(prediction$pred[1] -1.96*sqrt(var_prev_1),prediction$pred[2] -1.96*sqrt(var_prev_2))

date_pred <- c(2019+04/12, 2019+05/12) 
regression_diff_X_t <- lm(diff_X_t ~ dates[-1])
coef_reg <- regression_diff_X_t$coefficients

# Affichage des prévisions ainsi que l'intervalle de confiance à 95%
dev.off()
plot(NULL,NULL,xlim=c(2017+11/12, 2019+05/12),ylim=c(100,130),
     xlab="Années",ylab=expression(X[t]))
lines(X_t, type = "o", pch = 16)# les données de base

polygon(x=c(date_pred, rev(date_pred)),
        y=c(c(Bound_inf[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
            Bound_inf[2] + coef_reg[2]*date_pred[2] + coef_reg[1]),
            rev(c(Bound_sup[1] + coef_reg[2]*date_pred[1] + coef_reg[1], 
                  Bound_sup[2] + coef_reg[2]*date_pred[1] + coef_reg[1]))),col="grey",border=NA) #ça crée l'IC
lines(c(X_t,prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1]))
points(date_pred, prediction$pred + coef_reg[2]*(date_pred) + coef_reg[1], col="red", pch=16) # predictions


# Annexes
dev.off()
par(mfrow = c(1,2))
plot(X_t, xlab = "Années", ylab = expression(X[t])); acf(X_t, main = "")
plot(Z_t, xlab = "Années", ylab = expression(Z[t])); acf(Z_t, main = "")
plot(diff_X_t, xlab = "Années", ylab = expression(Delta*X[t])); acf(diff_X_t, main = "")


