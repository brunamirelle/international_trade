#Removes previously loaded objects
rm(list=ls())

#Sets working directory
setwd("~/Área de Trabalho/trade_projeto_1")

#Downloads and installs required packages
if(!require("xtable"))
{
  install.packages("xtable")
  library("xtable",character.only = T)
}

#theta de Costinot 
theta = 6.534

#Loads exporter-importer-industry  indices 
indices = read.csv("labelbook.csv", stringsAsFactors = F)

#Loads data on trade shares, consumption shares and employment
shares = read.csv("shares.csv", stringsAsFactors = F)

#Merges datasets
base = merge(indices, shares, by.x = c("expor_label","impor_label","indus_label"),
             by.y = c("expor", "impor", "indus"))

#Corrects exporter label -- ROU not ROM for Roumania
base$impor_label[base$impor_label=="ROM"] = "ROU"

#Works with millions of workers
base$L_i = base$L_i/10^(3)

#Retorna matriz A, tal que Aw = 0 é a condição de trade balance, onde w é um vetor com o salário
#de cada país. Argumentos da função são
#dados - Base de dados com trade share, consumption shares, e emprego em cada país
trade_balance <- function(dados)
{
  A = matrix(0, nrow = length(unique(dados$expor_label)), ncol = length(unique(dados$impor_label)))
  expor_nomes = unique(dados$expor_label)
  impor_nomes = unique(dados$impor_label)
  for(e_num in 1:length(expor_nomes))
    for(i_num in 1:length(impor_nomes))
    {
      ee = expor_nomes[e_num]
      ii = impor_nomes[i_num]
      
      base_rest = dados[dados$expor_label==ee&dados$impor_label==ii,]
      L_j = (dados$L_i[dados$expor_label==ii])[1]
      if(ee!=ii)
        A[e_num,i_num] = sum(base_rest$alpha_jk*base_rest$pi_ijk)*L_j else
          A[e_num,i_num] = (sum(base_rest$alpha_jk*base_rest$pi_ijk)-1)*L_j
    }

  colnames(A) = impor_nomes
  rownames(A) = expor_nomes
  return(A)
}

#Calcula sistema, supondo que trade shares gerados pelo modelo são iguais aos dados
A_dados = trade_balance(base)

#O vetor de salários que satisfaz trade balance é o autovetor associado ao autovalor 0 da matriz
salarios =  eigen(A_dados)$vectors[,which.min(abs(eigen(A_dados)$values))]

#Normaliza salário dos EUA para 1
salarios = salarios/salarios[40]

#Adiciona resultados
resultados = cbind("País" = unique(base$expor_label), "Salários" = round(salarios,digits=4))

#Adiciona salários dos exportadores na base
base_sal = merge(base, data.frame("expor_label"= unique(base$expor_label), "w_i" = salarios), by = "expor_label")

#Inclui produtivid_ade - computada do Stata
productivity = read.csv("zijk2.csv")
colnames(productivity)[3] = "zik"
productivity = productivity[order(productivity$exporter,productivity$industry),]

#Merges datasets
base_prod = merge(base_sal, productivity, by.x = c("indus","expor"), by.y = c("industry","exporter"))
base_prod = base_prod[order(base_prod$expor, base_prod$indus, base_prod$impor),]

#Calcula dijk segundo fórmula no trabalho

d = c()

for(x in 1:nrow(base_prod))
{
   print(x)
   i = base_prod$expor_label[x]
   j = base_prod$impor_label[x]
   k = base_prod$indus_label[x]
   pi_ijk = base_prod$pi_ijk[x]
   pi_jjk = base_prod$pi_ijk[base_prod$expor_label==j&base_prod$impor_label==j&base_prod$indus_label==k]
   w_i = base_prod$w_i[x]
   w_j = base_prod$w_i[base_prod$expor_label==j&base_prod$impor_label==j&base_prod$indus_label==k]
   z_ik = base_prod$zik[x]
   z_jk =  base_prod$zik[base_prod$expor_label==j&base_prod$impor_label==j&base_prod$indus_label==k]
 
   d = c(d,(pi_ijk/pi_jjk)^(-1/theta)*((w_j*z_jk)/(w_i*z_ik)))
}

#Adiciona barreiras a base, corrigindo "problemas"
base_prod = cbind(base_prod, "dijk" = d)
base_prod$dijk[is.nan(base_prod$dijk)&base$expor_label==base$impor_label] = 1
base_prod$dijk[is.nan(base_prod$dijk)&base$expor_label!=base$impor_label] = Inf 
base_prod$dijk[base_prod$dijk<=1]=1

#Salva resultados da calibragem
write.csv(base_prod,"calibragem.csv")

modelo_base = base_prod 


#Dado um vetor w (com a mesma ordenação da base dados passada) IxIxK, computa o trade share previsto pelo modelo,
#dados iceberg costs, produtividades e shares de consumo na base dados (ordenada)
trade_share_previsto <- function(w, dados)
{
  pi_pred = (w*dados$dijk/dados$zik)^(-theta)
  for(j in unique(dados$impor_label))
    for(k in unique(dados$impor))
    pi_pred[dados$impor_label==j&dados$indus_label==k] =  pi_pred[dados$impor_label==j&dados$indus_label==k]/sum(pi_pred[dados$impor_label==j&dados$indus_label==k])
  return(pi_pred)
}
#Função que dado um vetor (I-1) de salários (ordenados conforme a base) e uma base
#com os parâmetros, calcula a distância das equações do modelo em relação ao zero. O salário
#do último país é normalizado para 0
discrep_modelo <- function(w, dados)
{
  base_merge = data.frame("expor_label" = unique(dados$expor_label)[order(unique(dados$expor_label))], "w_iteracao" = c(w,1))
  dados_merged = merge(dados, base_merge, by = "expor_label")
  dados_merged = dados_merged[order(dados_merged$expor_label,dados_merged$indus_label,
                                    dados_merged$impor_label),]
  pi_ijk = trade_share_previsto(dados_merged$w_iteracao, dados_merged)
  dados_merged$pi_ijk = pi_ijk
  A = trade_balance(dados_merged)
  err = sum((A%*%c(w,1))^2)/10^20
  print(w)
  print(err)
  return(err)
}


####Contrafactuais####
#Cria base de contrafactuais
#Contrafactual_1 - reduzir custos de importação ao Brasil a 50% de seu valor
contra1 = modelo_base[, !(colnames(modelo_base)%in%c("w_i", "pi_ijk"))]
contra1$dijk[contra1$impor_label=="BRA"] = contra1$dijk[contra1$impor_label=="BRA"]/2

w_contra1 = optim(fn = discrep_modelo, par = salarios[-length(salarios)], method = "L-BFGS-B", dados = contra1)

