#Removes previously loaded objects
rm(list=ls())

#Sets working directory
setwd("~/Área de Trabalho/trade_projeto_1")

#Downloads and installs required packages
for(pack in c("xtable","pracma"))
if(!require(pack,character.only = T))
{
  install.packages(pack)
  library(pack,character.only = T)
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
resultados = cbind("País" = unique(base$expor_label), "Salários (calibragem)" = round(salarios,digits=4))

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
#Se infinito, vamos colocar 10^3 vezes o maior custo calibrado
base_prod$dijk[is.infinite(base_prod$dijk)] = NA
base_prod$dijk[is.na(base_prod$dijk)] = 10^3*max(base_prod$dijk,na.rm = T)

#Salva resultados da calibragem
write.csv(base_prod,"calibragem.csv")

modelo_base = base_prod 

####Resolvendo Modelo####
###Nós vamos resolver o modelo sob diferentes contrafactuais#### 
###Para isso, vamos trabalhar com arrays 3D###
#As dimensões em que trabalharemos são (exportador, produto, importador) = (i,k,j)

#Ordena base conforme essas dimensões
modelo_base = modelo_base[order(modelo_base$impor_label, modelo_base$indus_label, modelo_base$expor_label),]

#Caso base -- resolvemos modelo, agora para os parâmetros computados -- para ver se dá mt diferente
#usamos esses como base de comparação
modelo_caso = modelo_base

#Dimensões 
dim_exp = length(unique(modelo_caso$expor_label))
dim_indus = length(unique(modelo_caso$indus_label))
dim_imp = length(unique(modelo_caso$impor_label))
#Array com alpha_jk
alpha_jk = array(modelo_caso$alpha_jk, c(dim_exp, dim_indus, dim_imp)) 
#Array com alpha_ik
alpha_ik = array(apply(alpha_jk, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#Array com d_ijk
d_ijk = array(modelo_caso$dijk, c(dim_exp, dim_indus, dim_imp))
#Array com L_i
L_i = array(modelo_caso$L_i, c(dim_exp, dim_indus, dim_imp))
#Array com L_j
L_j =  array(apply(L_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#z_ik
z_ik = array(modelo_caso$zik, c(dim_exp,dim_indus, dim_imp))
#z_jk
z_jk = array(apply(z_ik, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#Função para calcular discrepância do modelo. Toma como argumento um vetor I-1
#de salários
resolve_modelo <- function(w)
{
  #Inclui último salário normalizado igual a 1
  w = c(w,1)
  
  #Matriz w_i
  w_i = array(w, c(dim_exp, dim_indus, dim_imp))
  
  #Matriz w_j
  w_j = array(apply(w_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  #Calula pi_ijk
  pi_ijk = ((w_i*d_ijk)/z_ik)^(-theta)
  norm_pi_ijk = array(apply(array(t(colSums(pi_ijk,dims=1)),c(dim_imp, dim_indus, dim_imp)),1,function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  pi_ijk = pi_ijk/norm_pi_ijk
  
  #Calcula condição de trade balance
  eq = (rowSums(pi_ijk*alpha_jk*w_j*L_j) - w*L_i[,1,1])
  #Remove última linha da equação (redundante)
  eq = eq[-length(eq)]
  print(max(abs(eq)))
  return(eq)
}

#Resolve modelo. Chute inicial é o equilíbrio do modelo anterior
base_sol = fsolve(resolve_modelo, x0 = salarios[-length(salarios)])

salarios = c(base_sol$x,1) 

resultados = cbind(resultados, "Salários (modelo no caso base)" = round(salarios,digits=4))


#Contra - 1 - reduzir barreiras brasileiras a 50%
modelo_caso = modelo_base
modelo_caso$dijk[modelo_caso$impor_label=="BRA"] = modelo_caso$dijk[modelo_caso$impor_label=="BRA"]/2
#Se barreiras deram menos que 1, forçá-las a serem 1
modelo_caso$dijk[modelo_caso$dijk<=1]=1

#Dimensões 
dim_exp = length(unique(modelo_caso$expor_label))
dim_indus = length(unique(modelo_caso$indus_label))
dim_imp = length(unique(modelo_caso$impor_label))
#Array com alpha_jk
alpha_jk = array(modelo_caso$alpha_jk, c(dim_exp, dim_indus, dim_imp)) 
#Array com alpha_ik
alpha_ik = array(apply(alpha_jk, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#Array com d_ijk
d_ijk = array(modelo_caso$dijk, c(dim_exp, dim_indus, dim_imp))
#Array com L_i
L_i = array(modelo_caso$L_i, c(dim_exp, dim_indus, dim_imp))
#Array com L_j
L_j =  array(apply(L_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#z_ik
z_ik = array(modelo_caso$zik, c(dim_exp,dim_indus, dim_imp))
#z_jk
z_jk = array(apply(z_ik, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#Função para calcular discrepância do modelo. Toma como argumento um vetor I-1
#de salários
resolve_modelo <- function(w)
{
  #Inclui último salário normalizado igual a 1
  w = c(w,1)
  
  #Matriz w_i
  w_i = array(w, c(dim_exp, dim_indus, dim_imp))
  
  #Matriz w_j
  w_j = array(apply(w_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  #Calula pi_ijk
  pi_ijk = ((w_i*d_ijk)/z_ik)^(-theta)
  norm_pi_ijk = array(apply(array(t(colSums(pi_ijk,dims=1)),c(dim_imp, dim_indus, dim_imp)),1,function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  pi_ijk = pi_ijk/norm_pi_ijk
  
  #Calcula condição de trade balance
  eq = (rowSums(pi_ijk*alpha_jk*w_j*L_j) - w*L_i[,1,1])
  #Remove última linha da equação (redundante)
  eq = eq[-length(eq)]
  print(max(abs(eq)))
  return(eq)
}

#Resolve modelo. Chute inicial é o equilíbrio do modelo anterior
contra_sol = fsolve(resolve_modelo, x0 = salarios[-length(salarios)])

salarios_contra = c(contra_sol$x,1)
wi_tilde = array(salarios_contra, c(dim_exp, dim_indus, dim_imp))
#wj_tilde = array(apply(wi_tilde, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#Vamos calcular a mudança no índice de preços
dijk_base = array(modelo_base$dijk, c(dim_exp, dim_indus, dim_imp))
zik_base = array(modelo_base$zik, c(dim_exp,dim_indus, dim_imp))
#zjk_base = array(apply(zik_base, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
wi_base = array(salarios, c(dim_exp,dim_indus, dim_imp))
#wj_base = array(apply(wi_base, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#tilde_pik/p_ik
pjk_change = (colSums(((wi_tilde*d_ijk)/z_ik)^(-theta),dims=1)/colSums(((wi_base*dijk_base)/zik_base)^(-theta),dims=1))^(-1/theta)

p_i_change = exp(rowSums(log(t(pjk_change))*alpha_ik[,,1]))

welfare_change = 100*((salarios_contra/salarios)*(1/p_i_change)-1)

resultados = cbind(resultados, "Contrafactual 1 - Ganho de bem-estar (%)" = round(welfare_change,digits=4))

#Contra - 2 - aumentar produtividade brasileira em 50%
modelo_caso = modelo_base
modelo_caso$zik[modelo_caso$expor_label=="BRA"] = modelo_caso$zik[modelo_caso$expor_label=="BRA"]*2

#Dimensões 
dim_exp = length(unique(modelo_caso$expor_label))
dim_indus = length(unique(modelo_caso$indus_label))
dim_imp = length(unique(modelo_caso$impor_label))
#Array com alpha_jk
alpha_jk = array(modelo_caso$alpha_jk, c(dim_exp, dim_indus, dim_imp)) 
#Array com alpha_ik
alpha_ik = array(apply(alpha_jk, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#Array com d_ijk
d_ijk = array(modelo_caso$dijk, c(dim_exp, dim_indus, dim_imp))
#Array com L_i
L_i = array(modelo_caso$L_i, c(dim_exp, dim_indus, dim_imp))
#Array com L_j
L_j =  array(apply(L_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
#z_ik
z_ik = array(modelo_caso$zik, c(dim_exp,dim_indus, dim_imp))
#z_jk
z_jk = array(apply(z_ik, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#Função para calcular discrepância do modelo. Toma como argumento um vetor I-1
#de salários
resolve_modelo <- function(w)
{
  #Inclui último salário normalizado igual a 1
  w = c(w,1)
  
  #Matriz w_i
  w_i = array(w, c(dim_exp, dim_indus, dim_imp))
  
  #Matriz w_j
  w_j = array(apply(w_i, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  #Calula pi_ijk
  pi_ijk = ((w_i*d_ijk)/z_ik)^(-theta)
  norm_pi_ijk = array(apply(array(t(colSums(pi_ijk,dims=1)),c(dim_imp, dim_indus, dim_imp)),1,function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
  pi_ijk = pi_ijk/norm_pi_ijk
  
  #Calcula condição de trade balance
  eq = (rowSums(pi_ijk*alpha_jk*w_j*L_j) - w*L_i[,1,1])
  #Remove última linha da equação (redundante)
  eq = eq[-length(eq)]
  print(max(abs(eq)))
  return(eq)
}

#Resolve modelo. Chute inicial é o equilíbrio do modelo anterior
contra_sol = fsolve(resolve_modelo, x0 = contra_sol$x)

salarios_contra = c(contra_sol$x,1)
wi_tilde = array(salarios_contra, c(dim_exp, dim_indus, dim_imp))
#wj_tilde = array(apply(wi_tilde, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#Vamos calcular a mudança no índice de preços
dijk_base = array(modelo_base$dijk, c(dim_exp, dim_indus, dim_imp))
zik_base = array(modelo_base$zik, c(dim_exp,dim_indus, dim_imp))
#zjk_base = array(apply(zik_base, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))
wi_base = array(salarios, c(dim_exp,dim_indus, dim_imp))
#wj_base = array(apply(wi_base, 1, function(x){t(x)}),c(dim_exp, dim_indus, dim_imp))

#tilde_pik/p_ik
pjk_change = (colSums(((wi_tilde*d_ijk)/z_ik)^(-theta),dims=1)/colSums(((wi_base*dijk_base)/zik_base)^(-theta),dims=1))^(-1/theta)

p_i_change = exp(rowSums(log(t(pjk_change))*alpha_ik[,,1]))

welfare_change = 100*((salarios_contra/salarios)*(1/p_i_change)-1)

resultados = cbind(resultados, "Contrafactual 2 - Ganho de bem-estar (%)" = round(welfare_change,digits=4))

print.xtable(xtable(resultados, caption = "Resultados do exercício"), file = "resultados.tex")
