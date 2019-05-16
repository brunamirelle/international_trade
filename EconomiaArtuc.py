#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  2 21:47:08 2019

@author: brunamirelle
"""

import numpy as np
from scipy.optimize import fsolve

#Resolve o modelo de Artuc (2008)
class EconomiaArtuc:
    
    """
    Inicializa economia
    alpha - parâmetro da função de produção (Cobb-Douglas)
    beta - fator de desconto dos trabalhadores
    C - custo de mudar de setor
    K_X e K_Y = dotações de capital em cada setor
    L_minus = dotação de trabalho na economia
    p0 = preço pré-liberalização do bem X
    """
    def __init__(self, alpha, beta, C, nu, K_X, K_Y, L_minus, p0):
        self.alpha = alpha
        self.beta = beta
        self.C = C
        self.nu = nu
        self.K_X = K_X
        self.K_Y = K_Y
        self.L_minus = L_minus
        self.p0 = p0
    
    #Função Omega no paper
    def Omega(self, mu):
        return self.nu*np.log(1 + np.exp(mu/self.nu))
    
    #CDF da distribuição logística (diferença dos choques)
    def G(self, mu):
        return np.exp(mu/self.nu)/(1 + np.exp(mu/self.nu))
    
    #Condição de primeira ordem para maximização dos lucros
    def FOC(self, L, K, p):
        return p*self.alpha*((L)**(self.alpha-1))*((K)**(1-self.alpha))
    
    #Deflator
    def phi(self,p):
        return p**(.5)
    
    #Retorna as equações necessárias para calcular steady-state
    #vec-guess - valores, salários e thresholds (vetor)
    #p - preço relativo do bem X no steady-state desejado
    def steady_discrep(self,vec_guess, p):
        V_X = vec_guess[0]
        V_Y = vec_guess[1]
        mu_X = vec_guess[2]
        mu_Y = vec_guess[3]
        L_X = vec_guess[4]
        L_Y = vec_guess[5]
        w_X = vec_guess[6]
        w_Y = vec_guess[7]
        discrep = np.zeros(8)
        discrep[0] = w_X - self.FOC(L_X, self.K_X, p)/self.phi(p)
        discrep[1] = w_Y - self.FOC(L_Y, self.K_Y, 1)/self.phi(p)
        discrep[2] = L_X + L_Y - self.L_minus
        discrep[3] = mu_X - self.beta*(V_Y-V_X) + self.C
        discrep[4] = mu_Y - self.beta*(V_X-V_Y) + self.C
        discrep[5] = V_X - w_X - self.beta*(V_X) - self.Omega(mu_X)
        discrep[6] = V_Y - w_Y - self.beta*(V_Y) - self.Omega(mu_Y)
        discrep[7] = L_X - (1-self.G(mu_X))*L_X - self.G(mu_Y)*L_Y
        #discrep[8] = L_Y - (1-self.G(mu_Y))*L_Y - self.G(mu_X)*L_X
        return discrep  
    
    #Calcula estado estacionário para um dado valor do preço externo
    def solve_steady(self, p):
        def steady_price(vec_guess):
            return self.steady_discrep(vec_guess, p)
        return fsolve(steady_price, np.ones(8))
    
    """
    Resolve o modelo para um dado valor pós-liberalização
    p_after = preço relativo de X pós-liberalização
    T_SS = número de períodos até novo steady-state (padrão 30, t=31 é steady-state)
    effective_at = a partir de que período vale a abertura (padrão é o zero)?
    tol = Tolerância para convergência do algoritmo
    """
    def solve_model(self, p_after, T_SS = 30, effective_at = 0, tol = 1e-5):
        #Calcula valores do steady-state inicial e final
        #Inicial
        V_X_0, V_Y_0, mu_X_0, mu_Y_0, L_X_0, L_Y_0, w_X_0, w_Y_0 = self.solve_steady(self.p0)
        #Final  
        V_X_T, V_Y_T, mu_X_T, mu_Y_T, L_X_T, L_Y_T, w_X_T, w_Y_T = self.solve_steady(p_after)
        
        
        #Chutes iniciais para a trajetória dos valores em cada setor
        z_X_t = np.linspace(V_X_0, V_X_T, T_SS + 1)
        z_Y_t = np.linspace(V_Y_0, V_Y_T, T_SS + 1)
        
        
        #Trajetória de preços
        p_t = np.concatenate((np.repeat(self.p0,effective_at),np.repeat(p_after,T_SS + 1 - effective_at)))
        
        #print(p_t)
        
        err = tol + 1
        
        while err > tol:
            #Calcula trajetórias para mu_X e mu_Y
            mu_X_t = self.beta*(np.concatenate((z_Y_t[1:],[V_Y_T])) -\
                                np.concatenate((z_X_t[1:],[V_X_T])) ) - self.C
            mu_Y_t = self.beta*(np.concatenate((z_X_t[1:],[V_X_T])) -\
                                np.concatenate((z_Y_t[1:],[V_Y_T])) ) - self.C
           
            L_X_t = np.empty(T_SS + 1)
            L_Y_t = np.empty(T_SS + 1)
            
            L_X_t[0] = L_X_0
            L_Y_t[0] = L_Y_0
            
            for tt in range(1, T_SS + 1):
                L_X_t[tt] = (1-self.G(mu_X_t[tt-1]))*L_X_t[tt-1] +\
                self.G(mu_Y_t[tt-1])*L_Y_t[tt-1]
                L_Y_t[tt] = self.L_minus - L_X_t[tt]
                
            w_Y_t = self.FOC(L_Y_t, self.K_Y, 1)/self.phi(p_t)
            w_X_t = self.FOC(L_X_t, self.K_X, p_t)/self.phi(p_t)
            
            z_tilde_X_t = w_X_t + self.beta*np.concatenate((z_X_t[1:], [V_X_T])) + self.Omega(mu_X_t)
            z_tilde_Y_t = w_Y_t + self.beta*np.concatenate((z_Y_t[1:], [V_Y_T])) + self.Omega(mu_Y_t)
            
            err = np.max((np.abs(z_tilde_X_t - z_X_t), np.abs(z_tilde_Y_t - z_Y_t)))
            
            z_X_t = z_tilde_X_t 
            z_Y_t = z_tilde_Y_t
            
            print(err)
        
        return [np.concatenate(([w_X_0],w_X_t,[w_X_T])),\
                np.concatenate(([w_Y_0],w_Y_t,[w_Y_T])),\
                np.concatenate(([V_X_0],z_X_t,[V_X_T])),\
                np.concatenate(([V_Y_0],z_Y_t,[V_Y_T])),\
                np.concatenate(([L_X_0],L_X_t,[L_X_T])),\
                np.concatenate(([L_Y_0],L_Y_t,[L_Y_T]))]

                                
    
            
        
    
        
        
    