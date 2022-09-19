import numpy as np
import pandas as pd
import torch
from math import factorial, log
     
#produces nicely-spaced s/tau_star values
class iSITH(torch.nn.Module):
    def __init__(self, tau_min=.1, tau_max=100., buff_max=None, k=50, ntau=50, dt=1, g=0.0,
                 ttype=torch.FloatTensor):
        super(iSITH, self).__init__()
        self.k = k
        self.tau_min = tau_min
        self.tau_max = tau_max
        if buff_max is None:
            buff_max = 3*tau_max
        self.buff_max = buff_max
        self.ntau = ntau
        self.dt = dt
        self.g = g

        self.c = (tau_max/tau_min)**(1./(ntau-1))-1
        self.tau_star = tau_min*(1+self.c)**torch.arange(ntau).type(torch.DoubleTensor)
        self.s = 1/self.tau_star

#Associative Memory Matrix M(s,a)
class mem_matrix(torch.nn.Module):
    def __init__(self, nf=2, S=torch.FloatTensor
                 #, A=torch.FloatTensor
                ):
        super(mem_matrix, self).__init__()
        #self.A = A
        self.S = S
        self.M = torch.FloatTensor().new_zeros((len(S)
                                                #,len(A)
                                                ,nf,nf))
    def time_update(a, f_IN: np.array, ctx: torch.FloatTensor):
        for s_idx, s in enumerate(self.S):
            #sum all predictions made by other stimuli
            predicted = np.dot(M[s_idx], (1-x for x in f_IN))
            #calculate observed stimuli subtracted by predicted stimuli
            f_bind = f_IN + s * np.dot(self.M[s_idx],f_IN) - predicted
            #bind observation and context, with learning factor a
            self.M[idx] += a * np.outer(f_bind, ctx)
            #set diagonal to zeros to prevent self-prediction
            self.M[idx].fill_diagonal_(0)
            
#Context F(s)
class context(torch.nn.Module):
    def __init__(self, nf=2, S=torch.FloatTensor):
        super(context, self).__init__()
        self.S = S
        self.F_pre = torch.FloatTensor().new_zeros(len(self.S), nf)
        self.F_post = torch.FloatTensor().new_zeros(len(self.S), nf)
    def time_update(f_IN: np.array):
        for idx, s in enumerate(self.S):
            self.F_pre = self.F_post
            self.F_post += - s*self.F_post + f_IN
            
#Stimulus Predictor P(s)          
class transition(torch.nn.Module):
    def __init__(self, nf=2, S=torch.FloatTensor):
        super(transition, self).__init__()
        self.S = S
        self.P = torch.FloatTensor().new_zeros(len(self.S), nf)
    def prob_input(f_IN: np.array, M=torch.FloatTensor):
        for idx, s in enumerate(self.S):
            self.P_IN[idx] = np.dot(M[idx], f_IN)
    def time_update():
        for idx, s in enumerate(self.S):
            self.P[idx] += s*self.P[idx] + self.P_IN[idx] - self.P[-1]
    