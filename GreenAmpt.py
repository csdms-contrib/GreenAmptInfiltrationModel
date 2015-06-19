#      
# Copyright (c) 2015, Peishi Jiang

#--------------------------------------------------------
# class GreenAmpt
#
#    __init__()
#    Fp
#    f
#    F
#    __EqnF
#    F_f
#
#--------------------------------------------------------
# This model is based on the Green-Ampt method described
# in <Applied Hydrology>, pp110.
# 
# Ref:
# Chow, Ven T., David R. Maidment, and Larry W. Mays. 
# Applied hydrology. 1988.
#--------------------------------------------------------

from scipy.optimize import newton
from pprint import pprint
from json import dump
from math import log
import numpy as np

class GreenAmpt(object):
    '''
    Green-Ampt Cumulative Infiltration 
    '''
        
    def __init__(self, K, dt, theta_i, theta_s, psi, i):
        """
        Constructor
        """
        self.K       = K                 # hydraulic conductivity
        self.dt      = dt                # time resolution
        self.theta_i = theta_i           # initial water content
        self.theta_s = theta_s           # saturated water content
        self.dtheta  = theta_s - theta_i # the change in the moisture content
        self.psi     = psi               # wetting front soil suction head
        if type(i) == list:
            self.i = i                   # rainfall intensity
        else:
            self.i = [i]
    
    def Fp(self, i):
        """
        Cumulative infiltration at the ponding time tp
        """
        return self.K * self.psi * self.dtheta / (i - self.K) 
    
    def F(self, F_t, dt_t):
        """
        Solve Equation of Green-Ampt Cumulative Infiltration __EqnF
        """
        F_t_next = lambda F: self.__EqnF(F_t, dt_t, F)
        return newton(F_t_next, 3)
    
    def f(self, F):
        """
        Generate Green-Ampt Infiltration Rate at time t
        """
        if F == 0:
            return -9999
        else:
            return self.K * (self.psi * self.dtheta / F + 1)    
        
    def __EqnF(self, F_t, dt_t, F):
        """
        Equation of Green-Ampt Cumulative Infiltration after ponding
        F:  Green-Ampt Cumulative Infiltration variable 
        """
        return F - F_t - self.K*dt_t - self.psi*self.dtheta*log((self.psi*self.dtheta+F)/(self.psi*self.dtheta+F_t))
    
    def F_f(self):
        """
        Generate the time series of cumulative infiltration and infiltration rate
        given the time series of rainfall intensity i
        """
        t_len = len(self.i)
        F_all = []; f_all = []; t_all = []
        # initial
        F_all.append(0)
        f_all.append(-9999)
        t_all.append(0)
        for ind in range(1, t_len+1):
            i_t = self.i[ind-1]
            f_t = f_all[ind-1]
            F_t = F_all[ind-1]
            if abs(f_t) <= i_t:
                # ponding occurs throught interval
                F_t_next = self.F(F_t, self.dt)
                f_t_next = self.f(F_t_next)
            elif abs(f_t) > i_t:
                # no ponding at the beginning of the interval
                F_t_next_temp = F_t + i_t*self.dt
                f_t_next_temp = self.f(F_t_next_temp)
                if abs(f_t_next_temp) > i_t:
                    # no ponding throughout interval
                    f_t_next = f_t_next_temp
                    F_t_next = F_t_next_temp
                elif abs(f_t_next_temp) <= i_t:
                    # ponding occurs during interval
                    Fp_t = self.Fp(i_t)
                    print i_t
                    dt_p = (Fp_t - F_t) / i_t
                    F_t_next = self.F(Fp_t, self.dt - dt_p)
                    f_t_next = self.f(F_t_next)        
            F_all.append(F_t_next)
            f_all.append(f_t_next)
            t_all.append(self.dt*(ind))        
        return F_all, f_all, t_all