# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 20:49:13 2020

@author: murad
"""

import numpy as np

def der (mat_pams, sig):
   # second derivatives of f
    df = np.zeros((6,6))

    #Material Parameter (0=E; 1=mue; 2=H; 3=M; 4=a; 5=b; 6=h; 7=G)
    M = mat_pams[3]
    a = mat_pams[4]
    b = mat_pams[5]
    h = mat_pams[6]
    sig11 = sig[0]
    sig22 = sig[1]
    sig12 = sig[5]
    
    # DIESE ABLIETUNGEN SIND NUR FÜR GERADE M KORREKT !!!!!!!
    
    df[0,0] = M*(2**M*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + 2**M*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + 2**M*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(-1/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - (-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + (-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2))/2    
    
    df[0,1] = M*(-2**M*h*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) - 2**M*h*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + 2**M*h*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(-h/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + h*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - h*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2))/2    
    
    df[0,5] = M*(2**M*b**2*sig12*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2) - 2**M*b**2*sig12*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) - a*b**2*sig12*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) - a*b**2*sig12*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(-h*sig22 + sig11)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) - a*b**2*sig12*(-h*sig22 + sig11)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))/2    
    
    
    
    df[1,0] = M*h*(-2**M*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) - 2**M*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - 2**M*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(-1/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - (-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + (-h*sig22 + sig11)*(h*sig22/4 - sig11/4)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2))/2    
    
    df[1,1] = M*h*(2**M*h*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + 2**M*h*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - 2**M*h*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(-h/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) + h*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)) - h*(-h*sig22 + sig11)*(-h*sig22/2 + sig11/2)/(4*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2))/2
    
    df[1,5] = M*h*(-2**M*b**2*sig12*(2 - a)*(M/2 - 1/2)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2) + 2**M*b**2*sig12*(2 - a)*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(M/2 - 1/2)*(-h*sig22 + sig11)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*b**2*sig12*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(1 - (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) - a*b**2*sig12*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) - a*b**2*sig12*(1 + (-h*sig22 + sig11)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) - a*b**2*sig12*(-h*sig22 + sig11)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + a*b**2*sig12*(-h*sig22 + sig11)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)))/2




    df[5,0] = M*b**2*sig12*(-a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) - a*(1/2 - (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(1/2 + (-h*sig22/4 + sig11/4)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + (2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a)*(M - 1)*(-h*sig22/4 + sig11/4)/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + M*b**2*sig12*(h*sig22/4 - sig11/4)*(-a*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + (2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a))/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)    
    
    df[5,1] = M*b**2*h*sig12*(-h*sig22/2 + sig11/2)*(-a*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + (2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a))/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2)) + M*b**2*sig12*(a*(M - 2)*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) - a*(M - 2)*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h/2 - h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) - a*(h/2 + h*(-h*sig22/2 + sig11/2)/(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) - h*(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a)*(M - 1)*(-h*sig22/2 + sig11/2)/(2*(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)))/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)
    
    df[5,5] = -M*b**4*sig12**2*(-a*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + (2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a))/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)**(3/2) + M*b**2*sig12*(a*b**2*sig12*(M - 2)*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(M - 2)*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + a*b**2*sig12*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2)/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + b**2*sig12*(2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a)*(M - 1)/(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2) + M*b**2*(-a*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 - np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + a*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))*(h*sig22/2 + sig11/2 + np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 2) + (2*np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2))**(M - 1)*(4 - 2*a))/np.sqrt(b**2*sig12**2 + (-h*sig22/2 + sig11/2)**2)

    
    return df