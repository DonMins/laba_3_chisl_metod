import math as mp
from symtable import Symbol
from scipy.misc import derivative
import numpy as np
from scipy.optimize import fmin
import scipy.optimize as cs
from sympy import *

EPS = 10**-5
a=4
b=6

def inputFunc(x):
    return mp.sin(3*x) / (mp.sqrt(x) + 5 )
def dc(x):
    return -9*sin(3*x)/(sqrt(x) + 5) + sin(3*x)/(2*x*(sqrt(x) + 5)**3) - 3*cos(3*x)/(sqrt(x)*(sqrt(x) + 5)**2) + sin(3*x)/(4*x**(3/2)*(sqrt(x) + 5)**2)

def maxi(a,b):

    PHI = (3 - sqrt(5)) / 2
    while (abs(b - a) >= 1):
            x1 = a + (b - a)  * PHI
            x2 = b - (b - a) * PHI

            if (dc(x1) >= dc(x2)):
                a = x1
            else :
                if (dc(x1) < dc(x2)):
                 b = x2
    return (a + b) / 2

def max_D2x(x,b):
    max = 0
    while (x<= b):
        tmp = abs(derivative(inputFunc, x, dx=1e-6,n=2))
        if (max < tmp):
            max = tmp
        x+=0.00001
    return max

def step():
    M2= 1.21617
    h = sqrt(12*EPS/(M2*(b-a)))
    return h

def trapez():
    h = step()
    n = int((b-a)/h)
    print (n)
    I=0
    for i in range(1,n,1):
        I+=inputFunc(a+h*i)
    res = ((inputFunc(a)+inputFunc(b))/2 + I)*h
    return res

if __name__ == "__main__":
    print(maxi(a,b))
