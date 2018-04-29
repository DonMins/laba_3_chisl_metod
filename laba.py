import math as mp
from sympy import *
import matplotlib.pyplot as plt
import numpy as np

EPS = 10**-5
PHI = (mp.sqrt(5)+1)/2

def inputFunc(x):
  # return mp.sin(3*x) / (mp.sqrt(x) + 5 )
   return mp.sin(2 * x) / (x**5 + x + 1)  # Таня

def inputFunc2(x,y):# задание 2
    return y/(1+x) - (0.5 * y**2)

def DxInputFunc2(x,y):
    return -y/((x+1)**2)


def Runge_Kuta3 ():
    x = 0
    y = 1
    h = 0.1
    n=0
    X = []
    Y= []
    while (n<3):
        k1 = h*inputFunc2(x, y)
        k2 = h*inputFunc2(x + h/2, y + k1/2)
        k3 = h*inputFunc2(x + h, y + 2*k2 - k1)
        x = x + h
        y = y + (k1 + 4*k2+k3)/6
        X.append(x)
        Y.append(y)
        n+=1
    return X, Y

def Runge_Kuta4():
    x = 0
    y = 1
    h = 0.1
    n = 0
    X = []
    Y = []
    while(n<3):

        k1 = h*inputFunc2(x,y)
        k2 = h*inputFunc2(x + h/3, y + k1/3)
        k3 = h*inputFunc2(x + 2*h/3, y - k1/3+k2)
        k4 = h*inputFunc2(x + h, y + k1+k3-k2)
        x= x + h
        y+=(k1 + 3*k2+3*k3+k4)/8
        n+=1
        X.append(x)
        Y.append(y)
    return X, Y

def prognoz():
    x = 0
    y = 1
    h = 0.1
    Y= []
    X=[]
    n=0
    while(n<3):
        yh2 = y+h**2*inputFunc2(x,y)
        y = y+h*inputFunc2(x,y)+1/2*(inputFunc2(x+h**2,yh2)-inputFunc2(x,y))
        n+=1
        x=x+h
        X.append(x)
        Y.append(y)
    return X,Y

def D_2X(x):
    #return -9*sin(3*x)/(sqrt(x) + 5) + sin(3*x)/(2*x*(sqrt(x) + 5)**3) - 3*cos(3*x)/(sqrt(x)*(sqrt(x) + 5)**2) + sin(3*x)/(4*x**(3/2)*(sqrt(x) + 5)**2)
    return -20*x**3*sin(2*x)/(x**5 + x + 1)**2 + (-10*x**4 - 2)*(-5*x**4 - 1)*sin(2*x)/(x**5 + x + 1)**3 \
          + 4*(-5*x**4 - 1)*cos(2*x)/(x**5 + x + 1)**2 - 4*sin(2*x)/(x**5 + x + 1) # ТАНЯ

def plot_func (a,b):
    x = []
    y = []

    while (a <= b):
        y.append(inputFunc(a))
        x.append(a)
        a += 0.001

    plt.figure("График интеграла ")
    plt.grid(True)
    plt.fill_between(x,y, color='green', alpha=0.25,hatch = 'xx' )

    plt.plot(x, y)

def plot_D_2X(a,b,maxi):
    x=[]
    y=[]
    point = maxi
    while (a<=b):
        y.append(D_2X(a))
        x.append(a)
        a+=0.001

    plt.figure("График 2 производной ")
    plt.grid(True)

    plt.plot(x, y, point[0], -point[1], 'o')


def maxi(a,b):
    while (abs(b - a) >= 1e-9):
        x2 = a + (b - a)/PHI
        x1 = b - (b - a)/PHI

        if (abs(D_2X(x1))) <= (abs(D_2X(x2))):
            a = x1
        else :
            if (abs(D_2X(x1))) > (abs((D_2X(x2)))):
                b = x2

    return ((a + b) / 2), abs(D_2X((a+b)/2))

def trapez(a,b,M2):
    n = sqrt(M2 * (b - a) ** 3 / (12 * EPS))
    n = round(n)
    h=(b-a)/n
    I=0
    for i in range(1,n):
        j=a+h*i
        I+=inputFunc(j)
    res = ((inputFunc(a)+inputFunc(b))/2 + I)*h

    return res,n

def Simpson(a,b,m):
    n=2*m
    h = (b-a)/n
    I2=0
    I4=inputFunc(a+h)
    for i in range(1,m):
        I4+=inputFunc(a+(2*i+1)*h)
        I2+=inputFunc(a+2*i*h)

    I = (inputFunc(a) + inputFunc(b)+4*I4+2*I2)*h/3

    return I
if __name__ == "__main__":

    a = 0
    b = 1
    # n = 2m
    m = 5

    M2 = maxi(a,b)
    plot_D_2X(a,b,M2)
    plot_func(a,b)

    trap = trapez(a,b,M2[1])
    simp =Simpson(a,b,m)

    print("Трапеция | ","n = ",trap[1]," I ~ ",trap[0])

    print("Cимпсон  |  n = 10","   I ~ ",simp)
    print("-------------------------------------------------------")

    rk = Runge_Kuta3()
    print("Рунге-Кутта - 3 порядка:")
    print("y1 = ",np.round(rk[1][0],5))
    print("y2 = ",np.round(rk[1][1],5))
    print("y3 = ",np.round(rk[1][2],5))
    rk2 = Runge_Kuta4()
    print("Рунге-Кутта - 4 порядка:")
    print("y1 = ",np.round(rk2[1][0],5))
    print("y2 = ",np.round(rk2[1][1],5))
    print("y3 = ",np.round(rk2[1][2],5))
    pk = prognoz()
    print("Прогноза коррекции:")
    print("y1 = ",np.round(pk[1][0],5))
    print("y2 = ",np.round(pk[1][1],5))
    print("y3 = ",np.round(pk[1][2],5))
    plt.figure("Графики")
    plt.grid(True)
    leg1,leg2,leg3 = plt.plot(rk[0],rk[1],rk2[0], rk2[1],pk[0],pk[1])
    plt.legend((leg1,leg2,leg3),("Рунге- Кутта 3 порядка","Рунге- Кутта 4 порядка","Прогноза коррекции"))
    plt.grid(True)
    plt.show()
