# -*- coding: utf-8 -*-
"""
Created on Wed Mar 31 09:24:27 2021

@author: Pol (Михайлов Павел, 316 группа, задача 1, вариант 32)
"""

import numpy as np
import pylab 
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm



"1. Зададим пределы области исследования по X, Y, T"
x1, x2 = 0, np.pi
y1, y2 = 0, 3
t1, t2 = 0, np.pi



"2. Зададим количество шагов сетки N и величину шага сетки h"
Nx = 100
Ny = 5
Nt = 200
hx = (x2-x1)/Nx
hy = (y2-y1)/Ny
tau = (t2-t1)/Nt



"3. Создадим массивы переменных"
x = np.linspace(x1,x2,Nx)
y = np.linspace(y1,y2,Ny)
t = np.linspace(t1,t2,Nt)
U_An = np.zeros((Nx,Ny,Nt*2+1))
U_Num = np.zeros((Nx,Ny,Nt*2+1))



"АНАЛИТИЧЕСКОЕ РЕШЕНИЕ"
"4. Вычислим значения аналитического решения в узлах сетки"
for n in range(0,Nx):
    for m in range(0,Ny):
        for k in range(0,Nt*2+1):
            U_An[n,m,k] = 1/17*np.cos(x[n])*(4*np.sin(tau*k/2)-np.cos(tau*k/2)+18*np.exp(-4*tau*k/2))





"ЧИСЛЕННОЕ РЕШЕНИЕ"
"5. Зададим необходимые коэффициенты для численного решения"
Ax = 2 * tau/hx/hx
Cx = 1 + 2 * Ax
Ay = 2 * tau/hy/hy
Cy = 1 + 2 * Ay



"6. Зададим функцию F из разностной схемы для х и у"
def Fx(a,b,c):
    return Ay*(U_Num[a,b-1,c-1] + U_Num[a,b+1,c-1]) + (1-2*Ay)*U_Num[a,b,c-1] + 0.5*tau*np.cos(x[a])*np.sin(tau*(c + 1)/2)


def Fy(a,b,c):
    return Ax*(U_Num[a-1,b,c-1] + U_Num[a+1,b,c-1]) + (1-2*Ax)*U_Num[a,b,c-1] + 0.5*tau*np.cos(x[a])*np.sin(tau*(c - 1)/2)



"7. Зададим функции, реализующие метод прогонки по x и y"
"Прогонка по x"
def SweepX(b,c):
    alpha = np.zeros(Nx)
    beta = np.zeros(Nx)
    alpha[1] = 1 #У нас ГУ Неймана при х = 0
    beta[1] = 0
    A = Ax
    B = Ax
    C = 1+2*Ax
    
    for n in range(1,Nx-1): #Прямая прогонка
        alpha[n+1] = B / (C - A*alpha[n])
        beta[n+1] = (A*beta[n] + Fx(n,b,c)) / (C - A*alpha[n])
    U_Num[Nx-1,b,c] = beta[Nx-1] / (1 - alpha[Nx-1])  # Из ГУ Неймана при х = \pi
    for n in range(Nx-1,0,-1): #Обратная прогонка
        U_Num[n-1,b,c] = alpha[n]*U_Num[n,b,c] + beta[n]
    for n in range(0,Nx): #Установка граничных значений
        U_Num[n,0,c] = U_Num[n,1,c]
        U_Num[n,Ny-1,c] = U_Num[n,Ny-2,c]


"Прогонка по y"
def SweepY(a,c):
    alpha = np.zeros(Ny)
    beta = np.zeros(Ny)
    alpha[1] = 1  #У нас ГУ Неймана при у = 0
    beta[1] = 0
    A = Ay
    B = Ay
    C = 1+2*Ay
    
    for m in range(1,Ny-1): #Прямая прогонка
        alpha[m+1] = B / (C - A*alpha[m])
        beta[m+1] = (A*beta[m] + Fy(a,m,c)) / (C - A*alpha[m])
    U_Num[a,Ny-1,c] = beta[Ny-1] / (1 - alpha[Ny-1]) # Из ГУ Неймана при у = 3
    for m in range(Ny-1,0,-1): #Обратная прогонка
        U_Num[a,m-1,c] = alpha[m]*U_Num[a,m,c] + beta[m]
    for m in range(0, Ny): #Установка граничных значений
        U_Num[0,m,c] = U_Num[1,m,c]
        U_Num[Nx-1,m,c] = U_Num[Nx-2,m,c]



"7. Применим метод прогонки"
for n in range(0,Nx): #Задание ГУ
    for m in range(0,Ny):
        U_Num[n, m, 0] = np.cos(x[n])

for k in range(1,2*Nt,2):
    for m in range(1,Ny-1):
        SweepX(m,k)
    for n in range(1,Nx-1):
        SweepY(n,k+1)
    
U_AnNum = np.zeros(((Nx,Ny,Nt*2+1)))
for n in range(0,Nx):
    for m in range(0,Ny):
        for k in range(0,Nt*2+1):
            U_AnNum[n,m,k] = np.abs(U_An[n][m][k] - U_Num[n][m][k])
        

        
"ГРАФИКИ"
dpi0 = 150
sk0 = int(Nt/5)
for k0 in range(0,Nt+1,sk0): 
    "8. Построение аналитического решения"
    fig = pylab.figure(dpi=dpi0)
    axes = Axes3D(fig)
    
    Y,X = np.meshgrid(y,x)
    
    surf = axes.plot_surface(X,Y, U_An[:,:,k0], rstride=4, cstride=4, cmap = cm.jet)
    axes.view_init(20, 45) #поворот графика
    
    pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
    pylab.ylabel('$y$', size=15, labelpad = 5)
    axes.set_zlabel('U_Analytical(x,y)', size = 15, labelpad = 10)
    
    pylab.title('График U_An (x, y, t) в момент t = ' +str(round(k0*tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
    cbar = fig.colorbar(surf, ax=axes) #boundaries = [-1.00, -0.75, -0.50, -0.25, 0, 0.25, 0.50, 0.75, 1.00]
    pylab.show()
    
    
    
   
    "9. Построение численного решения"
    fig_num = pylab.figure(dpi=dpi0)
    axes_num = Axes3D(fig_num)
    
    surf_num = axes_num.plot_surface(X,Y, U_Num[:,:,k0], rstride=4, cstride=4, cmap = cm.jet)
    axes_num.view_init(20, 45) #поворот графика
    
    pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
    pylab.ylabel('$y$', size=15, labelpad = 5)
    axes_num.set_zlabel('U_Numerical(x,y)', size = 15, labelpad = 10)
    
    pylab.title('График U_Num (x, y, t) в момент t = ' +str(round(k0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
    cbar = fig_num.colorbar(surf_num, ax=axes_num)
    pylab.show()
    
    
    
    
    "10. Построение разности численного и аналитического решений, погрешности"
    fig_annum = pylab.figure(dpi=dpi0)
    axes_annum = Axes3D(fig_annum)
    
    surf_annum = axes_annum.plot_surface(X,Y, U_An[:,:,k0] - U_Num[:,:,k0], rstride=4, cstride=4, cmap = cm.jet)
    axes_annum.view_init(20, 45) #поворот графика
    
    pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
    pylab.ylabel('$y$', size=15, labelpad = 5)
    axes_annum.set_zlabel('U_Error(x,y)', size = 15, labelpad = 10)
    
    pylab.title('График (U_An - U_Num) (x, y, t) в момент t = ' +str(round(k0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
    cbar = fig_annum.colorbar(surf_annum, ax=axes_annum)
    pylab.show()
    
    
    
    
    "11. Построение модуля разности численного и аналитического решений, погрешности"
    fig_annum = pylab.figure(dpi=dpi0)
    axes_annum = Axes3D(fig_annum)
    
    surf_annum = axes_annum.plot_surface(X,Y, np.abs(U_An[:,:,k0] - U_Num[:,:,k0]), rstride=4, cstride=4, cmap = cm.jet)
    axes_annum.view_init(20, 45) #поворот графика
    
    pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
    pylab.ylabel('$y$', size=15, labelpad = 5)
    axes_annum.set_zlabel('|U_Error(x,y)|', size = 15, labelpad = 10)
    
    pylab.title('График |U_An - U_Num| (x, y, t) в момент t = ' +str(round(k0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
    cbar = fig_annum.colorbar(surf_annum, ax=axes_annum)
    pylab.show()
    
    
    
"12. График аналитического решения от X и T"
spt = np.linspace(t1,t2,2*Nt+1)
m0 = 0
fig_annum = pylab.figure(dpi=dpi0)
axes_annum = Axes3D(fig_annum)

SPT, X = np.meshgrid(spt, x)

surf_annum = axes_annum.plot_surface(X, SPT, U_An[:,m0,:], rstride=4, cstride=4, cmap = cm.jet)
axes_annum.view_init(10, 40) #поворот графика

pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
pylab.ylabel('$t$', size=15, labelpad = 5)
axes_annum.set_zlabel('U_An(x,t)', size = 15, labelpad = 10)

pylab.title('График U_An (x, y, t) в момент y = ' +str(round(m0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
cbar = fig_annum.colorbar(surf_annum, ax=axes_annum)
pylab.show()



"13. График численного решения от X и T"
fig_annum = pylab.figure(dpi=dpi0)
axes_annum = Axes3D(fig_annum)

surf_annum = axes_annum.plot_surface(X, SPT, U_Num[:,m0,:], rstride=4, cstride=4, cmap = cm.jet)
axes_annum.view_init(10, 40) #поворот графика

pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
pylab.ylabel('$t$', size=15, labelpad = 5)
axes_annum.set_zlabel('U_Num(x,t)', size = 15, labelpad = 10)

pylab.title('График U_Num (x, y, t) в момент y = ' +str(round(m0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
cbar = fig_annum.colorbar(surf_annum, ax=axes_annum)
pylab.show()



"14. График модуля ошибки решения от X и T"
fig_annum = pylab.figure(dpi=dpi0)
axes_annum = Axes3D(fig_annum)

surf_annum = axes_annum.plot_surface(X, SPT, np.abs(U_An[:,m0,:] - U_Num[:,m0,:]), rstride=4, cstride=4, cmap = cm.jet)
axes_annum.view_init(10, 40) #поворот графика

pylab.xlabel('$x$', size=15, labelpad = 5) #labelpad: отступы от оси для подписей
pylab.ylabel('$t$', size=15, labelpad = 5)
axes_annum.set_zlabel('|U_Error(x,t)|', size = 15, labelpad = 10)

pylab.title('График |U_An - U_Num| (x, y, t) в момент y = ' +str(round(m0 * tau,2)) + ' \n Nt=' +str(Nt)+'  Nx='+str(Nx)+'  Ny='+str(Ny))
cbar = fig_annum.colorbar(surf_annum, ax=axes_annum)
pylab.show()