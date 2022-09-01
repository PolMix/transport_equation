# -*- coding: utf-8 -*-
"""
Created on Sat Apr 10 22:47:19 2021

@author: Pol(Михайлов Павел 316 группа, задача 2 вариант 16)


"""

from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pylab 



"1. Построение характеристик для задачи"
x1, x2 = -2, 1
t1, t2 = 0, 2
NC=2


plt.ylim(t1,t2)
plt.xlim(x1,x2)
plt.ylabel('t', size = 15)
plt.xlabel('x', size = 15)
plt.grid(True)
plt.title('Семейство характеристик')
for C in range(-NC-1,NC+2,1):
    plt.plot([x1,x2], [-x1 + C/NC, -x2 + C/NC], label='C='+str(round(C/NC,2)))
plt.axvline(x=-1.0, linewidth = 2.5, linestyle = '--')
plt.axvline(x=0.0, linewidth = 2.5, linestyle = '--')
plt.legend()
mpl.rcParams['figure.dpi'] = 150
plt.show()



"2. Зададим число шагов по временной и пространственной координатам"
x1, x2 = -1, 0
t1, t2 = 0, 2
Nx = 100
Nt = 400
e = 0.000001 #Погрешность для метода Ньютона

h = -(x2 - x1)/Nx
tau = (t2 - t1)/Nt
x = np.linspace(x2, x1, Nx)
t = np.linspace(t1, t2, Nt)



"3. Зададим массив для нашей сеточной функции, определим граничные/начальные условия"
y = np.zeros((Nx, Nt))
for n in range(0,Nx): #Начальное условие
    y[n][0] = np.log(2/(1 - 2*x[n]))
    
for k in range(0,Nt): #Граничное условие
    y[0][k] = np.log(2/(np.exp(t[k]) - 2*t[k]*np.exp(-t[k])))

    
    
"4. Задание функций f и f'  "
def f(n,k):
    return 1/tau * (y[n+1][k+1] - y[n+1][k]) - 1/h * (y[n+1][k+1] - y[n][k+1]) + np.exp(t[k+1] + y[n+1][k+1]) - 1

def fd(n,k):
    return 1/tau - 1/h + np.exp(t[k+1] + y[n+1][k+1])



"5. Реализация метода Ньютона"
for n in range(0,Nx-1):
    for k in range(0,Nt-1):
        y[n+1][k+1] = 0.5*(y[n+1][k] + y[n][k+1]) #Приближение для следующего шага
        while(np.abs(f(n,k) / fd(n,k)) > e):
            y[n+1][k+1] = y[n+1][k+1] - f(n,k) / fd(n,k)


  
"6. Аналитическое решение"
U = np.zeros((Nx,Nt))
for n in range(0,Nx):
    for k in range(0,Nt):
        U[n][k] = np.log(2/(np.exp(t[k]) - 2*(x[n]+t[k])*np.exp(-t[k])))

        

"7. Построение аналитического решения"
T,X = np.meshgrid(t,x)
fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')

surf = axes.plot_surface(X, T, U[:,:], rstride=4, cstride=4, cmap = cm.jet)
pylab.xlabel('$x$', size=15, labelpad = 5)
plt.xticks(rotation = 0)
pylab.ylabel('$t$', size=15, labelpad = 5)
plt.yticks(rotation = 0)
 
axes.set_zlabel('U(x,t)', size = 15, labelpad = 10)

axes.view_init(20, 160) #поворот графика
fig.colorbar(surf, ax = [axes], pad = 0.05, location='left')
plt.title('График аналитического решения U(x, t) \n Nx=' +str(Nx) + '  Nt=' 
          +str(Nt)+ '  tau/h=' +str(np.abs(tau/h))+ '  e=' +str(e))
plt.show()
    

        
"8. Построение численного решения"
T,X = np.meshgrid(t,x)
fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')

surf = axes.plot_surface(X, T, y[:][:], rstride=4, cstride=4, cmap = cm.jet)
pylab.xlabel('$x$', size=15, labelpad = 5)
plt.xticks(rotation = 0)
pylab.ylabel('$t$', size=15, labelpad = 5)
plt.yticks(rotation = 0)

axes.set_zlabel('y(x,t)', size = 15, labelpad = 10)


axes.view_init(20, 155) #поворот графика
fig.colorbar(surf, ax = [axes], pad = 0.05, location='left')
plt.title('График численного решения y(x, t) \n Nx=' +str(Nx) + '  Nt=' 
          +str(Nt)+ '  tau/h=' +str(np.abs(tau/h))+ '  e=' +str(e))
plt.show()



"9. Построение погрешности решения"
T,X = np.meshgrid(t,x)
fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')

surf = axes.plot_surface(X, T, np.abs(U[:][:] - y[:][:]), rstride=4, cstride=4, cmap = cm.jet)
pylab.xlabel('$x$', size=15, labelpad = 5)
plt.xticks(rotation = 0)
pylab.ylabel('$t$', size=15, labelpad = 5)
plt.yticks(rotation = 0)

axes.set_zlabel('y_Error(x,t)', size = 15, labelpad = 10)



axes.view_init(20, 120) #поворот графика
plt.colorbar(surf, ax = [axes], pad = 0.05, location='left')
plt.title('График погрешности |U(x,t) - y(x, t)| \n Nx=' +str(Nx) + '  Nt=' 
          +str(Nt)+ '  tau/h=' +str(np.abs(tau/h))+ '  e=' +str(e))
plt.show()
print(np.abs(U[Nx-1][Nt-1] - y[Nx-1][Nt-1]))
print(np.sqrt(tau*tau+h*h))



"Построение зависимости прогрешности в узле от величины шагов сетки"
Nx_1 = 40
Nx_2 = 2500
Nx_N = 14

Nt_1 = 40
Nt_2 = 2500
Nt_N = 14

Nx = np.linspace(Nx_1,Nx_2,Nx_N)
Nt = np.linspace(Nt_1,Nt_2,Nt_N)
h = np.zeros((Nx_N))
tau = np.zeros((Nt_N))

for a in range(0,Nx_N):
    h[a] = (x2 - x1)/Nx[a]
for b in range(0,Nt_N):
    tau[b] = (t2 - t1)/Nt[b]

Ry = np.zeros((Nx_N,Nt_N))
rU = np.zeros((Nx_N,Nt_N))

for a in range(0,Nx_N):
    for b in range(0,Nt_N):
        Ex = np.linspace(x2,x1,int(Nx[a]))
        Et = np.linspace(t1,t2,int(Nt[b]))
        Ey = np.zeros((int(Nx[a]), int(Nt[b])))
        for n in range(0,int(Nx[a])): #Начальное условие
            Ey[n][0] = np.log(2/(1 - 2*Ex[n]))
            
        for k in range(0,int(Nt[b])): #Граничное условие
            Ey[0][k] = np.log(2/(np.exp(Et[k]) - 2*Et[k]*np.exp(-Et[k])))
            
            
        "Задание функций f и f'  "
        def Ef(n,k,a,b):
            return 1/tau[b]*(Ey[n+1][k+1] - Ey[n+1][k]) + 1/h[a]*(Ey[n+1][k+1] - Ey[n][k+1]) + np.exp(Et[k+1] + Ey[n+1][k+1]) - 1
        
        def Efd(n,k,a,b):
            return 1/tau[b] + 1/h[a] + np.exp(Et[k+1] + Ey[n+1][k+1])
        
        
        "Реализация метода Ньютона"
        for n in range(0,int(Nx[a])-1):
            for k in range(0,int(Nt[b])-1):
                Ey[n+1][k+1] = 0.5*(Ey[n+1][k]+Ey[n][k+1]) #Приближение для следующего шага
                while(np.abs(Ef(n,k,a,b) / Efd(n,k,a,b)) > e):
                    Ey[n+1][k+1] = Ey[n+1][k+1] - Ef(n,k,a,b) / Efd(n,k,a,b)
                    
                    
        "Аналитическое решение"
        eU = np.zeros((int(Nx[a]),int(Nt[b])))
        for n in range(0,int(Nx[a])):
            for k in range(0,int(Nt[b])):
                eU[n][k] = np.log(2/(np.exp(Et[k]) - 2*(Ex[n]+Et[k])*np.exp(-Et[k])))
        
        "Выбираем узел, близкий к самой далекой точке от ГУ и НУ"
        Ry[a][b] = Ey[int(Nx[a]) - 1][int(Nt[b]) - 1]
        rU[a][b] = eU[int(Nx[a]) - 1][int(Nt[b]) - 1]


"9. Построение ошибки узла"
NT,NX = np.meshgrid(Nt,Nx)
fig = plt.figure()
axes = fig.add_subplot(111, projection='3d')
surf = axes.plot_surface(NX, NT, np.abs(rU[:][:] - Ry[:][:]), rstride=1, cstride=1, cmap = cm.jet)
pylab.xlabel('$h$', size=15, labelpad = 5)
plt.xticks(rotation = 0)
pylab.ylabel('$tau$', size=15, labelpad = 5)
plt.yticks(rotation = 0)

axes.set_zlabel('|RU - Ry|', size = 15, labelpad = 10)
axes.view_init(20, 130) #поворот графика
plt.colorbar(surf, ax = [axes], pad = 0.05, location='left')
plt.title('График погрешности |RU(x,t) - Ry(x, t)| в узле \n (' +str(Nx_1) + ',' 
          +str(Nx_2)+','+str(Nx_N)+') ; (' +str(Nt_1)+ ',' +str(Nt_2)+ ','
          +str(Nt_N)+');  e=' +str(e))
plt.show()


TAU,H = np.meshgrid(tau,h)
plt.contourf(H, TAU, np.abs(rU[:][:] - Ry[:][:]), cmap = cm.jet)
plt.xlabel('$h$', size=15, labelpad = 5)
plt.xticks(rotation = 0)
plt.ylabel('$tau$', size=15, labelpad = 5)
plt.yticks(rotation = 0)

axes.set_zlabel('|RU - Ry|', size = 15, labelpad = 10)
plt.colorbar(pad = 0.12)
plt.title('График погрешности |RU(x,t) - Ry(x, t)| \n (' +str(Nx_1) + ',' 
          +str(Nx_2)+','+str(Nx_N)+') ; (' +str(Nt_1)+ ',' +str(Nt_2)+ ','
          +str(Nt_N)+');  e=' +str(e))
plt.show()