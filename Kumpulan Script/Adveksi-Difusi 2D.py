#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import sys
#%%
def percentage(part, whole):
        percentage = 100 * float(part)/float(whole)
        return str(round(percentage,2)) + "%"
#%%
#MODUL 3 : ADVEKSI-DIFUSI 2D
#SKENARIO I
#%%
#input parameter awal
C = 1.30 #Arus/Konstanta Adveksi
ad = 1.30 #Konstanta Difusi
theta = 0 + 30 #Arah Arus (Geographics Convertion (0 degree N))
#parameter lanjutan
q = 0.95 #syarat kestabilan
x = 500 #Jumlah Grid x
y = 500 #Jumlah Grid y
dx = 5
dy = 5
Tend = 100 + 3 #Waktu Simulasi
dt = 0.5
#Polutan
px = 250 #Grid x polutan dibuang
py = 230 + 3 #Grid x polutan dibuang
Ka = 1000 + 30 #Jumlah Polutan (Ic)
#Perhitungan U & V (dipecah menjadi sumbu x dan y)
u = C * np.sin(theta*np.pi/180)
v = C * np.cos(theta*np.pi/180)
dt_count = 1/((abs(u)/(q*dx))+(abs(v)/(q*dy))+(2*ad/(q*dx**2))+(2*ad/(q*dx**2)))

#%%

Nx = int(x/dx) #jumlah mesh dalam arah x
Ny = int(x/dy) #jumlah mesh dalam arah x
Nt = int(Tend/dt) #Jumlah timestep
#perhitungan titik lepas polutan
px1 = int (px/dx)
py1 = int (py/dy)

#penyederhanaan fungsi
lx = u*dt/dx
ly = v*dt/dy
ax = ad*dt/dx**2
ay = ad*dt/dy**2
cfl = (2*ax + 2*ay + abs(lx) + abs(ly))
#perhitungan cfl
if cfl >= q:
    print('CFL Violated, Please use dt :'+str(round(dt_count,4)))
    sys.exit()
#%%
#Pembuatan Grid
x_grid = np.linspace(0-dx, x+dx, Nx+2) #Untuk GhostNode pada Ujung Boundary
y_grid = np.linspace(0-dy, y+dx, Ny+2) #Untuk GhostNode pada Ujung Boundary
t = np.linspace(0, Tend, Nt+1)
x_mesh,y_mesh = np.meshgrid(x_grid, y_grid)
F = np.zeros((Nt+1,Ny+2,Nx+2))
#Kondisi Awal
F[0,py1,px1] = Ka #Kondisi Awal
#%%
#Iterasi
for n in range (0, Nt):
    #F[n,py1,px1] = Ka #Gunakan code ini jika ingin melakukan kondisi kontinu
    for i in range (1,Ny+1):
        for j in range(1,Nx+1):
            F[n+1,i,j]=((F[n,i,j]*(1-abs(lx)-abs(ly)))  +                         (0.5*(F[n,i-1,j]*(ly+abs(ly)))) +                         (0.5*(F[n,i+1,j]*(abs(ly)-ly))) +                         (0.5*(F[n,i,j-1]*(lx+abs(lx)))) +                         (0.5*(F[n,i,j+1]*(abs(lx)-lx))) +                         (ay*(F[n,i+1,j]-2*(F[n,i,j])+F[n,i-1,j])) +                         (ax*(F[n,i,j+1]-2*(F[n,i,j])+F[n,i,j-1])))  #Diskritisasi
    #Kondisi Batas (Dirichlet Condition)
    F[n+1,0,:] = 0 #Batas Bawah
    F[n+1,:,0] = 0 #Batas Kiri
    F[n+1,Ny+1,:] = 0 #Batas Atas
    F[n+1,:,Nx+1] = 0 #Batas Kanan

#%%
    #Output Gambar
    plt.clf()
    plt.pcolor(x_mesh, y_mesh, F[n+1,:,:], cmap = 'seismic', shading = 'flat', edgecolor = 'k')
    cbar = plt.colorbar(orientation = 'vertical', shrink = 0.95, extend = 'both')
    cbar.set_label (label = 'Konsentrasi', size = 9)
    #plt.clim(0,100)
    plt.title('Ariestia Finola Damanik_26050120120003 - Skenario 1 \n t='+str(round(dt*(n+1),3))+', Kondisi Awal='+str(Ka),fontsize=10)
    plt.xlabel('x_grid',fontsize=9)
    plt.ylabel('y_grid',fontsize=9)
    plt.axis([0,x,0,y])
    #plt.pause(0.01)
    plt.savefig(str(n+1)+'.png', dpi = 300)
    plt.close()
    print('running timestep ke :' + str(n+1) + ' dari :' + str(Nt) + ' ('+ percentage(n+1,Nt)+')')


# In[ ]:




