# Tugas Akhir Prk.Pemodelan Oseanografi (Team 11)🗂
Repositori ini dibuat untuk memenuhi Tugas Akhir Praktikum Pemodelan Oseanografi 2022. _Repository_ ini memuat file berupa _script phyton (py)_ yang dapat memproses beberapa pemodelan oseanografi seperti Adveksi-Difusi dan Hidrodinamika. Pengerjaan untuk _repository_ kali ini, menggunakan bahasa pemrograman _python_ yang dapat dilakukan pada beberapa platform seperti _Google Colaboratory_ dan _Jupyter Notebook_. Sedangkan untuk _library_ yang digunakan kali ini adalah _Numpy, Matplotlib,_ dan _NDBC_ (dari siphon.simplewebservice.ndbc). Seluruh _script_ yang dibuat adalah hasil Team 11 Oseanografi 2020. Semoga dapat bermanfaat!

## 1. AUTHORS (TEAM 11)👩🏻‍🧑🏻‍💻 
1. Ariestia Finola Damanik 26050120120003 A
2. Barriel Jimly Al Ariefanzah 26050120140122 B
3. Dimas Sukma Hadi 26050120140170 B
4. Fahri Rahmalia 26050120130037 A
5. Nisa Uswatun Khasanah 26050120130119 A
6. Pilipus Najuita Nduru 26050120120008 A
7. Yvetty Zilla Durhan 26050120120021 A

## 2. Cara Penggunaan _Script_ 💬
1. Pengguna dapat membuka folder **Kumpulan Script** pada repository ini
2. Di dalamnya terdapat 3 _script_, bisa dibuka salah satu
![image](https://user-images.githubusercontent.com/89583653/169198940-5e831637-1628-4be3-b0be-a87c5226ee7f.png)
3. Nanti akan muncul _script_-nya seperti gambar di bawah ini
![image](https://user-images.githubusercontent.com/89583653/169198974-f7dc70ad-7edf-4b67-85f1-3685a3914558.png)
4. _Script_ dapat di-_copy_ untuk nantinya dipergunakan dan disesuaikan di _Jupyter Notebook_

## 3. Metode Pengerjaan🛠️
1. Adveksi-Difusi 2D
2. Hidrodinamika 1D
3. Hidrodinamika 2D

### 3.1. Adveksi-Difusi 2D
```
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
```

Script Perhitungan
```
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
```

Hasil / Output Gambar
```
#Output Gambar
    plt.clf()
    plt.pcolor(x_mesh, y_mesh, F[n+1,:,:], cmap = 'seismic', shading = 'flat', edgecolor = 'k')
    cbar = plt.colorbar(orientation = 'vertical', shrink = 0.95, extend = 'both')
    cbar.set_label (label = 'Konsentrasi', size = 9)
    #plt.clim(0,100)
    plt.title('Team 11 - Skenario 1 \n t='+str(round(dt*(n+1),3))+', Kondisi Awal='+str(Ka),fontsize=10)
    plt.xlabel('x_grid',fontsize=9)
    plt.ylabel('y_grid',fontsize=9)
    plt.axis([0,x,0,y])
    #plt.pause(0.01)
    plt.savefig(str(n+1)+'.png', dpi = 300)
    plt.close()
    print('running timestep ke :' + str(n+1) + ' dari :' + str(Nt) + ' ('+ percentage(n+1,Nt)+')')

```
Contoh dari hasil akan berada setelah bagian Tutorial pengerjaan Script dan Skenario Mode.

📝Dasar Teori📝
- Adveksi 2D
  Adveksi 2D merupakan proses pergerakan substansi oleh fluida yang dipengaruhi oleh gaya-gaya tertentu dalam arah vertikal dan horizontal.
  Berikut adalah persamaannya:
 
  ![image](https://user-images.githubusercontent.com/90328732/169702448-0982b791-22e9-4992-afac-e32613a8516a.png)
- Difusi 2D
  Difusi 2D merupakan proses penyebaran substansi dari substansi yang lebih tinggi ke yang lebih rendah dipengaruhi oleh gaya-gaya tertentu dalam arah  vertikal dan horizontal.
  Berikut adalah persamaannya:
  
  ![image](https://user-images.githubusercontent.com/90328732/169702605-8ee91b16-5398-40f9-a741-01ce1eb7e372.png)

📝Persamaan Pembangun Model Adveksi-Difusi 2D📝
  
  Model adveksi-difusi 2D adalah model matematika yang menggambarkan proses transportasi suatu zat yang dipengaruhi gaya gravitasi dan penyebaran sekaligus.
  Persamaan adveksi-difusi 2D adalah model matematika yang menggambarkan proses transportasi suatu zat yang dipengaruhi gaya dalam dua dimensi. 
  Berikut adalah persamaannya:
  
  ![image](https://user-images.githubusercontent.com/90328732/169703024-b449b4e6-4413-426e-8980-a57c0c301fcf.png)

📝Deskritisasi Persamaan Adveksi-Difusi 2D📝
 
 Untuk membentuk suatu persamaan model 2D yang mendekati proses kejadian di alam maka perlu adanya deskritisasi terhadap persamaan tersebut. Deskritisasi merupakan suatu metode untuk mencari solusi persamaan secara numerik dari suatu persamaan matematika sehingga dapat dinyatakan baik dalam dimensi ruang ataupun waktu. Proses deksritisasi model 2D pada bagian atau suku adveksi umumnya menggunakan metode eksplisit upstream. Metode yang sama juga berlaku untuk deskritisasi suku difusi.
  Metode eksplisit upstream (pada model 2D adveksi) merupakan metode eksplisit dimana persamaan beda hingga dengan metode ini menggunakan pendekatan beda maju untuk turunan waktu, sedangkan untuk turunan terhadap ruang dilakukan dengan melihat arah kecepatan u. Jika u > 0 maka turunan terhadap ruang menggunakan pendekatan beda mundur, sebaliknya jika u < 0 digunakan pendekatan beda maju. Persamaan dari metode diskritisasi untuk suku adveksi 2D adalah sebagai berikut :
  
  ![image](https://user-images.githubusercontent.com/90328732/169703189-1a4ade5c-e16f-4930-b4f3-96e122854d5e.png)

  Model 2D untuk mekanisme transpor difusi dapat menggunakan pendekatanbeda maju untuk turunan waktu dan beda pusat untuk turunan ruang. Indeks n untuk waktu, indeks i untuk ruang, dan koefisiesn difusi AD dianggap konstan terhadap ruang dan waktu. Persamaan diskritisasi untuk model 2D difusi adalah sebagai berikut:
  
  ![image](https://user-images.githubusercontent.com/90328732/169703312-8e9401fa-f320-455c-a0c9-e33ee09a2bee.png)

  Pada model 2D untuk proses adveksi dan difusi yang telah digabung maka deskritisasi persamaan adalah menggabungkan dua suku yakni suku adveksi dan suku difusi. Persamaan diskritisasi untuk model adveksi difusi 2D adalah sebagai berikut :
  
  ![image](https://user-images.githubusercontent.com/90328732/169703349-94281e2d-5645-4db8-9075-f886819c4527.png)

📝Penentuan Syarat dan Nilai Batas📝

Syarat batas merupakan suatu kondisi yang menggambarkan kondisi di batas (ruang maupun waktu) dari model yang dibangun. Pada model 2D, syarat batas dari metode eksplisit upstream diberikan pada nilai awal (hulu) dan nilai akhir (hilir). Syarat batas di hulu dan di hilir adalah sebagai berikut :

![image](https://user-images.githubusercontent.com/90328732/169703476-c748878d-2f45-4bdf-8e93-a9a88a1ce399.png)

Iterasi kemudian akan dihentikan apabila telah mencapai batas:

![image](https://user-images.githubusercontent.com/90328732/169703501-5e61afa9-5111-46a3-a7a0-171b97d94e47.png)

 Syarat awal yang digunakan dalam skenario model 2D adveksi-difusi ini adalahdengan memberikan harga 0 disemua titik konsentrasi polutan kecuali di titik-titiksumber yang tersebar dan sumber bersifat tidak kontinu.

📝Kriteria Kestabilan📝

Kriteria kestabilan merupakan suatu metode untuk menentukan seberapa besarnilai stabilitas dari model yang dibangun. Kriteria kestabilan yang digunakan untuk menyelesaikan pemodelan 2D adveksi difusi ini adalah sebagai berikut:

![image](https://user-images.githubusercontent.com/90328732/169703585-428001bb-a5cb-48f5-942b-a08be0cfb25a.png)


📄Tutorial Pengerjaan Script dan Skenario Model📄
1. _Mandatory library python matploblib_ dimasukkan supaya dapat diberikan efek visual pada grafik, kemudian _numpy_ untuk numerik, dan _sys_ untuk mengakses konfigurasi interpreter pada saat _runtime_. Pendefinisian juga dimasukkan.
2. Parameter perhitungan polutan awal dimasukkan.
3. _Script_ perhitungan dibuat untuk mengetahui persebaran polutan.
4. _Script_ untuk _grid_ dibuat sebagai alat bantu untuk melihat persebaran polutan.
5. Iterasi dilakukan hingga seluruh syarat batas terpenuhi.
6. _Script output_ gambar dibuat untuk mendapatkan gambar penyebaran polutan.
7. _Script_ dapat di-_run_.
-

📄Hasil / Output gambar📄

![1](https://user-images.githubusercontent.com/105839721/169779931-9d28d4f4-dce0-460c-b8fe-e453967d42e0.jpg)

![101](https://user-images.githubusercontent.com/105839721/169779967-c2fefb79-0f70-42c3-960a-ac26a0b122b1.jpg)

![200](https://user-images.githubusercontent.com/105839721/169779983-e2037773-9479-40e4-b324-25b400f48f94.jpg)


### 3.2. Hidrodinamika 1D
```
#!/usr/bin/env python
# coding: utf-8

# In[1]:


import matplotlib.pyplot as plt
import numpy as np

#
# Proses Awal
# 

p = 5000 #Panjang Grid
T = 1200 #Waktu Simulasi
A = 0.5 #Amplitudo
D = 15 #Depth/Kedalaman
dt = 2
dx = 100
To = 300 #Periode

g = 9.8
pi = np.pi
C = np.sqrt(g*D) #Kecepatan Arus
s = 2*pi/To #Kecepatan Sudut Gelombang
L = C*To #Panjang Gelombang
k = 2*pi/L #Koefisien Panjang Gelombang
Mmax = int(p//dx)
Nmax = int(T//dt)
```

Script Perhitungan
```
zo = [None for _ in range(Mmax)]
uo = [None for _ in range (Mmax)]

hasilu = [None for _ in range(Nmax)]
hasilz = [None for _ in range (Nmax)]

for i in range(1, Mmax+1) :
    zo[i-1] = A*np.cos(k*(i)*dx)
    uo[i-1] = A*C*np.cos(k*((i)*dx+(0.5)*dx))/(D+zo[i-1])
for i in range(1, Nmax+1) :
    zb = [None for _ in range(Mmax)]
    ub = [None for _ in range (Mmax)]
    zb[0] = A*np.cos(s*(i)*dt)
    ub[-1] = A*C*np.cos(k*L-s*(i)*dt)/(D+zo[-1])
    for j in range(1, Mmax) :
        ub[j-1] = uo[j-1]-g*(dt/dx)*(zo[j]-zo[j-1])
    for k in range (2, Mmax+1) :
        zb[k-1] = zo[k-1]-(D+zo[k-1])*(dt/dx)*(ub[k-1]-ub[k-2])
        hasilu[i-1] = ub
        hasilz[i-1] = zb
    for p in range (0, Mmax) :
        uo[p] = ub[p]
        zo[p] = zb[p]
```

Hasil Gambar
```
#
# PEMBUATAN GRAFIK
#
def rand_col_hex_string():
    return f'#{format(np.random.randint(0,16777215), "#08x")[2:]}'
  
hasilu_np = np.array(hasilu)
hasilz_np = np.array(hasilz)

fig0, ax0 = plt.subplots(figsize=(12,8))
for i in range (1, 16) :
    col0 = rand_col_hex_string()
    line, = ax0.plot(hasilu_np[:,i-1], c=col0, label=f'n={i}')
    ax0.legend()

    ax0.set(xlabel='Waktu', ylabel='Kecepatan Arus',
            title=''' Team 11_Praktikum Pemodelan Oseanografi 2022
            Perubahan Kecepatan Arus Dalam Grid Tertentu di Sepanjang Waktu''')
    ax0.grid()
```
![1](https://user-images.githubusercontent.com/105837184/169705728-28e87ab8-1f48-4008-8ce1-9ddf725e809d.png)

-
-

### 3.3. Hidrodinamika 2D
📝Dasar Teori📝

Model Hidrodinamika adalah simulasi suatu aliran yang didasarkan pada persamaan matematika dengan menggambarkan fenomena fisik aliran dan penyelesaian persamaan matematika secara numerik. Model Hidrodinamika 2D sendiri merupakan pemodelan hidrodinamika yang menggunakan lebih dari 1 parameter. Sesuai dengan namanya, 2D, tentu memiliki 2 bentuk, entah panjang maupun lebar. Dalam model hidrodinamika 2D sendiri dapat diartikan dengan penggambaran maupun tiruan dari suatu fenomena, kejadian maupun proses pada kurun waktu tertentu dengan memperhatikan suatu parameter dan anomali yang terjadi.

📝Parameter dan Anomali Hidrodinamika 2D📝
- Parameter

Parameter yang digunakan pada pemodelan hidrodinamika 2D meliputi lebih dari satu parameter. Contohnya, untuk memodelkan gelombang digunakan parameter arah angin dan juga kecepatan angin ataupun tekanan.
- Anomali

Anomali atau penyimpangan yang terjadi saat memodelkan hidrodinamika 2D biasanya diakibatkan oleh kondisi lapangan sesungguhnya dari parameter yang digunakan. Contohnya, dalam melakukan pemodelan gelombang dengan parameter angin. Tiupan angin yang berbeda-beda di berbagai situasi dapat mengakibatkan pemodelan tidak 100% sesuai dengan kondisi aslinya. 

📝Karakteristik Pemodelan Hidrodinamika 2D📝
- Medan dipresentasikan sebagai hasil permukaan yang kontinu (x,y)
- Kedalaman air tidak dianggap seragam
- Kecepatan tidak dianggap seragam
- Baik digunakan untuk gradien yang curam

Dalam pelaksanaan Praktikum pemodelan Oseanografi yang telah dilakukan, dipelajari penggambaran serta interpretasi model di wilayah perairan, baik berupa aliran panjang dan dengan kedalaman tertentu. Dengan penggunaan model hidrodinamika 2 dimensi dapat diketahui identifikasi ketebalan persebaran _oil spill_ di suatu perairan. Tidak hanya itu, model hidrodinamika 2d juga dapat digunakan untuk memodelkan dan menganalisis suatu fenomena seperti pemodelan gelombang akibat gaya pembangkit angin, pemodelan sampah plastik di laut dan pemodelan _coastal dynamics_ dan sedimentasi pantai.  Untuk lebih jelasnya, dapat dilihat contoh model data gelombang _National Buoy Data Center (NDBC)_ di Peraian Hilo, Hawaii dekat Samudera Pasifik Utara dalam _script_ berikut:

**Contoh Pemodelan Hidrodinamika 2D**
1. _Script_ dapat diambil melalui folder yang ada di _repository_ ini.
2. Pada awal _script_ diberikan keterangan terlebih dahulu, dimana data-data diambil.
```
# In[1]:


# Copyright (c) 2018 Siphon Contributors.
# Distributed under the terms of the BSD 3-Clause License
#SPDX-License-Identifier: BSD-3-Clause
"""
NBDC Buoy Meteorological Data Request
=====================================
The NDBC keeps a 45-day recent rolling file for each buoy. This examples shows how to access
the basic meteorological data from a buoy and make a simple plot.
"""
```
3. Dilakukan penginputan data dan _library_.
```
import matplotlib.pyplot as plt

from siphon.simplewebservice.ndbc import NDBC

###################################################
# Get a pandas data framw of all of the observation, meteorological data is the default
# observation set the query.
df = NDBC.realtime_observations('51003') #Station ID
df.head()
```
4. Membuat _plot time series_ untuk _pressure, wind speed, gust, direction,_ dan _water temperature_.
```
###################################################
# Let's make a simple time series plot to checkout what the data look like
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10))
ax2b = ax2.twinx()

# Pressure
ax1.plot(df['time'], df['pressure'], color='black')
ax1.set_ylabel('pressure [hPa]')
fig.suptitle('Team 11_Praktikum Pemodelan Oseanografi 2022', fontsize=18)


# Wind speed, gust, direction
ax2.plot(df['time'], df['wind_speed'], color='tab:orange')
ax2.plot(df['time'], df['wind_gust'], color='tab:olive', linestyle='--')
ax2b.plot(df['time'], df['wind_direction'], color='tab:blue', linestyle='-')
ax2.set_ylabel('Wind Speed [m/s]')
ax2.set_ylabel('Wind Direction')


#Water temperature
ax3.plot(df['time'], df['water_temperature'], color='tab:brown')
ax3.set_ylabel('Water Temperature [degC]')

plt.show()
```
5. Hasil _running script_.

![image](https://user-images.githubusercontent.com/89583653/169735133-878f2913-9057-4fc2-91c5-ba99444edcb1.png)

**Penjelasan**

Dari grafik yang ada, dapat diinterpretasikan hasil bahwa dalam rentang 45 hari model bersifat fluktuatif mulai dari terjadi naik ataupun turunnya parameter hidro-oseanografi yang bekerja. Seperti pada keberadaan kecepatan angin dan nilai tekanan cenderung berubah-ubah setiap harinya. Grafik yang dihasilkan antara angin dan tekanan tidak terlalu jauh berbeda, karena diketahui bahwa angin merupakan udara yang bergerak dikarenakan perbedaan tekanan di permukaan bumi. Sehingga perbedaan tekanan juga mempengaruhi perubahan kecepatan dan arah angin di permukaan laut. Dalam rentang 45 hari pengamatan, yakni bulan April-Mei diketahui memasuki musim Peralihan I dari musim kemarau ke hujan. Apabila dikaitkan dengan variabilitas musiman gelombang dan arus yang bekerja di Perairan Samudera Pasifik bagian Utara maka pengaruh musim kemarau masih mendominasi. Didapatkan pula nilai tertinggi suhu permukaan laut di bulan April-Mei yaitu 26,4 derajat C dan nilai minimum yaitu 25,4 derajat C. Tinggi rendahnya nilai suhu permukaan laut diduga akibat pengaruh pergerakan angin yang bekerja. Suhu minimal terjadi akibat keberadaan angin muson tenggara yang bertiup dan menyebabkan terjadinya Transpor Ekman, sehingga terjadi kekosongan yang berakibat naiknya air (upwelling) yang membawa suhu rendah dari bawah menuju ke lapisan permukaan.


## 4. Penutup✨
Secara terperinci, berdasarkan praktikum Pemodelan Oseanografi sebanyak 4 modul yang telah dilakukan, dapat ditarik ringkasan dan simpulan sebagai berikut:
1. Penggabungan persamaan adveksi-difusi 2D pada dasarnya dapat dipergunakan dalam analisis pergerakan maupun perluasan polutan. Dalam mengetahui pergerakannya, perlu dilakukan diskritisasi terhadap persamaan yang ada dan dimasukkan dalam bahasa pemrograman berupa Phyton maupun Google Colab. Dalam proses perluasan serta penyebaran polutan yang terjadi, diketahui pula faktor penyebab polutan menyebar akibat faktor hidro-oseanografi (arah angin yang bekerja, gelombang dan arus) maupun karakteristik wilayah perairan yang ada.
2. Model hidrodinamika 1D yang dipergunakan untuk memperhitungkan kondisi suatu aliran yang dihasilkan dari distribusi garam, suhu, tunduk pada berbagai kondisi gaya dan batas tertentu. Nantinya akan diketahui hubungan keberadaan waktu yang rendah serta akurasi simulasi dapat dikontrol dalam pemodelan yang digunakan. Mulai dari skema urutan yang kecil sampai lebih tinggi. Nilai serta simulasi skema yang dihasilkan juga dapat mengalami error akibat terlalu banyak data input atau keberadaan nilai parameter yang lebih kompleks. Oleh sebab itu, penggambaran grafik yang dihasilkan tidak selalu mulus. 
3. Model hidrodinamika 2d sendiri dapat diartikan dengan penggambaran maupun tiruan dari suatu fenomena/kejadian/proses pada suatu waktu tertentu. Dengan penggunaan hidrodinamika 2d dapat digunakan untuk memodelkan dan menganalisis suatu fenomena seperti pemodelan gelombang akibat gaya pembangkit angin, pemodelan sampah plastik di laut dan pemodelan coastal dynamics dan sedimentasi pantai. Untuk praktik pemodelan yang telah dilakukan, diketahui hubungan atau korelasi antara kecepatan angin dan tinggi gelombang yang sifatnya linier. Artinya, semakin besar angin yang bekerja (kecepatan dalam knots) pada suatu wilayah periaran maka semakin besar pula gelombang ataupun arus yang akan terbentuk. 

## Ucapan Terima Kasih
Demikian pemenuhan tugas akhir praktikum Pemodelan Oseanografi ini kami buat. 
Seluruh authors memohon maaf apabila masih terdapat kesalahan dan kekurangan secara dalam tugas akhir ini. Tak lupa, Tim 11 selaku _author_ dari _repository_ kali ini juga mengucapkan terimakasih kepada:
1. Seluruh dosen pengampu mata kuliah Pemodelan Oseanografi yang memberikan gambaran umum terkait setiap materi;
    Dr. Aris Ismanto, S.Si, M.Si.
    Prof. Dr. Denny Nugroho Sugianto S.T., M.Si.
    Dr. Elis Indrayanti S.T., M.Si.
    Rikha Widiaratih, S.Si, M.Si. 
2. Seluruh Tim Asisten Praktikum Pemodelan Oseanografi 2022 yang mendampingi dalam praktikum sampai dengan pengerjaan tugas akhir;
    Sekar Adiningsih
    Gamma Haqqul Fikriawan
    Clara Clarita
    Yoas Heryanto
4. Rekan-rekan Oseanografi 2020 yang turut membantu dan mendukung penyelesaian tugas akhir.

Tidak menutup kemungkinan, bahwa _repository_ ini masih jauh dari kata sempurna. Oleh karena itu, kritik dan saran masih sangat diperlukan agar dapat menjadikannya lebih baik. Diharapkan pula, dengan adanya _repository_ ini dapat memberikan kebermanfaatan untuk semua.

Terima kasih.

                                                                                                                 Semarang, 25 Mei 2022
                                                                                                                         Penulis
                                                                                                                                          
                                                                                                                          Tim 11
