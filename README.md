# Tugas Akhir Prk.Pemodelan Oseanografi (Team 11)
Repositori ini dibuat untuk memenuhi Tugas Akhir Praktikum Pemodelan Oseanografi 2022. _Repository_ ini memuat file berupa _script phyton (py)_ yang dapat memproses beberapa pemodelan oseanografi seperti Adveksi-Difusi dan Hidrodinamika. Pengerjaan untuk _repository_ kali ini, menggunakan bahasa pemrograman _python_ yang dapat dilakukan pada beberapa platform seperti _Google Colaboratory_ dan _Jupyter Notebook_. Sedangkan untuk _library_ yang digunakan kali ini adalah _Numpy, Matplotlib, IPython, Scipy,_ dan _Pprint_. Seluruh _script_ yang dibuat adalah hasil Team 11 Oseanografi 2020. Semoga dapat bermanfaat!
## 1. AUTHORS (TEAM 11)
1. Ariestia Finola Damanik 26050120120003
2. Barriel Jimly Al Ariefanzah 26050120140122
3. Dimas Sukma Hadi 26050120140170
4. Fahri Rahmalia 26050120130037
5. Nisa Uswatun Khasanah 26050120130119
6. Pilipus Najuita Nduru 26050120120008
7. Yvetty Zilla Durhan 26050120120021

## 2. Cara Penggunaan _Script_
1. Pengguna dapat membuka folder **Kumpulan Script** pada repository ini
2. Di dalamnya terdapat 3 _script_, bisa dibuka salah satu
![image](https://user-images.githubusercontent.com/89583653/169198940-5e831637-1628-4be3-b0be-a87c5226ee7f.png)
3. Nanti akan muncul _script_-nya seperti gambar di bawah ini
![image](https://user-images.githubusercontent.com/89583653/169198974-f7dc70ad-7edf-4b67-85f1-3685a3914558.png)
4. _Script_ dapat di-_copy_ untuk nantinya dipergunakan dan disesuaikan di _Jupyter Notebook_

## 3. Metode Pengerjaan
1. Adveksi-Difusi 2D
2. Hidrodinamika 1D
3. Hidrodinamika 2D

### 3.1. Adveksi-Difusi 2D
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
### 3.2. Hidrodinamika 1D
-
-
-

### 3.3. Hidrodinamika 2D
(ini aku cuma nyoba-nyoba ya)
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


-

## 4. Penutup
