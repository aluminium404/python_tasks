from scipy.special import spherical_jn as jn #Сферическая функция Бесселя первого рода
from scipy.special import spherical_yn as yn #Сферическая функция Бесселя второго рода
import matplotlib.pyplot as plt
import urllib.request as url 
from numpy import pi 
from re import split 
import numpy as np
import os
# начальные значения для расчета ЭПР
def hn(l, z):
    return jn(l, z) + 1j * yn(l, z)
def an(l, z):
    return jn(l, z) / hn(l, z)
def bn(l, z):
    return (z * jn(l - 1, z) - l * jn(l, z)) \
 / (z * hn(l - 1, z) - l * hn(l, z))
# вывод исходных данных
URL = 'https://jenyay.net/uploads/Student/Modelling/task_02.txt'
file = url.urlopen(URL)
list = file.readlines()
my_string = list[0].decode("utf-8") # принимает первую строку за 0
values=split("; ",my_string)
values1=split("=", values[0])
values2=split("=", values[1])
values3=split("=", values[2])
D=float(values1[1])
fmin=float(values2[1])
fmax=float(values3[1])
print ('D =', D, 'fmin =', fmin, 'fmax =', fmax)
Z = 10000 # число точек на отрезке
r = 0.5 * D
f = np.linspace(fmin, fmax, Z) # указываем количество элементов
lmd = 3e8 / f
k = 2 * pi / lmd
Sum_arr = [((-1) ** n) * (n + 0.5) * (an(n, k * r) - bn(n, k * r)) \
    for n in range(1, 50)] # массив элементов под знаком суммы
Sum = np.sum(Sum_arr, axis=0) # суммирование элементов массива
Sig = (lmd ** 2) / pi * (np.abs(Sum) ** 2) # вычисление ЭПР
plt.plot(f/0.01e9, Sig) # построение графика
plt.xlabel('$f, МГц$')
plt.ylabel('$\sigma, м^2$')
plt.grid() 
plt.show()
# указываем путь

try:
    os.mkdir('results')
except OSError:
    print ("Папка 'results' уже существует")

complete_file = os.path.join('results', 'task_02_307b_Alimov_1.txt')
# преобразуем nparray в list, для облегчения работы с многомерным массивом
fn = f/1e9
fntl = fn.tolist() 
Stl = Sig.tolist()
f = open(complete_file, 'w')
f.write('f, ГГц         Sigma, м^2\n')
#расчёт значений для заданного количества точек
for i in range(Z):
 
    f.write(str("%.2f" % fntl[i])+'           '+str("%.15f" % Stl[i])+"\n")
f.close()

