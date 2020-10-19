import numpy as np
import matplotlib.pyplot as plt
import os.path

A = 10
x = np.linspace(-5.12,5.12,200)

def f(x):
    return A + x**2 - A*np.cos(2*np.pi*x)

y = f(x)

print(y)

if not os.path.exists('results'):
    os.mkdir('results')

with open("results/task_01_307b_Alimov_1.txt",'w') as txt:
    
    i = 0
    
    for _x in x:
        
        txt.write(str(_x)+'    '+str(y[i])+'\n')
        
        i+=1


plt.plot(x,y)
plt.grid()
plt.xlabel("x")
plt.ylabel("f(x)")
plt.show()
