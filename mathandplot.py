import os
import csv
import random as rand
import matplotlib.pyplot as plt
x = range(1,10)
y = [(i+rand.random())**2 for i in x]
z = [(i+rand.random()*2)**2 for i in x]
print y

plt.figure()
font ={'family':'Computer Modern Roman', 'size':16}
plt.rc('font', **font)
plt.rc('text', usetex= True)
plt.plot(x,y, 'o')
plt.plot(x,z, 'o')
plt.title('y and z versus x' )
plt.xlabel('x sec')
plt.ylabel('x**2')
plt.show()
dirName = '\Users\Khatere\Documents\sampython'
fileName = 'randomdata.dat'

if not os.path.exists(dirName):
    os.makedirs(dirName)
    
dataList = list()
[dataList.append([x[i] ,y[i], z[i]]) for i in range(len(x))]
print dataList
with open(os.path.join(dirName, fileName),'w') as csvfile:
    writer = csv.writer(csvfile, delimiter = ',')
    writer.writerow(['x', 'y', 'z'])
    writer.writerows(dataList)