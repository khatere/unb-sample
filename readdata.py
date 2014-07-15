import os
import csv
import numpy as np
import matplotlib.pyplot as plt


#Helper function

def ourFit(x,y,deg,xi= None):
 if xi== None:
     xi =x
 fitCoeff = np.polyfit(x, y, deg)
 fitPoly = np.poly1d(fitCoeff)
 return fitPoly(xi)


dirName = '\Users\Khatere\Documents\sampython'
fileName = 'randomdata.dat'

with open(os.path.join(dirName,fileName),'r') as csvfile:
    dataReader = csv.reader(csvfile, delimiter = ',')
    data= list()
    for row in dataReader:
        data.append(row)
        
    #to get rid of the the column titles (x,y,z) to be able to work on the data as numbers
    print data
    data.pop(0)
    print data
    
    #data = data[1:]
    #print data

    
    x = [int(data[i][0]) for i in range(len(data))]
    y = [float(data[i][1]) for i in range(len(data))]
    z=  [float(data[i][2]) for i in range(len(data))]
#q = [y[i]*z[i] for i in range(len(data))]
#print q

deg = 2
yfitValues = ourFit(x, y, deg)
xrandom= [4,5,6,5,8,5]
yfitrandom = ourFit(x,y,deg,xrandom)

#yFitCoeff = np.polyfit(x,y,deg)
#yFitPoly = np.poly1d(yFitCoeff)
#yFitValues = yFitPoly(x)

#using the Helper function instead
yFitValues = ourFit(x,y,deg)


#print yFitCoeff
#print yFitPoly
#print yFitPoly(yFitValues)
plt.figure()
plt.plot(x,y,'o')
plt.plot(x,yFitValues)
plt.xlim(0,10)
plt.show()



