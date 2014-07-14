
def dostuff(x):
    y = [i**2+4 for i in x]
    return y 

x1 = range(1,4)
x2 = range (3,7)

y1 = [i**2 + 4 for i in x1]
y2 = [i**2 + 4 for i in x2]
z1 = dostuff(x1)
z2 = dostuff(x2)
print x1
print y1
print y2
print z1
print z2