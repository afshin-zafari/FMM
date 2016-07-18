import sys,math
def load_result(fname):
  f = open(fname,'rb')
  y  = list()
  for v in f:
    a = v.find(',')
    b = v[0:a]
    y.append(float(b))
  
  f.close()
  return y

f1 = sys.argv[1]
f2 = sys.argv[2]

y1=load_result(f1)
y2=load_result(f2)
print len(y1),len(y2)
e=0;n1=0;n2=0
for i in range(len(y1)):
  e += (y1[i] - y2[i] ) **2
  n1+=y1[i]**2
  n2+=y2[i]**2

es =math.sqrt(e)
n1s=math.sqrt(n1)
n2s=math.sqrt(n2)
print "Error norm: ", es,es/n1s,es/n2s