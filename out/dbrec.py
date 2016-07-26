import sys
if len(sys.argv) < 9 :
  sys.exit()
N=int(sys.argv[1])
L=int(sys.argv[2])
P=int(sys.argv[3])
Q=int(sys.argv[4])
T=int(sys.argv[5])
F=sys.argv[6]
f1=sys.argv[7]
f2=sys.argv[8]
f =open(f1,'rb')
D=0
FMM=0
E=0

for line in f:
  if '#Data' in line:
    a=line.find(':')
    b=line.find('(',a+5)
    D=int(line[a+1:b])
  if 'FMM' in line:
    a=line.find(':')
    b=line.find('(',a+1)
    FMM=float(line[a+1:b])
  if 'Error norm:' in line:
    flds=line.split()  
    E=float(flds[3])
f.close()

t=[[0,0],
   [0,0],
   [0,0],
   [0,0],
   [0,0],
   [0,0]]
def get_type(s):
  types=['Leaves->Second','Upward','Translation','Downward','Second','NearField']
  for i in range(len(types)):
    if types[i]==s:
      return i
  return -1
f=open(f2,'rb')
i=0

for line in f :
  if line[0] =='a': continue
  if line[0] =='s': continue
  flds=line.split()
  i=get_type(flds[4])
  if i>=0:
    t[i][0]=float(flds[1])
    t[i][1]=int(flds[3])
  #i +=1
f.close()
tT = sum([tt[1] for tt in t])
if 's' in F: T=0
print N,L,P,Q,T,F,FMM,tT,D,E,
for i in range(len(t)):
  print t[i][0],t[i][1],
print
