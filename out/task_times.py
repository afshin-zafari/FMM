import sys
def load_file(fname):
  out=list()
  fl = open(fname,'rb')
  for line in fl:
    fld = line.split(',')
    if len(fld) != 5: continue
    s=float(fld[0])/3e9
    r=float(fld[1])/3e9
    f=float(fld[2])/3e9    
    d=float(fld[3])/3e9
    t=float(fld[4])
    out.append([r-s,f-r,d-f,t])
  fl.close()
  return out
def avg(l):
  if len(l) ==0 : return 0,0
  return sum(l) / len(l)*1e6,len(l)
recs=load_file(sys.argv[1])
tasks=['','Leaves->Second Last','Upward','Translation','Second Last->Leaves','--','Downward','NearField','--']
print "average distance between tasks' events in microseconds"
print 'submit-run \t run-finish \t finish-destr \t   count \t type '
for type in range(1,9):
  t=filter(lambda x: x[-1] == type,recs)  
  t1=avg([tt[0] for tt in t])
  if t1[0]==0:continue
  t2=avg([tt[1] for tt in t])
  t3=avg([tt[2] for tt in t])
  print '%9.2f \t %9.2f \t %9.2f \t %7d \t %s' %(t1[0],t2[0],t3[0],t1[1],tasks[type])