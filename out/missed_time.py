import sys
def load_file(fname):  
  f=open(fname,'rb')
  for line in f:
    if 'Mv_near' in line:
      a=line.find(':')
      b=line.find('(',a+1)
      mvnear=float(line[a+1:b])
    if 'Mv_far' in line:
      a=line.find(':')
      b=line.find('(',a+1)
      mvfar=float(line[a+1:b])
    if 'Barrier' in line:
      a=line.find(':')
      b=line.find('(',a+1)
      if b==-1:continue
      barrier=float(line[a+1:b])
    if 'FMM' in line:
      a=line.find(':')
      b=line.find('(',a+1)
      fmm=float(line[a+1:b])
  f.close()  
  return mvnear,mvfar,barrier,fmm

fname=sys.argv[1]
n,f,b,m=load_file(fname)
print
print "Near+Far+Barrier: ",n+f+b
print "Diff: ",m-b-f-n