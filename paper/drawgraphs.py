import pylab as plt
import sys

#0              5                       10                        15                       20
_N,_L,_P,_Q,_T,_FLAGS,_FMM,_nT,_nD,_Err,_aT1,_nT1,_aT2,_nT2,_aT3,_nT3,_aT4,_nT4,_aT5,_nT5,_aT6,_nT6=range(22)
intFlds=[_N,_L,_P,_Q,_T,_nT,_nD,_nT1,_nT2,_nT3,_nT4,_nT5,_nT6]
taskTypes=['Leaves->Second Last','Upward','Translation','Second Last->Leaves','Downward','NearField']
def load_file(fname):
  out = list()
  f=open(fname,'rb')
  for line in f:
    flds=line.split()
    a=list()
    for i in range(len(flds)):
      if i in intFlds:
        v=int(flds[i])
      elif i!= 5:
        v = float(flds[i])
      else:
        v=flds[i]        
      if  i == _T and v>1:
        v /=2
      a.append(v)
    out.append(a)
  f.close()
  return out
#--------------------------------------------------------------
def draw_time_threads_Q(recs,Q,fname):
  rf=filter(lambda x: x[_Q] == Q, recs)
  g1=filter(lambda x: x[_P] == 100, rf)  
  g2=filter(lambda x: x[_P] == 200, rf)
  g3=filter(lambda x: x[_P] == 400, rf)
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  
  y1=[g[_FMM] for g in g1]
  y2=[g[_FMM] for g in g2]  
  y3=[g[_FMM] for g in g3]
  
  lg=[None,None,None]
  plt.figure()
  lg[0],=plt.plot(x1,y1,'r^-',lw=3)
  lg[1],=plt.plot(x2,y2,'gs-',lw=3)
  lg[2],=plt.plot(x3,y3,'bo-',lw=3)
  plt.xlabel('Cores')
  plt.ylabel('Time(s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['P=100','P=200','P=400'])
  plt.grid(True)
  plt.title('FMM Execution time, Q=%d'%Q)
  plt.savefig(fname)
  plt.close()
  

def draw_time_threads_P(recs,P,fname):
  rf=filter(lambda x: x[_P] == P, recs)
  g1=filter(lambda x: x[_Q] == 10, rf)  
  g2=filter(lambda x: x[_Q] == 20, rf)
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  y1=[g[_FMM] for g in g1]
  y2=[g[_FMM] for g in g2]
  
  lg=[None,None]
  plt.figure()
  lg[0],=plt.plot(x1,y1,'r^-',lw=3)
  lg[1],=plt.plot(x2,y2,'gs-',lw=3)
  plt.xlabel('Cores')
  plt.ylabel('Time(s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['Q=10','Q=20'])
  plt.grid(True)
  plt.title('FMM Execution time, P=%d'%P)
  plt.savefig(fname)
  plt.close()
  
#------------------------------------------------
def get_max_time(recs):
  if len(recs)==0:
    return 1
  return max (recs)
def draw_speedup_Q(recs,Q,fname):
  rf=filter(lambda x: x[_Q] == Q, recs)
  print "----- SpeedUp Q"
  for r in rf:
    print r
  
  g0=filter(lambda x: x[_P] ==  50, rf)
  print "====="
  for g in g0:
    print g
  g1=filter(lambda x: x[_P] == 100, rf)  
  g2=filter(lambda x: x[_P] == 200, rf)
  g3=filter(lambda x: x[_P] == 300, rf)
  g4=filter(lambda x: x[_P] == 400, rf)
  x0=[g[_T] for g in g0]
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  x4=[g[_T] for g in g4]
  
  t0=get_max_time([g[_FMM] for g in g0 if g[_T]==1])
  t1=get_max_time([g[_FMM] for g in g1 if g[_T]==1])
  t2=get_max_time([g[_FMM] for g in g2 if g[_T]==1])
  t3=get_max_time([g[_FMM] for g in g3 if g[_T]==1])
  t4=get_max_time([g[_FMM] for g in g4 if g[_T]==1])
  
  y0=[t0/g[_FMM] for g in g0]
  y1=[t1/g[_FMM] for g in g1]
  y2=[t2/g[_FMM] for g in g2]  
  y3=[t3/g[_FMM] for g in g3]
  y4=[t4/g[_FMM] for g in g4]
  
  lg=[None,None,None,None,None,None]
  plt.figure()
  lg[0],=plt.plot(x0,y0,'r^-',lw=3)
  lg[1],=plt.plot(x1,y1,'mv-',lw=3)
  lg[2],=plt.plot(x2,y2,'gs-',lw=3)
  lg[3],=plt.plot(x3,y3,'co-',lw=3)
  lg[4],=plt.plot(x4,y4,'bo-',lw=3)
  lg[5],=plt.plot([0,1,2,4,8],[1,1,2,4,8],'k-',lw=3)
  plt.plot([1,1],[0,12],'k--',lw=3)
  plt.xlabel('Threads')
  plt.ylabel('Speedup')
  plt.ylim([0,10])
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'4', '8','16'])
  plt.legend(lg,['P=50','P=100','P=200','P=300','P=400','Ideal'],loc='upper left')
  plt.grid(True)
  plt.title('N= %d'%(g1[0][_N]))
  
  plt.savefig(fname)
  plt.close()
  
#------------------------------------------------
def draw_speedup_P(recs,P,fname):
  rf=filter(lambda x: x[_P] == P, recs)
  print "----- SpeedUp P"
  for r in rf:
    print r
  g1=filter(lambda x: x[_Q] == 10, rf)  
  g2=filter(lambda x: x[_Q] == 20, rf)
  g3=filter(lambda x: x[_Q] == 50, rf)
  g4=filter(lambda x: x[_Q] ==100, rf)
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  x4=[g[_T] for g in g4]
  
  t1=get_max_time([g[_FMM] for g in g1 if g[_T] == 1])
  t2=get_max_time([g[_FMM] for g in g2 if g[_T] == 1])
  t3=get_max_time([g[_FMM] for g in g3 if g[_T] == 1])
  t4=get_max_time([g[_FMM] for g in g4 if g[_T] == 1])
  
  y1=[t1/g[_FMM] for g in g1]
  y2=[t2/g[_FMM] for g in g2]  
  y3=[t3/g[_FMM] for g in g3]
  y4=[t4/g[_FMM] for g in g4]
  
  lg=[None,None,None,None,None,None]
  plt.figure()
  lg[0],=plt.plot(x1,y1,'r^-',lw=3)
  lg[1],=plt.plot(x2,y2,'gs-',lw=3)
  lg[2],=plt.plot(x3,y3,'bo-',lw=3)
  lg[3],=plt.plot(x4,y4,'mv-',lw=3) 
  lg[4],=plt.plot([0,1,2,4,8],[1,1,2,4,8],'k-',lw=3)
  lg[5],=plt.plot([1,1],[0,10],'k--',lw=3)
  
  plt.xlabel('Threads')
  plt.ylabel('Speedup')
  plt.ylim([0,10])
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'4', '8','16'])
  plt.legend(lg,['Q=10','Q=20','Q=50','Q=100','Ideal',''],loc='upper left')
  plt.grid(True)
  plt.title('N= %d'%(g1[0][_N]))
  
  plt.savefig(fname)
  plt.close()
#------------------------------------------------
def draw_speedup_merge(recs,fname):
  rf1=filter(lambda x: x[_P] ==  50, recs)
  for r in rf1:
    print r
  rf2=filter(lambda x: x[_P] == 300, recs)
  g1=filter(lambda x: x[_Q] == 10, rf1)  
  g2=filter(lambda x: x[_Q] == 20, rf1)
  g3=filter(lambda x: x[_Q] == 50, rf1)
  g4=filter(lambda x: x[_Q] ==100, rf1)

  g12=filter(lambda x: x[_Q] == 10, rf2)  
  g22=filter(lambda x: x[_Q] == 20, rf2)
  g32=filter(lambda x: x[_Q] == 50, rf2)
  g42=filter(lambda x: x[_Q] ==100, rf2)

  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  x4=[g[_T] for g in g4]

  x12=[g[_T] for g in g12]
  x22=[g[_T] for g in g22]
  x32=[g[_T] for g in g32]
  x42=[g[_T] for g in g42]

  
  t1=get_max_time([g[_FMM] for g in g1 if g[_T] == 1])
  t2=get_max_time([g[_FMM] for g in g2 if g[_T] == 1])
  t3=get_max_time([g[_FMM] for g in g3 if g[_T] == 1])
  t4=get_max_time([g[_FMM] for g in g4 if g[_T] == 1])
  
  t12=get_max_time([g[_FMM] for g in g12 if g[_T] == 1])
  t22=get_max_time([g[_FMM] for g in g22 if g[_T] == 1])
  t32=get_max_time([g[_FMM] for g in g32 if g[_T] == 1])
  t42=get_max_time([g[_FMM] for g in g42 if g[_T] == 1])

  y1=[t1/g[_FMM] for g in g1]
  y2=[t2/g[_FMM] for g in g2]  
  y3=[t3/g[_FMM] for g in g3]
  y4=[t4/g[_FMM] for g in g4]

  y12=[t12/g[_FMM] for g in g12]
  y22=[t22/g[_FMM] for g in g22]  
  y32=[t32/g[_FMM] for g in g32]
  y42=[t42/g[_FMM] for g in g42]
  
  lg=[None,None,None,None,None,None,None,None,None,None]
  lbls=['Q=10, P=50' ,'Q=20, P=50' ,'Q=50, P=50' ,'Q=100, P=50' ,
        'Q=10, P=300','Q=20, P=300','Q=50, P=300','Q=100, P=300',
        'Ideal','']  
  fig=plt.figure()
  ax = fig.add_subplot(111)

  lg[4],=ax.plot(x12,y12,'r^-',lw=3,label=lbls[4])
  lg[0],=ax.plot(x1,y1,'r^--',lw=3,label=lbls[0])

  lg[6],=ax.plot(x32,y32,'bo-',lw=3,label=lbls[6])
  lg[2],=ax.plot(x3,y3,'bo--',lw=3,label=lbls[2])

  lg[7],=ax.plot(x42,y42,'mv-',lw=3,label=lbls[7]) 
  lg[3],=ax.plot(x4,y4,'mv--',lw=3,label=lbls[3]) 

  lg[1],=ax.plot(x2,y2,'gs--',lw=3,label=lbls[1])
  lg[5],=ax.plot(x22,y22,'gs-',lw=3,label=lbls[5])

  lg[8],=ax.plot([0,1,2,4,8],[1,1,2,4,8],'k-',lw=3,label=lbls[8])
  #lg[9],=ax.plot([1,1],[0,10],'k--',lw=3,label=lbls[9])
  
  plt.xlabel('Threads')
  plt.ylabel('Speedup')
  plt.ylim([0,4])
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'4', '8','16'])
  box = ax.get_position()
  ax.set_position([box.x0, box.y0+.2, box.width * 0.8, box.height*.8])
  
  plt.legend(bbox_to_anchor=(-0.01,-.4,1.03,.1), loc='lower left',
             ncol=3, mode="expand")  
  plt.grid(True)
  fig.subplots_adjust(top=.9,bottom=.25,hspace=.75)
  ax.set_title('N=%d'%(g1[0][_N]))
  
  plt.savefig(fname)
  plt.close()
def draw_speedup_mixed(rf,fname):  
  print "----- SpeedUp P and Q "
  for r in rf:
    print r
  g1=filter(lambda x: x[_Q] == 10 and x[_P]== 50, rf)  
  g2=filter(lambda x: x[_Q] ==100 and x[_P]== 300, rf)
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  
  t1=get_max_time([g[_FMM] for g in g1 if g[_T] == 1])
  t2=get_max_time([g[_FMM] for g in g2 if g[_T] == 1])
  
  y1=[t1/g[_FMM] for g in g1]
  y2=[t2/g[_FMM] for g in g2]  
  
  lg=[None,None,None]
  plt.figure()
  lg[0],=plt.plot(x1,y1,'r^-',lw=3)
  lg[1],=plt.plot(x2,y2,'bs-',lw=3)
  lg[2],=plt.plot([0,1,2,4,8],[1,1,2,4,8],'k-',lw=3)
  plt.plot([1,1],[0,12],'k--',lw=3)
  
  plt.xlabel('Threads')
  plt.ylabel('Speedup')
  plt.ylim([0,8])
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'4', '8','16'])
  plt.legend(lg,['Q=10, P=50','Q=100, P=300','Ideal'],loc='upper left')
  plt.grid(True)
  plt.title('N=%d'%(g1[0][_N]))
  
  plt.savefig(fname)
  plt.show()
  plt.close()
#------------------------------------------------
def draw_speedup_far_field(recs,Q,P,fname):
  rf=filter(lambda x: x[_Q] == Q and x[_P] == P and x[_T] != 0, recs)
  print " ---- "
  for r in rf:
    print r
  r1=[r[_aT1]*r[_nT1] for r in rf]
  r2=[r[_aT2]*r[_nT2] for r in rf]
  r3=[r[_aT3]*r[_nT3] for r in rf]
  r4=[r[_aT4]*r[_nT4] for r in rf]
  r5=[r[_aT5]*r[_nT5] for r in rf]
  
  x=[r[_T] for r in rf]  
  
  y1=[r1[0]/r for r in r1]
  y2=[r2[0]/r for r in r2]
  y3=[r3[0]/r for r in r3]
  y4=[r4[0]/r for r in r4]
  y5=[r5[0]/r for r in r5]
  
  lg=[None,None,None,None,None,None,None]
  plt.figure()
  lg[0],=plt.plot(x,y1,'r^-',lw=3)
  lg[1],=plt.plot(x,y2,'gs-',lw=3)
  lg[2],=plt.plot(x,y3,'bo-',lw=3)
  lg[3],=plt.plot(x,y4,'c*-',lw=3)
  lg[4],=plt.plot(x,y5,'mv-',lw=3)
  lg[5],=plt.plot([1,2,4],[1,2,4],'k-',lw=3)
  lg[6],=plt.plot([1,2,4],[1,4,8],'k--',lw=3)
  plt.xlabel('Cores')
  plt.ylabel('Speed up')
  plt.xticks([1,2,4,8],['1' ,'2', '4','8'])
  plt.legend(lg,[taskTypes[0],taskTypes[1],taskTypes[2],taskTypes[3],taskTypes[4],'Ideal (FPU)','Ideal (CPU)'],loc='upper left')
  plt.grid(True)
  plt.title('Speed up , N= %d, Q=%d, P=%d'%(rf[0][_N],Q,P))
  
  plt.savefig(fname)  
  plt.close()
  
#------------------------------------------------
def draw_time_tasks_P(recs,P,fname):
  rf=filter(lambda x: x[_P] == P, recs)
  g1=filter(lambda x: x[_Q] == 10, rf)
  g2=filter(lambda x: x[_Q] == 20, rf)
  f1=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g1 ]
  n1=[g[_aT6] for g in g1]
  f2=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g2 ]
  n2=[g[_aT6] for g in g2]
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  
  lg=[None,None,None,None]  
  plt.figure()
  lg[0],=plt.plot(x1,n1,'rs-',lw=3)
  lg[1],=plt.plot(x2,n2,'rs--',lw=3)
  lg[2],=plt.plot(x1,f1,'bo-',lw=3)
  lg[3],=plt.plot(x2,f2,'bo--',lw=3)
  
  plt.xlabel('Cores')
  plt.ylabel('Time($\mu$s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['Near Field, Q=10','Near Field, Q=20','Far Field, Q=10','Far Field, Q=20'],loc='upper left')
  plt.grid(True)
  plt.title('Tasks Execution times, P=%d'%P)

  plt.savefig(fname)
  plt.close()
  
#------------------------------------------------
def draw_time_tasks_P_far(recs,P,fname):
  rf=filter(lambda x: x[_P] == P, recs)
  g1=filter(lambda x: x[_Q] == 10, rf)
  g2=filter(lambda x: x[_Q] == 20, rf)
  f1=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g1 ]
  f2=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g2 ]
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  
  lg=[None,None]  
  plt.figure()
  lg[0],=plt.plot(x1,f1,'bo-',lw=3)
  lg[1],=plt.plot(x2,f2,'bo--',lw=3)
  
  plt.xlabel('Cores')
  plt.ylabel('Time($\mu$s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['Far Field, Q=10','Far Field, Q=20'],loc='upper left')
  plt.grid(True)
  plt.title('Tasks Execution times for far fields, P=%d'%P)

  plt.savefig(fname)
  plt.close()
  
#-------------------------------------------------------------------
def draw_time_tasks_Q(recs,Q,fname):
  rf=filter(lambda x: x[_Q] == Q, recs)
  g1=filter(lambda x: x[_P] == 100, rf)
  g2=filter(lambda x: x[_P] == 200, rf)
  g3=filter(lambda x: x[_P] == 400, rf)
  f1=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g1 ]
  n1=[g[_aT6] for g in g1]
  f2=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g2 ]
  n2=[g[_aT6] for g in g2]
  f3=[(g[_aT1]+g[_aT2]+g[_aT3]+g[_aT4]+g[_aT5]) for g in g3 ]
  n3=[g[_aT6] for g in g3]
  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  
  lg=[None,None,None,None,None,None]  
  plt.figure()
  lg[0],=plt.plot(x1,n1,'rs-',lw=3)
  lg[1],=plt.plot(x2,n2,'rs--',lw=3)
  lg[2],=plt.plot(x3,n3,'rs:',lw=3)
  lg[3],=plt.plot(x1,f1,'bo-',lw=3)
  lg[4],=plt.plot(x2,f2,'bo--',lw=3)
  lg[5],=plt.plot(x3,f3,'bo:',lw=3)
  
  plt.xlabel('Cores')
  plt.ylabel('Time($\mu$s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['Near Field, P=100','Near Field, P=200','Near Field, P=400','Far Field, P=100','Far Field, P=200','Far Field, P=400'],loc='upper left')
  plt.grid(True)
  plt.title('Tasks Execution times, Q=%d'%Q)

  plt.savefig(fname)
  plt.close()
  
#------------------------------------------------
def draw_time_tasks_Q_near(recs,Q,fname):
  rf=filter(lambda x: x[_Q] == Q, recs)
  g1=filter(lambda x: x[_P] == 100, rf)
  g2=filter(lambda x: x[_P] == 200, rf)
  g3=filter(lambda x: x[_P] == 400, rf)
  n1=[g[_aT6] for g in g1]
  n2=[g[_aT6] for g in g2]
  n3=[g[_aT6] for g in g3]

  x1=[g[_T] for g in g1]
  x2=[g[_T] for g in g2]
  x3=[g[_T] for g in g3]
  
  lg=[None,None,None]  
  plt.figure()
  lg[0],=plt.plot(x1,n1,'rs-',lw=3)
  lg[1],=plt.plot(x2,n2,'rs--',lw=3)
  lg[2],=plt.plot(x3,n3,'rs:',lw=3)
  
  plt.xlabel('Cores')
  plt.ylabel('Time($\mu$s)')
  plt.xticks([0,1,2,4,8],['Sequential', '1' ,'2', '4','8'])
  plt.legend(lg,['Near Field, P=100','Near Field, P=200','Near Field, P=400'],loc='upper left')
  plt.grid(True)
  plt.title('Tasks Execution times for near fields, Q=%d'%Q)

  plt.savefig(fname)
  plt.close()
#--------------------------------------------------------------------------------  
def group_N(records,N):
  recs = filter(lambda x: x[_N] == N,records)
  if len(recs) ==0: return  
  for r in recs:print r
  if False:
    draw_time_threads_Q(recs,10,'N%d_Q10.pdf'%N)
    draw_time_threads_Q(recs,20,'N%d_Q20.pdf'%N)
    
    draw_time_threads_P(recs,100,'N%d_P100.pdf'%N)
    draw_time_threads_P(recs,200,'N%d_P200.pdf'%N)
    draw_time_threads_P(recs,400,'N%d_P400.pdf'%N)
  
    draw_time_tasks_P(recs,100,'N%d_tasks_P100.pdf'%N)
    draw_time_tasks_P(recs,200,'N%d_tasks_P200.pdf'%N)
    draw_time_tasks_P(recs,400,'N%d_tasks_P400.pdf'%N)
    
    draw_time_tasks_Q(recs,10,'N%d_tasks_Q10.pdf'%N)
    draw_time_tasks_Q(recs,20,'N%d_tasks_Q20.pdf'%N)
    
    draw_time_tasks_P_far(recs,100,'N%d_tasks_far_only_P100.pdf'%N)
    draw_time_tasks_P_far(recs,200,'N%d_tasks_far_only_P200.pdf'%N)
    draw_time_tasks_P_far(recs,400,'N%d_tasks_far_only_P400.pdf'%N)
    
    draw_time_tasks_Q_near(recs,10,'N%d_tasks_near_only_Q10.pdf'%N)
    draw_time_tasks_Q_near(recs,20,'N%d_tasks_near_only_Q20.pdf'%N)
  
  
    recsF=filter(lambda x: 'ftwSF' == x[_FLAGS],recs)  
    draw_speedup_merge(recsF,'run_far_field_only_N%d_speedup_merged.pdf'%N)
    recsNF=filter(lambda x: 'ftwSNF' == x[_FLAGS],recs)  
    draw_speedup_mixed(recsNF,'run_mixed_N%d_speedup.pdf'%N)
  recsN=filter(lambda x: 'ftwSN' == x[_FLAGS],recs)  
  draw_speedup_Q(recsN,10 ,'run_near_field_only_N%d_speedup_Q10.pdf'%N)
  draw_speedup_Q(recsN,100 ,'run_near_field_only_N%d_speedup_Q100.pdf'%N)
fname=sys.argv[1]
recs=load_file(fname)
recs=filter(lambda x: x[_FMM] != 0 , recs)
recs=sorted(recs,key=lambda x:x[_T])
group_N(recs,1e5)
#group_N(recs,1e6)
#group_N(recs,2e6)
#plt.show()
print "Finished"