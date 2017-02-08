#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

#include <dirent.h>

#include "util.hpp"
#ifndef OMP_TASKS
  #include "sg/superglue.hpp"
  #include "sg/option/instr_trace.hpp"
  SuperGlue<Options> *sgEngine;
#endif //OMP_TASKS

Time::TimeUnit exTime;
FMM::Matrix   *pts_ptr,*c_p,*q_p;
FMM::SGMatrix *Y,*C_p,*Q_p;
FMM::Tree     *OT_ptr;
//FMM::MemoryPool *pool;
//FMM::Config fmm_config;
//int N,L;


/*---------------------------------------------------------------*/
void init(FMM::Matrix &pts){
  int M = FMM::N;
  double theta = 0.0, theta_mod =0.0;
  for ( int i=1; i<=M;i++){
    pts(i,1) = 2* cos(theta)*(1.0 + .1*cos(theta_mod));
    pts(i,2) = 1* sin(theta)*(1.0 + .1*cos(theta_mod));
    pts(i,3)=0.0;
    theta += 2.*M_PI /FMM::N;
    theta_mod += 32*M_PI /FMM::N;
  }
}
/*---------------------------------------------------------------*/
void fmm_solver(){
  OT_ptr =new FMM::Tree;
  FMM::Tree &OT=*OT_ptr;    
  pts_ptr= new FMM::Matrix (FMM::N,3);
  FMM::Matrix &pts=*pts_ptr;

  FMM::Reader rdr(FMM::config.tree,FMM::config.ops,OT);
  fprintf(stdout,"Reading tree...\n");
  rdr.read();
  fprintf(stdout,"Reading operators ...\n");
  rdr.read_op();
  OT.Q = FMM::config.Q;
  fprintf(stdout,"Initializing points...\n");
  init(pts);

  fprintf(stdout,"Compute NearFields...\n");
  FMM::compute_near_field(OT,pts);
    
	
  c_p = new FMM::Matrix (FMM::N,1,1.0);
  FMM::Matrix &c = * c_p;
    
  q_p = new FMM::Matrix (FMM::N,1,0.0);
  FMM::Matrix &q = * q_p;

    
  C_p =new FMM::SGMatrix(c);
  FMM::SGMatrix &C = *C_p;
    
  Q_p=new FMM::SGMatrix(q);
  FMM::SGMatrix &Q = *Q_p;
  Y=&Q;
  exTime= Time::getTime();
    
  cout << Time::getTime() << endl;
      timeval tv;
    int e=gettimeofday(&tv,NULL);
    cout << "TOD: " << tv.tv_sec <<"," << tv.tv_usec << "," << e << endl;
    #pragma omp parallel 
    {
      #pragma omp single
      {
	fprintf(stdout,"MatVec near field...\n");
	if(FMM::config.NF)
	  FMM::mv_near_field(OT,C,Q);            
      
	fprintf(stdout,"MatVec far field...\n");
	if(FMM::config.FF)
	  FMM::MatVec(OT,C,Q);
      }
    }
#pragma omp taskwait 
#pragma omp barrier 
  cout << Time::getTime() << endl;
     e=gettimeofday(&tv,NULL);
			     cout << "TOD: " << tv.tv_sec <<"," << tv.tv_usec << "," << e << endl;
#ifndef OMP_TASKS
  if (FMM::config.t){
    sgEngine->barrier(); // Wait until all tasks finished
  }
  double execTime = ((Time::getTime() - exTime)*1.0)/3000000000.0;
#else  
  double execTime = ((Time::getTime() - exTime)*1.0);
#endif  
  cout << " Program Finished. Time(s): " <<  execTime << endl;
    
}
/*---------------------------------------------------------------*/
int main(int argc , char *argv[])
{
  if ( argc <9){
    fprintf(stderr,"Usage %s N L Q S T tree_file ops_file nfstSOwa \n",argv[0]);
    exit(-1);
  }
  FMM::parse_args(argc,argv);
  FMM::stats.t = 0;
  FMM::pool = new FMM::MemoryPool(1e9);
#ifndef OMP_TASKS
  if(FMM::config.t)
    sgEngine = new SuperGlue<Options>(FMM::config.cores);
#endif

  fmm_solver();

  char res[25];
  sprintf(res,"results_%c_%c_%c_%c.txt",
	  FMM::config.n?'n':'f',
	  FMM::config.s?'s':'t',
	  FMM::config.S?'S':'O',
	  FMM::config.w?'w':'a'
	  );
	
  Y->get_matrix()->export_data(res);
    
  char trace[25];
  sprintf(trace,"trace_%c_%c_%c_%c.txt",
	  FMM::config.n?'n':'f',
	  FMM::config.s?'s':'t',
	  FMM::config.S?'S':'O',
	  FMM::config.w?'w':'a'
	  );
#ifndef OMP_TASKS
  Trace<Options>::dump(trace);
#endif
  delete FMM::pool;
}
