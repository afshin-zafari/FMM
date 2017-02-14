#include "types.hpp"

namespace FMM{
  /*----------------------------------------------------------------------------*/
  extern Config config;
  Statistics stats;
  vector<SGTaskGemv *> tlist;
  /*----------------------------------------------------------------------------*/
  void cblas_dgemv(const int layout,
		   const bool TransA,
		   const int M, const int N,
		   const double alpha, const double *A, const int lda,
		   const double *X, const int incX,
		   const double beta, double *Y, const int incY){
    if (!TransA){
      for ( int i =0 ; i< M ; i++){
	for ( int j=0;j<N;j++){
	  Y[i] = beta * Y[i] + alpha * A[j*lda+i] * X[j];
	}
      }
    }
    else{
      for ( int j=0; j<M ; j++){
	for ( int i=0; i<N ; i++){
	  Y[i] = beta * Y[i] + alpha * A[i*lda+j] * X[j];
	}
      }
    }
  }
  /*----------------------------------------------------------------------------*/
  void submit(SGTaskGemv *t){
    if ( config.omp){
      t->run();
      return;
    }
    cout << "OMP not used." << endl;
    if (config.s and !config.t){
      t->run();
      return;
    }
    #ifndef OMP_TASKS  
      stats.t++;
      if (config.l){
	tlist.push_back(t);
      }
      else {
	Time::TimeUnit tt = Time::getTime();
	t->s = tt;
	sgEngine->submit(t);
	tt = Time::getTime() - tt;
	submit_time += tt;
      }
    #endif
  }
  /*----------------------------------------------------------------------------*/
  void submit_all(){
#ifndef OMP_TASKS
    for ( auto t: tlist){
      sgEngine->submit(t);
    }
#endif
    fprintf(stdout,"submit all called.\n");
  }
  /*----------------------------------------------------------------------------*/
  void gemv_leaves(SGMatrix &V, SGMatrix &x, SGMatrix &S, int group){
    SGTaskGemv *t= new SGTaskGemv(V,x,group,S,group);
    assert(V.get_matrix());
    assert(S.get_part(group).get_matrix());
    assert(x.get_part(group).get_matrix());
    t->type=1;
    t->set_name("Vx->s");
    submit(t);
  }
  /*----------------------------------------------------------------------------*/
  void gemv_upward(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
    assert(M2M.get_matrix());
    assert(v.get_part(i1).get_matrix());
    assert(y.get_part(i2).get_matrix());
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
    t->type=2;
    t->set_name("Ms->s");
    submit(t);
  }
  /*----------------------------------------------------------------------------*/
  void gemv_translation(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
    assert(M2M.get_matrix());
    assert(v.get_part(i1).get_matrix());
    assert(y.get_part(i2).get_matrix());
    Time::TimeUnit tt = Time::getTime();
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
    tt = Time::getTime() - tt;
    stats.dur[0] += tt;
    t->type=3;
    t->set_name("Ms->o");
    submit(t);

  }
  /*----------------------------------------------------------------------------*/
  void gemv_final(SGMatrix &A, SGMatrix &v, int i, SGMatrix &y){
    assert(A.get_matrix());
    assert(v.get_part(i).get_matrix());
    assert(y.get_part(i).get_matrix());
    SGTaskGemv *t= new SGTaskGemv(A,v,i,y,i);
    t->transA = true;
    t->type=4;
    t->set_name("V^To->y");
    submit(t);    
  }
  /*----------------------------------------------------------------------------*/
  void gemv_downward(SGMatrix &A, SGMatrix &x, SGMatrix &y,int i ){
    assert(A.get_matrix());
    assert(x.get_part(i).get_matrix());
    assert(y.get_part(i).get_matrix());
    y.get_part(i).get_matrix()->print();
    SGTaskGemv *t= new SGTaskGemv(A,x,i,y,i);
    t->type=5;
    t->set_name("Vx->S");
    submit(t);
  }
  /*----------------------------------------------------------------------------*/
  void gemv_down_obs(SGMatrix &A, SGMatrix &v, int i1, SGMatrix &y, int i2){
    assert(A.get_matrix());
    assert(v.get_part(i1).get_matrix());
    assert(y.get_part(i2).get_matrix());
    SGTaskGemv *t= new SGTaskGemv(A,v,i1,y,i2);
    t->transA = true;
    t->type=6;
    t->set_name("M^To->o");
    submit(t);
  }
  /*----------------------------------------------------------------------------*/
  void MatVec(  Tree & OT,SGMatrix &x , SGMatrix &y){

    int nLev = length(OT);
    int Q = OT.Q ;
    bool dbg=!true;

    EventLog ev("Init.");

    // INIT
    for (int iLev = 3;iLev<= nLev;iLev++){
      int nGroups = length(OT(iLev).group) ;
      Matrix *S = new Matrix (Q,nGroups);
      Matrix *O = new Matrix (Q,nGroups);
      double *sm = S->get_data_memory();
      double *om = O->get_data_memory();
      OT.Levels(iLev).AuxSources.level = iLev;
      OT.Levels(iLev).AuxObs    .level = iLev;
      OT.Levels(iLev).AuxSources.set_part_size(nGroups);
      OT.Levels(iLev).AuxObs    .set_part_size(nGroups);
      for ( int g=1;g<=nGroups;g++){
	OT.Levels(iLev).AuxSources.set_part(g,Q,1,sm);
	OT.Levels(iLev).AuxObs    .set_part(g,Q,1,om);
	sm += Q;
	om += Q;
      }
    }
    ev.End();
    // Radiation
    // Equivalent sources distributions are evaluated for each group at each
    // level (if method is 'h', we always map actual (original) sources to equivalent
    // sources of each level; if method is 'h2' we recursively map equivalent sources of
    // level iLev to equivalent sources at level iLev-1 (parent level) )
    EventLog *eR = new EventLog("Radiation.");
    for ( int iLev = nLev; iLev >=3;iLev--){
      Level &ThisLev = OT(iLev) ;
      int nGroups = length(ThisLev.group) ;
      for (int iGroup = 1;iGroup <= nGroups;iGroup++){
	if ((iLev == nLev) || (tolower(OT.method) == 'h')){
	  // original sources to equivalent sources mapping
#ifdef FMM_3D
	  for (int jLev = iLev+1;jLev <= nLev;jLev++){
	    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
	  }
#endif // FMM_3D
	  if(dbg)fprintf(stdout,"Leaves\t\t\t Level(%d).V(%d) \t\t\t C(%d) \t\t Level(%d).S(%d)\n",
			 iLev,iGroup,iGroup,iLev,iGroup);
	  gemv_leaves ( OT(iLev).group(iGroup).V , x , OT.Levels(iLev).AuxSources,iGroup);
	}
	else{
	  // child to parent sources mapping
	  for ( int jj = 1;jj<= length(OT.Levels(iLev).group(iGroup).child); jj++){
	    int jGroup = OT(iLev).group(iGroup).child(jj) ;
	    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
	    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
	    if(dbg)fprintf(stdout,"Upward\t\t\t Level(%d).M2M(%d,%d) \t\t Level(%d).S(%d) \t\t Level(%d).S(%d)\n",
			   iLev,(int)i3d(1),(int)i3d(2)  ,iLev+1,jGroup    ,iLev,iGroup);
	    gemv_upward(OT(iLev).M2M(i3d(1),i3d(2)),OT.Levels(iLev+1).AuxSources,jGroup,OT.Levels(iLev).AuxSources,iGroup);
	  }
	}
      }
    }
    delete eR;
    // Translation
    // Translators transform the equivalent source density on source group to an
    // equivalent field density on the observation group: this operation is
    // performed, at each level, for all group pairs which are in far field
    // interaction list of each other
    EventLog *eT = new EventLog("Translation");
    for (int iLev = 3;iLev< nLev;iLev++){
      Level &ThisLev = OT.Levels(iLev) ;
      int nGroups = length(ThisLev.group) ;
      for (int iGroup = 1; iGroup <= nGroups;iGroup++){
	for ( int jj = 1; jj<=length(OT.Levels(iLev).group(iGroup).intList);jj++){
	  Time::TimeUnit tt = Time::getTime();
	  int jGroup = OT.Levels(iLev).group(iGroup).intList(jj) ;
	  double t1 = std::round(
				 ( OT.Levels(iLev).group(jGroup).groupcenter(1) -
				   OT.Levels(iLev).group(iGroup).groupcenter(1) ) /
				 OT.Levels(iLev).group(iGroup).cubelength ) ;
	  double t2 = std::round(
				 ( OT.Levels(iLev).group(jGroup).groupcenter(2) -
				   OT.Levels(iLev).group(iGroup).groupcenter(2) ) /
				 OT.Levels(iLev).group(iGroup).cubelength ) ;
	  int ix = t1;//t(1) ; 
	  int jy = t2;//t(2) ;
	  if(dbg)fprintf(stdout,"Translation \t\t Level(%d).T(%d,%d) \t\t Level(%d).S(%d) \t\t Level(%d).O(%d)\n",
			 iLev,4+ix,4+jy  ,iLev,iGroup    ,iLev,jGroup);
	  tt = Time::getTime() - tt;
	  stats.dur[1] += tt;

	  gemv_translation(OT.Levels(iLev).T(4+ix,4+jy) , OT.Levels(iLev).AuxSources,iGroup, OT.Levels(iLev).AuxObs,jGroup);
				

	}
      }
    }
    delete eT;
    // Receiving
    // The receiving step is reciprocal to the radiation step
    EventLog *eC = new EventLog ("Receiving");
    for (int iLev = 3;iLev <= nLev; iLev++){
      Level &ThisLev = OT(iLev) ;
      int nGroups = length(ThisLev.group) ;
      for ( int iGroup = 1;iGroup <= nGroups;iGroup ++){
	if ((iLev == nLev) || (tolower(OT.method) == 'h')){
	  //                Matrix &LeavesList = OT(iLev).group(iGroup).child ;
#ifdef FMM_3D
	  for ( int jLev = iLev+1;jLev<= nLev;jLev++){
	    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
	  }
#endif // FMM_3D
	  gemv_final(OT(iLev).group(iGroup).V,OT.Levels(iLev).AuxObs,iGroup,y);
	  if(dbg)fprintf(stdout,"Final\t\t\t\t Level(%d).V^T(%d) \t\t Level(%d).O(%d) \t\t Q(%d)\n",
			 iLev,iGroup  ,iLev,iGroup    ,iGroup);

	}
	else{
	  for ( int jj = 1; jj<=length(OT(iLev).group(iGroup).child);jj++){
	    int jGroup = OT(iLev).group(iGroup).child(jj) ;
	    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
	    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
	    gemv_down_obs(OT(iLev).M2M(i3d(1),i3d(2)) , OT.Levels(iLev).AuxObs,iGroup,OT.Levels(iLev+1).AuxObs,jGroup);
	    if(dbg)fprintf(stdout,"DownObs\t\t\t Level(%d).M2M^T(%d,%d) \t\t Level(%d).O(%d) \t\t Level(%d).O(%d)\n",
			   iLev,(int)i3d(1),(int)i3d(2), iLev,iGroup  ,iLev+1,jGroup);
	  }
	}
      }
    }

    delete eC;
  }

  /*====================================================================*/
  void gemv(SGMatrix &A,SGMatrix &x, SGMatrix &y){
    SGTaskGemv *t = new SGTaskGemv(A,x,y);
    t->type=7;
    t->set_name("NFx->y");
    submit(t);
  }
  /*====================================================================*/
  void gemvt(SGMatrix &A,SGMatrix &x, SGMatrix &y){
    SGTaskGemv *t = new SGTaskGemv(A,x,y);
    t->transA = true;
    t->type=8;
    t->set_name("NF^Tx->y");
    submit(t);
  }
  /*====================================================================*/
  void mv_near_field(Tree &OT,SGMatrix &C, SGMatrix &Q){
    bool dbg = !true;
    Matrix &c =* C.get_matrix();
    Matrix &q =* Q.get_matrix();
    double *qm = q.get_data_memory();
    double *cm = c.get_data_memory();
    int last = OT.num_levels();
    int ng=OT.Levels(last).group.size();
    C.set_part_size(ng);
    Q.set_part_size(ng);
    cout << "C part size: " << ng << endl;
    for(int g=1; g<=ng; g++)
      {
	int n = OT.Levels(last).group(g).child.size();
	int ofs=OT.Levels(last).group(g).child(1)-1;
	C.set_part(g,n,1,cm+ofs);
	Q.set_part(g,n,1,qm+ofs);
      }

    if (config.NF){
      int nfCount = OT.NearField.size();
      for (int i=0;i<nfCount;i++)
	{
	  NearFieldBlock *nf=OT.NearField[i];
	  gemv (*nf->sgMat,C.get_part(nf->j),Q.get_part(nf->i));
	  if(dbg)fprintf(stdout,"Nfld(%d,%d) \n",nf->j,nf->i);
	}
      if(dbg)fprintf(stdout,"Near Field Tasks submitted\n");
    }
  }

  /*=======================================================================================*/
  void compute_near_field(Tree &OT,Matrix &pts){
    /*
      for g in last level
      for n in nbr(g)
      for p1 in  pts(n.child)
      for p2 in pts(g.child)
      mat(i,j) = kernel(sqrt(sum((p1-p2)^2)))
    */
    int last = OT.num_levels();
    int ng=OT.Levels(last).group.size();
    for(int g=1; g<=ng;g++)
      {
        int nChildGroup=OT.Levels(last).group(g).child.size();
        int nNbr=OT.Levels(last).group(g).neighboursList.size();
        for (int n=1;n<=nNbr;n++)
	  {
            int Nbr = OT.Levels(last).group(g).neighboursList(n);
            int nChildNbr= OT.Levels(last).group(Nbr).child.size();
            Matrix &mat=*new Matrix ( nChildGroup, nChildNbr);
            int p = OT.Levels(last).group(Nbr).child(1);//starting index in the pts array
            NearFieldBlock *nfb  = new NearFieldBlock(g,Nbr,nChildGroup,nChildNbr,p,&mat);
            OT.NearField.push_back(nfb);
            for ( int p1=1;p1<=nChildGroup; p1++)
	      {
                for ( int p2=1;p2<=nChildNbr; p2++)
		  {
                    int i1=OT.Levels(last).group(g).child(p1);
                    int i2=OT.Levels(last).group(Nbr).child(p2);
                    if (i1==i2)
		      {
                        mat(p1,p2)=0;
                        continue;
		      }
                    double x1=pts(i1,1);
                    double y1=pts(i1,2);
                    double z1=pts(i1,3);

                    double x2=pts(i2,1);
                    double y2=pts(i2,2);
                    double z2=pts(i2,3);
                    double R = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
                    R =std::sqrt(R);
                    R =-std::log(R);
                    mat(p1,p2)=R;
		  }
	      }
	  }
      }
  }
  /*=======================================================================================*/
}//namespace FMM
