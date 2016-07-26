#include "fmm.hpp"

extern Config config;
Statistics stats;

Matrix &MatVec_old(Matrix &x ,  Tree & OT){
    Matrix &y= * new Matrix(size(x,1),size(x,2));
    int nLev = length(OT);
    int Q = OT.Q ;

    Tree &MatVecTemp=*new Tree;
    // INIT
    for (int iLev = 3;iLev<= nLev;iLev++){
        int nGroups = length(OT(iLev).group) ;
        MatVecTemp(iLev).AuxSources = zeros(Q,nGroups) ;    // at each level we store the auxiliary sources for each group
        MatVecTemp(iLev).AuxObs = zeros(Q,nGroups) ;        // at each level we store the auxiliary field density for each group
    }

    // Radiation
    // Equivalent sources distributions are evaluated for each group at each
    // level (if method is 'h', we always map actual (original) sources to equivalent
    // sources of each level; if method is 'h2' we recursively map equivalent sources of
    // level iLev to equivalent sources at level iLev-1 (parent level) )
    for ( int iLev = nLev; iLev >=3;iLev--){
        Level &ThisLev = OT(iLev) ;
        int nGroups = length(ThisLev.group) ;
        for (int iGroup = 1;iGroup <= nGroups;iGroup++){
            if ((iLev == nLev) || (tolower(OT.method) == 'h')){
                // original sources to equivalent sources mapping
                Matrix &LeavesList = OT.Levels(iLev).group(iGroup).child ;
                #ifdef FMM_3D
                for (int jLev = iLev+1;jLev <= nLev;jLev++){
                    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
                }
                #endif // FMM_3D
                SGMatrix &sg_x =  *new SGMatrix;
                sg_x = x.pick_rows(LeavesList);
                MatVecTemp(iLev).AuxSources.column(iGroup) = OT(iLev).group(iGroup).V * sg_x ;
            }
            else{
                // child to parent sources mapping
                for ( int jj = 1;jj<= length(OT.Levels(iLev).group(iGroup).child); jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
                    MatVecTemp(iLev).AuxSources.column(iGroup) = MatVecTemp(iLev).AuxSources.column(iGroup) +
                                OT(iLev).M2M(i3d(1),i3d(2)) * MatVecTemp(iLev+1).AuxSources.column(jGroup) ;
                }
            }
        }
    }

    // Translation
    // Translators transform the equivalent source density on source group to an
    // equivalent field density on the observation group: this operation is
    // performed, at each level, for all group pairs which are in far field
    // interaction list of each other
    for (int iLev = 3;iLev<= nLev;iLev++){
        Level &ThisLev = OT.Levels(iLev) ;
        int nGroups = length(ThisLev.group) ;
        for (int iGroup = 1; iGroup <= nGroups;iGroup++){
            for ( int jj = 1; jj<=length(OT.Levels(iLev).group(iGroup).intList);jj++){
                int jGroup = OT.Levels(iLev).group(iGroup).intList(jj) ;
                Matrix &t = round(
                                  ( OT.Levels(iLev).group(jGroup).groupcenter -
                                    OT.Levels(iLev).group(iGroup).groupcenter ) /
                                  OT.Levels(iLev).group(iGroup).cubelength ) ;
                int ix = t(1) ; int jy = t(2) ;
                MatVecTemp(iLev).AuxObs.column(jGroup) = MatVecTemp.Levels(iLev).AuxObs.column(jGroup) +
                    OT.Levels(iLev).T(4+ix,4+jy) * MatVecTemp.Levels(iLev).AuxSources.column(iGroup) ;
            }
        }
    }

    // Receiving
    // The receiving step is reciprocal to the radiation step
    for (int iLev = 3;iLev <= nLev; iLev++){
        Level &ThisLev = OT(iLev) ;
        int nGroups = length(ThisLev.group) ;
        for ( int iGroup = 1;iGroup <= nGroups;iGroup ++){
            if ((iLev == nLev) || (tolower(OT.method) == 'h')){
                Matrix &LeavesList = OT(iLev).group(iGroup).child ;
                #ifdef FMM_3D
                for ( int jLev = iLev+1;jLev<= nLev;jLev++){
                    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
                }
                #endif // FMM_3D
                SGMatrix &sg_x = *new SGMatrix;
                sg_x = x.pick_rows(LeavesList);
                MatVecTemp(iLev).AuxSources.column(iGroup) = OT(iLev).group(iGroup).V *  sg_x ;
                SGMatrix &sg_y = *new SGMatrix;
                sg_y = y.pick_rows(LeavesList);
                sg_y = sg_y + OT(iLev).group(iGroup).V.transpose() * MatVecTemp(iLev).AuxObs.column(iGroup) ;
            }
            else{
                for ( int jj = 1; jj<=length(OT(iLev).group(iGroup).child);jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
                    MatVecTemp(iLev+1).AuxObs.column(jGroup) = MatVecTemp(iLev+1).AuxObs.column(jGroup) +
                        OT(iLev).M2M(i3d(1),i3d(2)).transpose() * MatVecTemp(iLev).AuxObs.column(iGroup) ;
                }
            }
        }
    }

    return y;
}
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
vector<SGTaskGemv *> tlist;

void submit(SGTaskGemv *t){
    if (config.s and !config.t){
        t->run();
		return;
	}
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
}
void submit_all(){
	for ( auto t: tlist){
		sgEngine->submit(t);
	}
	fprintf(stdout,"submit all called.\n");
}
void gemv_leaves(SGMatrix &V, SGMatrix &x, SGMatrix &S, int group){
    SGTaskGemv *t= new SGTaskGemv(V,x,group,S,group);
	assert(V.get_matrix());
	assert(S.get_part(group).get_matrix());
	assert(x.get_part(group).get_matrix());
	t->type=1;
	t->set_name("Vx->s");
    submit(t);
}
void gemv_upward(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
	assert(M2M.get_matrix());
	assert(v.get_part(i1).get_matrix());
	assert(y.get_part(i2).get_matrix());
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
	t->type=2;
	t->set_name("Ms->s");
    submit(t);
}
void gemv_translation(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
	TL;
	assert(M2M.get_matrix());
	assert(v.get_part(i1).get_matrix());
	assert(y.get_part(i2).get_matrix());
	Time::TimeUnit tt = Time::getTime();
	TL;
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
	TL;
	tt = Time::getTime() - tt;
	TL;
	stats.dur[0] += tt;
	TL;
	t->type=3;
	t->set_name("Ms->o");
	TL;
    submit(t);

}
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
void MatVec(  Tree & OT,SGMatrix &x , SGMatrix &y){

    int nLev = length(OT);
    int Q = OT.Q ;
	bool dbg=true;

	EventLog ev("Init.");
	TL;

    // INIT
    for (int iLev = 3;iLev<= nLev;iLev++){
        int nGroups = length(OT(iLev).group) ;
        Matrix *S = new Matrix (Q,nGroups);
        Matrix *O = new Matrix (Q,nGroups);
        double *sm = S->get_data_memory();
        double *om = O->get_data_memory();
        OT.Levels(iLev).AuxSources.set_part_size(nGroups);
        OT.Levels(iLev).AuxObs    .set_part_size(nGroups);
	TL;
        for ( int g=1;g<=nGroups;g++){
            OT.Levels(iLev).AuxSources.set_part(g,Q,1,sm);
            OT.Levels(iLev).AuxObs    .set_part(g,Q,1,om);
            sm += Q;
            om += Q;
        }
    }
	TL;
	ev.End();
    // Radiation
    // Equivalent sources distributions are evaluated for each group at each
    // level (if method is 'h', we always map actual (original) sources to equivalent
    // sources of each level; if method is 'h2' we recursively map equivalent sources of
    // level iLev to equivalent sources at level iLev-1 (parent level) )
	EventLog *eR = new EventLog("Radiation.");
	TL;
    for ( int iLev = nLev; iLev >=3;iLev--){
        Level &ThisLev = OT(iLev) ;
        int nGroups = length(ThisLev.group) ;
	TL;
        for (int iGroup = 1;iGroup <= nGroups;iGroup++){
            if ((iLev == nLev) || (tolower(OT.method) == 'h')){
	TL;
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
	TL;
                // child to parent sources mapping
                for ( int jj = 1;jj<= length(OT.Levels(iLev).group(iGroup).child); jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
//                    OT.Levels(iLev).AuxSources.column(iGroup) = OT.Levels(iLev).AuxSources.column(iGroup) +
//                                OT(iLev).M2M(i3d(1),i3d(2)) * OT.Levels(iLev+1).AuxSources.column(jGroup) ;
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
	TL;
    for (int iLev = 3;iLev< nLev;iLev++){
        Level &ThisLev = OT.Levels(iLev) ;
        int nGroups = length(ThisLev.group) ;
	TL;
        for (int iGroup = 1; iGroup <= nGroups;iGroup++){
            for ( int jj = 1; jj<=length(OT.Levels(iLev).group(iGroup).intList);jj++){
				Time::TimeUnit tt = Time::getTime();
                int jGroup = OT.Levels(iLev).group(iGroup).intList(jj) ;
                double t1 = round(
                                  ( OT.Levels(iLev).group(jGroup).groupcenter(1) -
                                    OT.Levels(iLev).group(iGroup).groupcenter(1) ) /
                                  OT.Levels(iLev).group(iGroup).cubelength ) ;
                double t2 = round(
                                  ( OT.Levels(iLev).group(jGroup).groupcenter(2) -
                                    OT.Levels(iLev).group(iGroup).groupcenter(2) ) /
                                  OT.Levels(iLev).group(iGroup).cubelength ) ;
                /*Matrix &t = round(
                                  ( OT.Levels(iLev).group(jGroup).groupcenter -
                                    OT.Levels(iLev).group(iGroup).groupcenter ) /
                                  OT.Levels(iLev).group(iGroup).cubelength ) ;*/
                int ix = t1;//t(1) ; 
				int jy = t2;//t(2) ;
//                OT.Levels(iLev).AuxObs.column(jGroup) = OT.Levels(iLev).AuxObs.column(jGroup) +
//                    OT.Levels(iLev).T(4+ix,4+jy) * OT.Levels(iLev).AuxSources.column(iGroup) ;
                if(dbg)fprintf(stdout,"Translation \t\t Level(%d).T(%d,%d) \t\t Level(%d).S(%d) \t\t Level(%d).O(%d)\n",
                        iLev,4+ix,4+jy  ,iLev,iGroup    ,iLev,jGroup);
	TL;
				tt = Time::getTime() - tt;
				stats.dur[1] += tt;
	TL;

                gemv_translation(OT.Levels(iLev).T(4+ix,4+jy) , OT.Levels(iLev).AuxSources,iGroup, OT.Levels(iLev).AuxObs,jGroup);
	TL;
				

            }
        }
    }
	TL;
	delete eT;
    // Receiving
    // The receiving step is reciprocal to the radiation step
	EventLog *eC = new EventLog ("Receiving");
    for (int iLev = 3;iLev <= nLev; iLev++){
        Level &ThisLev = OT(iLev) ;
        int nGroups = length(ThisLev.group) ;
	TL;
        for ( int iGroup = 1;iGroup <= nGroups;iGroup ++){
            if ((iLev == nLev) || (tolower(OT.method) == 'h')){
//                Matrix &LeavesList = OT(iLev).group(iGroup).child ;
                #ifdef FMM_3D
                for ( int jLev = iLev+1;jLev<= nLev;jLev++){
                    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
                }
                #endif // FMM_3D
//                SGMatrix &sg_x = *new SGMatrix;
//                sg_x = x.pick_rows(LeavesList);
//                OT.Levels(iLev).AuxSources.column(iGroup) = OT(iLev).group(iGroup).V *  sg_x ;
//                gemv_downward( OT(iLev).group(iGroup).V ,  x ,OT.Levels(iLev).AuxSources, iGroup );
//                fprintf(stdout,"Downward\t\t\t Level(%d).V(%d) \t\t\t C(%d) \t\t\t Level(%d).S(%d)\n",
//                        iLev,iGroup  ,iGroup    ,iLev,iGroup);
//                SGMatrix &sg_y = *new SGMatrix;
//                sg_y = y.pick_rows(LeavesList);
//                sg_y = sg_y + OT(iLev).group(iGroup).V.transpose() * OT.Levels(iLev).AuxObs.column(iGroup) ;
                gemv_final(OT(iLev).group(iGroup).V,OT.Levels(iLev).AuxObs,iGroup,y);
                if(dbg)fprintf(stdout,"Final\t\t\t\t Level(%d).V^T(%d) \t\t Level(%d).O(%d) \t\t Q(%d)\n",
                        iLev,iGroup  ,iLev,iGroup    ,iGroup);

            }
            else{
	TL;
                for ( int jj = 1; jj<=length(OT(iLev).group(iGroup).child);jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
//                    OT.Levels(iLev+1).AuxObs.column(jGroup) = OT.Levels(iLev+1).AuxObs.column(jGroup) +
//                        OT(iLev).M2M(i3d(1),i3d(2)).transpose() * OT.Levels(iLev).AuxObs.column(iGroup) ;
                    gemv_down_obs(OT(iLev).M2M(i3d(1),i3d(2)) , OT.Levels(iLev).AuxObs,iGroup,OT.Levels(iLev+1).AuxObs,jGroup);
                    if(dbg)fprintf(stdout,"DownObs\t\t\t Level(%d).M2M^T(%d,%d) \t\t Level(%d).O(%d) \t\t Level(%d).O(%d)\n",
                            iLev,(int)i3d(1),(int)i3d(2), iLev,iGroup  ,iLev+1,jGroup);
                }
            }
        }
    }
	TL;

	delete eC;
}

/*====================================================================*/
void gemv(SGMatrix &A,SGMatrix &x, SGMatrix &y){
	TL;
    SGTaskGemv *t = new SGTaskGemv(A,x,y);
	t->type=7;
	t->set_name("NFx->y");
	TL;
    submit(t);
	TL;
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
	bool dbg = true;
    Matrix &c =* C.get_matrix();
    Matrix &q =* Q.get_matrix();
    double *qm = q.get_data_memory();
    double *cm = c.get_data_memory();
    int last = OT.num_levels();
    int ng=OT.Levels(last).group.size();
    C.set_part_size(ng);
    Q.set_part_size(ng);
	cout << "C part size: " << ng << endl;
	TL;
    for(int g=1; g<=ng; g++)
    {
        int n = OT.Levels(last).group(g).child.size();
        int ofs=OT.Levels(last).group(g).child(1)-1;
	TL;
        C.set_part(g,n,1,cm+ofs);
	TL;
        Q.set_part(g,n,1,qm+ofs);
    }
	TL;

	if (config.NF){
	  int nfCount = OT.NearField.size();
	  TL;
	  for (int i=0;i<nfCount;i++)
	    {
	      NearFieldBlock *nf=OT.NearField[i];
	      TL;
	      gemv (*nf->sgMat,C.get_part(nf->j),Q.get_part(nf->i));
	      if(dbg)fprintf(stdout,"Nfld(%d,%d) \n",nf->j,nf->i);
	    }
	  TL;
	  if(dbg)fprintf(stdout,"Near Field Tasks submitted\n");
	  TL;
	}
}
