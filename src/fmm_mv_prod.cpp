#include "fmm.hpp"
#include "sg/superglue.hpp"
#include "sg/option/instr_trace.hpp"

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
/*
class SGTaskGemv : public Task<Options, 3> {
private:
    SGMatrix *A,*v,*y;
public:
    bool transA;
    enum{
		Read=ReadWriteAdd::read,
		Write=ReadWriteAdd::write,
		Add=ReadWriteAdd::add};
    enum{COL_MAJOR,ROW_MAJOR};
	std::string get_name() { return "A"; }
    SGTaskGemv(SGMatrix &A_, SGMatrix &v_, int i1,SGMatrix &Y_, int i2){
        A = &A_;
        v = &v_.get_part(i1);
        y = &Y_.get_part(i2);
        register_args();
    }
    SGTaskGemv(SGMatrix &A_, SGMatrix &v_, SGMatrix &Y_){
        A = &A_;
        v = &v_;
        y = &Y_;
        register_args();
    }
    void register_args(){
        SGHandle &hA = A->get_handle();
        SGHandle &hv = v->get_handle();
        SGHandle &hy = y->get_handle();
        register_access(ReadWriteAdd::read, hA);
        register_access(ReadWriteAdd::read, hv);
        if ( config.w )
            register_access(ReadWriteAdd::write, hy);
        else
            register_access(ReadWriteAdd::add, hy);
        transA = false;

    }
    //void register_access(int axs, SGHandle &h){}
    void run(){
        const int M = A->get_matrix()->rows();
        const int N = A->get_matrix()->cols();
        const int lda = M;
        const double *Mat = A->get_matrix()->get_data_memory();
        const double *X   = v->get_matrix()->get_data_memory();
        double *Y   = y->get_matrix()->get_data_memory();
        assert(Mat);
        assert(X);
        assert(Y);
        int mx= v->get_matrix()->rows();
        int nx= v->get_matrix()->cols();
        int my= y->get_matrix()->rows();
        int ny= y->get_matrix()->cols();
        //fprintf(stdout,"A %dx%d X %dx%d Y%dx%d\n",M,N,mx,nx,my,ny);
        assert(nx*ny==1);
        if ( !transA){
            assert(M==my );
            assert(N==mx );
        }else{
            assert(N==my );
            assert(M==mx );
        }
        cblas_dgemv(COL_MAJOR,transA,M, N, 1.0, Mat, lda, X, 1, 1.0, Y, 1);

    }
};
*/
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
		sgEngine->submit(t);
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
    submit(t);
}
void gemv_upward(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
    submit(t);
}
void gemv_translation(SGMatrix &M2M, SGMatrix &v, int i1,SGMatrix &y,int i2){
    SGTaskGemv *t= new SGTaskGemv(M2M,v,i1,y,i2);
    submit(t);

}
void gemv_final(SGMatrix &A, SGMatrix &v, int i, SGMatrix &y){
    SGTaskGemv *t= new SGTaskGemv(A,v,i,y,i);
    t->transA = true;
    submit(t);    
}
void gemv_downward(SGMatrix &A, SGMatrix &x, SGMatrix &y,int i ){
    y.get_part(i).get_matrix()->print();
    SGTaskGemv *t= new SGTaskGemv(A,x,i,y,i);
    submit(t);
}
void gemv_down_obs(SGMatrix &A, SGMatrix &v, int i1, SGMatrix &y, int i2){
    SGTaskGemv *t= new SGTaskGemv(A,v,i1,y,i2);
    t->transA = true;
    submit(t);
}
void MatVec(  Tree & OT,SGMatrix &x , SGMatrix &y){

    int nLev = length(OT);
    int Q = OT.Q ;

    // INIT
    for (int iLev = 3;iLev<= nLev;iLev++){
        int nGroups = length(OT(iLev).group) ;
        Matrix *S = new Matrix (Q,nGroups);
        Matrix *O = new Matrix (Q,nGroups);
        double *sm = S->get_data_memory();
        double *om = O->get_data_memory();
        OT.Levels(iLev).AuxSources.set_part_size(nGroups);
        OT.Levels(iLev).AuxObs    .set_part_size(nGroups);
        for ( int g=1;g<=nGroups;g++){
            OT.Levels(iLev).AuxSources.set_part(g,Q,1,sm);
            OT.Levels(iLev).AuxObs    .set_part(g,Q,1,om);
            sm += Q;
            om += Q;
        }
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
                #ifdef FMM_3D
                for (int jLev = iLev+1;jLev <= nLev;jLev++){
                    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
                }
                #endif // FMM_3D
                gemv_leaves ( OT(iLev).group(iGroup).V , x , OT.Levels(iLev).AuxSources,iGroup);
                if(0)fprintf(stdout,"Leaves\t\t\t Level(%d).V(%d) \t\t\t C(%d) \t\t Level(%d).S(%d)\n",
                        iLev,iGroup,iGroup,iLev,iGroup);
            }
            else{
                // child to parent sources mapping
                for ( int jj = 1;jj<= length(OT.Levels(iLev).group(iGroup).child); jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
//                    OT.Levels(iLev).AuxSources.column(iGroup) = OT.Levels(iLev).AuxSources.column(iGroup) +
//                                OT(iLev).M2M(i3d(1),i3d(2)) * OT.Levels(iLev+1).AuxSources.column(jGroup) ;
                    if(0)fprintf(stdout,"Upward\t\t\t Level(%d).M2M(%d,%d) \t\t Level(%d).S(%d) \t\t Level(%d).S(%d)\n",
                            iLev,(int)i3d(1),(int)i3d(2)  ,iLev+1,jGroup    ,iLev,iGroup);
                    gemv_upward(OT(iLev).M2M(i3d(1),i3d(2)),OT.Levels(iLev+1).AuxSources,jGroup,OT.Levels(iLev).AuxSources,iGroup);
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
//                OT.Levels(iLev).AuxObs.column(jGroup) = OT.Levels(iLev).AuxObs.column(jGroup) +
//                    OT.Levels(iLev).T(4+ix,4+jy) * OT.Levels(iLev).AuxSources.column(iGroup) ;
                if(0)fprintf(stdout,"Translation \t\t Level(%d).T(%d,%d) \t\t Level(%d).S(%d) \t\t Level(%d).O(%d)\n",
                        iLev,4+ix,4+jy  ,iLev,iGroup    ,iLev,jGroup);
                gemv_translation(OT.Levels(iLev).T(4+ix,4+jy) , OT.Levels(iLev).AuxSources,iGroup, OT.Levels(iLev).AuxObs,jGroup);

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
                if(0)fprintf(stdout,"Final\t\t\t\t Level(%d).V^T(%d) \t\t Level(%d).O(%d) \t\t Q(%d)\n",
                        iLev,iGroup  ,iLev,iGroup    ,iGroup);

            }
            else{
                for ( int jj = 1; jj<=length(OT(iLev).group(iGroup).child);jj++){
                    int jGroup = OT(iLev).group(iGroup).child(jj) ;
                    TempDoubleVector &t = OT(iLev+1).group(jGroup).groupcenter - OT(iLev).group(iGroup).groupcenter ;
                    Matrix &i3d = round( ( sign(t) + 3 ) / 2 ) ;
//                    OT.Levels(iLev+1).AuxObs.column(jGroup) = OT.Levels(iLev+1).AuxObs.column(jGroup) +
//                        OT(iLev).M2M(i3d(1),i3d(2)).transpose() * OT.Levels(iLev).AuxObs.column(iGroup) ;
                    gemv_down_obs(OT(iLev).M2M(i3d(1),i3d(2)) , OT.Levels(iLev).AuxObs,iGroup,OT.Levels(iLev+1).AuxObs,jGroup);
                    if(0)fprintf(stdout,"DownObs\t\t\t Level(%d).M2M^T(%d,%d) \t\t Level(%d).O(%d) \t\t Level(%d).O(%d)\n",
                            iLev,(int)i3d(1),(int)i3d(2), iLev,iGroup  ,iLev+1,jGroup);
                }
            }
        }
    }


}

/*====================================================================*/
void gemv(SGMatrix &A,SGMatrix &x, SGMatrix &y){
    SGTaskGemv *t = new SGTaskGemv(A,x,y);
    submit(t);
//    y.get_matrix()->print();
}
/*====================================================================*/
void gemvt(SGMatrix &A,SGMatrix &x, SGMatrix &y){
    SGTaskGemv *t = new SGTaskGemv(A,x,y);
    t->transA = true;
    submit(t);
}
/*====================================================================*/
void mv_near_field(Tree &OT,SGMatrix &C, SGMatrix &Q){
    Matrix &c =* C.get_matrix();
    Matrix &q =* Q.get_matrix();
    double *qm = q.get_data_memory();
    double *cm = c.get_data_memory();
    int last = OT.num_levels();
    int ng=OT.Levels(last).group.size();
    C.set_part_size(ng);
    Q.set_part_size(ng);
    for(int g=1; g<=ng;g++)
    {
        int n = OT.Levels(last).group(g).child.size();
        int ofs=OT.Levels(last).group(g).child(1)-1;
        C.set_part(g,n,1,cm+ofs);
        Q.set_part(g,n,1,qm+ofs);
    }


    int nfCount = OT.NearField.size();
    bool *worked = new bool[ng];
    for(int i=0;i<ng;i++)
        worked[i]=false;
    for (int i=0;i<nfCount;i++)
    {
        NearFieldBlock *nf=OT.NearField[i];
//        if ( worked[nf->i] and worked[nf->j] )
//            continue;
        worked[nf->i] = worked[nf->j]=true;
//        fprintf(stdout,"Nfld(%d,%d) \n",nf->i,nf->j);
        gemv (*nf->sgMat,C.get_part(nf->j),Q.get_part(nf->i));
//        fprintf(stdout,"Nfld(%d,%d)^T \n",nf->i,nf->j);
//        gemvt(*nf->sgMat,C.get_part(nf->i),Q.get_part(nf->j));
    }
    //fprintf(stdout,"Near Field Tasks submitted\n");
}
