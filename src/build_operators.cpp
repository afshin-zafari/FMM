#include "fmm.hpp"
/*=======================================================================================*/
void compute_translator(Tree &OT,KernelFcn kernel){
    double theta;
    int cnt;
    Matrix TauSurf(OT.Q+1,3);
    for ( theta = 0.0,cnt=1; theta < 2* M_PI && cnt <= OT.Q+1; theta +=2*M_PI /(OT.Q+1),cnt++){
        TauSurf(cnt,1)=0.5 * sqrt(2)*cos(theta);
        TauSurf(cnt,2)=0.5 * sqrt(2)*sin(theta);
        TauSurf(cnt,3)=0   * theta;//FMM_3D
    }
    int nLev = OT.num_levels();
    for ( int iLev = 3; iLev < nLev; iLev++){
        Level & ThisLevel = OT.Levels(iLev);
        Matrix &SourcePoints=ThisLevel.group(1).cubelength * TauSurf;
        for ( int ix=-3;ix<=3;ix++){
            for(int jy =-3;jy<=3;jy++){
#               ifdef FMM_3D
                for(int kz=-3;kz<=3;kz++){
#               endif // FMM_3D
                Matrix V(1,3);
                V(1,1)=ix * ThisLevel.group(1).cubelength;
                V(1,2)=jy * ThisLevel.group(1).cubelength;
#               ifdef FMM_3D
                V(1,3)=kz * ThisLevel.group(1).cubelength;
                #else
                V(1,3)=0;
#               endif
                Matrix &TestPoints=SourcePoints + repmat(V,OT.Q,1);
                Matrix &xtest = repmat(TestPoints.column(1),1,SourcePoints.rows()) ;
                Matrix &ytest = repmat(TestPoints.column(2),1,SourcePoints.rows()) ;
                Matrix &ztest = repmat(TestPoints.column(3),1,SourcePoints.rows()) ;

                Matrix &xsource = repmat(SourcePoints.column(1),1,TestPoints.rows()) ;
                Matrix &ysource = repmat(SourcePoints.column(2),1,TestPoints.rows()) ;
                Matrix &zsource = repmat(SourcePoints.column(3),1,TestPoints.rows()) ;
                Matrix &X = xtest-xsource.transpose();
                Matrix &Y = ytest-ysource.transpose();
                Matrix &Z = ztest-zsource.transpose();
                Matrix &R = sqrt( (X^2) + (Y^2) + (Z^2)) ;
                Matrix &T = kernel(R) ;
                OT.Levels(iLev).T(4+ix,4+jy) = T ;
#           ifdef FMM_3D
                }
#           endif // FMM_3D
            }
        }
    }
}
/*=======================================================================================*/
void compute_radiation( Tree &OT, Matrix &pts, char method, KernelFcn kernel){
    double theta;
    int cnt;
    Matrix TauSurf(OT.Q+1,3);
    for ( theta = 0.0,cnt=1; theta < 2* M_PI && cnt <= OT.Q+1; theta +=2*M_PI /(OT.Q+1),cnt++){
        TauSurf(cnt,1)=0.5 * sqrt(2)*cos(theta);
        TauSurf(cnt,2)=0.5 * sqrt(2)*sin(theta);
        TauSurf(cnt,3)=0   * theta;//FMM_3D
    }
    Matrix SigmaSurf(2*OT.Q+1,3);
    for ( theta = 0.0,cnt=1; theta < 2* M_PI && cnt <= 2*OT.Q+1; theta +=M_PI /(2*OT.Q+1),cnt++){
        SigmaSurf(cnt,1)=1.5 * cos(theta);
        SigmaSurf(cnt,2)=1.5 * sin(theta);
        SigmaSurf(cnt,3)=0   * theta;//FMM_3D
    }
    int nLev = OT.num_levels();
    for ( int iLev = 3; iLev < nLev; iLev++){
        Level &ThisLevel = OT.Levels(iLev);
        int nGroups = length(ThisLevel.group);
        Matrix &SourceFieldPoints = TauSurf*ThisLevel.group(1).cubelength ;
        Matrix &TestFieldPoints = SigmaSurf*ThisLevel.group(1).cubelength ;
        Matrix &xtest = repmat(TestFieldPoints.column(1),1,size(SourceFieldPoints,1)) ;
        Matrix &ytest = repmat(TestFieldPoints.column(2),1,size(SourceFieldPoints,1)) ;
        Matrix &ztest = repmat(TestFieldPoints.column(3),1,size(SourceFieldPoints,1)) ;

        Matrix &xsource = repmat(SourceFieldPoints.column(1),1,size(TestFieldPoints,1)) ;
        Matrix &ysource = repmat(SourceFieldPoints.column(2),1,size(TestFieldPoints,1)) ;
        Matrix &zsource = repmat(SourceFieldPoints.column(3),1,size(TestFieldPoints,1)) ;
        Matrix &X = xtest-xsource.transpose() ;
        Matrix &Y = ytest-ysource.transpose() ;
        Matrix &Z = ztest-zsource.transpose() ;
        Matrix &R = sqrt( (X^2) + (Y^2) + (Z^2)) ;
        Matrix &B = kernel(R) ;
        Matrix U,S,V;
        SVD(B,U,S,V);
        Matrix &Sd = diag(S) ;
        int Ns = numel(Sd> (svd_threshold*Sd(1))) ;
        Sd = 1.0/Sd.row_slice(1,Ns) ;
        Matrix &Binv = V.col_slice(1,Ns) * diag(Sd) * U.col_slice(1,Ns).transpose() ;

        if ((iLev == nLev) || (tolower(method) == 'h')){
            // in the H2 method bases are nested, i.e. expressed
            // recursively in terms of the bases at leaf level (nLev); only
            // bases at leaf level need be computed and stored
            for (int iGroup = 1;iGroup <= nGroups;iGroup++){
                GroupType &TestGroup = ThisLevel.group(iGroup) ;
                TestFieldPoints = TestFieldPoints + repmat(TestGroup.groupcenter,2*OT.Q,1) ;
                Matrix &LeavesList = TestGroup.child ;
                #ifdef FMM_H_METHOD
                for (int jLev = iLev+1;j<=nLev;jLev ++){
                    LeavesList = vertcat(OT(jLev).group(LeavesList).child) ;
                }
                #endif
                Matrix &SourcePoints = pts.pick_rows(LeavesList) ;

                Matrix &xtest = repmat(TestFieldPoints.column(1),1,size(SourcePoints,1)) ;
                Matrix &ytest = repmat(TestFieldPoints.column(2),1,size(SourcePoints,1)) ;
                Matrix &ztest = repmat(TestFieldPoints.column(3),1,size(SourcePoints,1)) ;

                Matrix &xsource = repmat(SourcePoints.column(1),1,size(TestFieldPoints,1)) ;
                Matrix &ysource = repmat(SourcePoints.column(2),1,size(TestFieldPoints,1)) ;
                Matrix &zsource = repmat(SourcePoints.column(3),1,size(TestFieldPoints,1)) ;
                Matrix &X = xtest-xsource.transpose() ;
                Matrix &Y = ytest-ysource.transpose() ;
                Matrix &Z = ztest-zsource.transpose() ;
                Matrix &R = sqrt( (X^2) + (Y^2) + (Z^2)) ;
                Matrix &W = kernel(R) ;

                TestFieldPoints = TestFieldPoints - repmat(TestGroup.groupcenter,2*OT.Q,1) ;

                // V = B^(-1) * W is the radiation pattern mapping actual
                // sources to equivalent sources, by matching the field on
                // the testing surface
                OT.Levels(iLev).group(iGroup).V = Binv * W ;
            }
        }
        else{
            // Compute transfer matrix to ascend/descend the tree
            // Given a parent group, there are at most 4 children (8 for 3D
            // problems); the relative location parent/son does not depend
            // on the specific group, there are only 4 (8) transfer
            // operators at each level.
            for (int ix = 1;ix<=2;ix++){
                for(int jy = 1;jy<= 2;jy++){
                    #ifdef FMM_3D
                    for ( int kz = 1; kz<=2 ; kz++){     // 3D problems
                    #endif // FMM_3D
                        // remove translations along z axis (2D problem => quad-tree)
                        Matrix &Temp=*new Matrix(1,3);
                        Temp(1,1)=ix -1.5;
                        Temp(1,2)=jy -1.5;
                        Temp(1,3)=0;
                        Matrix &Shift = 0.5*ThisLevel.group(1).cubelength * Temp;//([ix jy 1.5] - 1.5) ;     // vector joining barycenters of parent and son groups
                        Matrix &SourcePoints = 0.5*TauSurf*ThisLevel.group(1).cubelength + repmat(Shift,OT.Q,1) ;
                        Matrix &xtest = repmat(TestFieldPoints.column(1),1,size(SourcePoints,1)) ;
                        Matrix &ytest = repmat(TestFieldPoints.column(2),1,size(SourcePoints,1)) ;
                        Matrix &ztest = repmat(TestFieldPoints.column(3),1,size(SourcePoints,1)) ;

                        Matrix &xsource = repmat(SourcePoints.column(1),1,size(TestFieldPoints,1)) ;
                        Matrix &ysource = repmat(SourcePoints.column(2),1,size(TestFieldPoints,1)) ;
                        Matrix &zsource = repmat(SourcePoints.column(3),1,size(TestFieldPoints,1)) ;
                        Matrix &X = xtest-xsource.transpose() ;
                        Matrix &Y = ytest-ysource.transpose() ;
                        Matrix &Z = ztest-zsource.transpose() ;
                        Matrix &R = sqrt( (X^2) + (Y^2) + (Z^2) );
                        Matrix &W = kernel(R) ;

// V = B^(-1) * W is the transfer operator mapping equivalent sources at
// level iLev to equivalent sources at level iLev-1 (parent level), by
// matching the field on the testing surface of parent level
//                         A.Level(iLev).M2M{ix,jy,kz} = Binv * W ;  % 3D M2M operator
                        OT.Levels(iLev).M2M(ix,jy) = Binv * W ;
                        #ifdef FMM_3D
                    }
                    #endif // FMM_3D
                }
            }

        }


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
