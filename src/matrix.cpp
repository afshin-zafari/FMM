#include <assert.h>
#include <vector>
#include <iomanip>
#include <fstream>

#include "matrix.hpp"

using namespace std;
Matrix &repmat(Matrix &M, int r, int c)
{
    assert(r >0 && c >0);
    if ( M.is_locked()){
        return M.do_lock_op<Matrix &>("repmat");
    }
    Matrix &X=* new Matrix (r*M.rows(),c*M.cols());

    for ( int i=0; i<X.rows(); i++)
    {
        for ( int j=0; j<X.cols(); j++)
        {
            X(i+1,j+1) = M(i%M.rows()+1,j%M.cols()+1);
        }
    }
    return X;
}
Matrix &sqrt(Matrix &M)
{
    return M.sqrt();
}
Matrix &log_kernel(Matrix &M)
{
    return M.log();
}
int size(Matrix &M, int rc)
{
    assert(rc >=1 && rc <=3);
    if ( rc ==1)
        return M.rows();
    return M.cols();
}
Matrix &diag (Matrix &M)
{
    return M.diag();
}
Matrix &round(Matrix &M)
{
    return M.round();
}
int numel(Matrix &M)
{
    return M.numel();
}


void SVD(Matrix &A, Matrix &U, Matrix &S, Matrix &V)
{
    int m=A.rows(),n=A.cols();
    printf("-------------------\n");
    A.print();
    Matrix B(m<n?m:n,1);
    ElementType *work = B.get_data_memory();

    U.rebuild(m,n);
    S.rebuild(n,1);
    V.rebuild(n,n);

    int v_ldu = V.rows();
	int ret;

    /*
    LAPACKE_dgesvd( matrix_layout ,char jobu, char jobvt,
                   lapack_int m, lapack_int n, double* a,
                   lapack_int lda, double* s, double* u, lapack_int ldu,
                   double* vt, lapack_int ldvt, double* superb );
   */
	/*
    lapack_int ret=LAPACKE_dgesvd( LAPACK_COL_MAJOR,'A','A',
                           m,n, A.get_data_memory(),m, S.get_data_memory(), U.get_data_memory(), U.rows(),
                           V.get_data_memory(), v_ldu, work);
						   */
   if ( ret >0 )
        fprintf(stderr,"SVD failed: %d\n",ret);
    printf("-------------------\n");
   S.print();
    printf("-------------------\n");
   U.print();
    printf("-------------------\n");
   V.print();
    printf("-------------------\n");


}
Matrix &zeros(int r, int c)
{
    assert(r>0 && c>0);
    Matrix &X = * new Matrix (r,c);
    return X;
}
Matrix &Matrix::max(int dim){
    int m = (dim ==1)?1:M;
    int n = (dim ==2)?1:N;
    if(locked){
        return do_lock_op<Matrix &>("max (Mat)");
    }
    Matrix &X=*new Matrix(m,n);
    if ( m==1){
        for ( int j=0;j<N;j++){
            ElementType mx=-1e-100;
            for ( int i=0;i<=M;i++){
                if ( data[j*M+i]>mx)
                    mx = data[j*M+i];
            }
            X(1,j+1)=mx;
        }
    }
    if ( n==1){
        for ( int i=0;i<=M;i++){
            ElementType mx=-1e-100;
            for ( int j=0;j<N;j++){
                if ( data[j*M+i]>mx)
                    mx = data[j*M+i];
            }
            X(i+1,1)=mx;
        }
    }
    return X;
}
Matrix &Matrix::min(int dim){
    int m = (dim ==1)?1:M;
    int n = (dim ==2)?1:N;
    if ( locked)
        return do_lock_op<Matrix &>("min(Mat)");
    Matrix &X=*new Matrix(m,n);
    if ( m==1){
        for ( int j=0;j<N;j++){
            ElementType mx=+1e+100;
            for ( int i=0;i<=M;i++){
                if ( data[j*M+i]<mx)
                    mx = data[j*M+i];
            }
            X(1,j+1)=mx;
        }
    }
    if ( n==1){
        for ( int i=0;i<=M;i++){
            ElementType mx=-1e-100;
            for ( int j=0;j<N;j++){
                if ( data[j*M+i]<mx)
                    mx = data[j*M+i];
            }
            X(i+1,1)=mx;
        }
    }
    return X;
}
Matrix &Matrix::floor(){
    if(locked)
        return do_lock_op<Matrix &>("floor(Mat)");
    Matrix &X=*new Matrix(M,N);
    for(int i=0;i<M;i++)
        for(int j=0;j<N;j++)
            X(i+1,j+1) = std::floor(data[j*M+i]);
    return X;
}
Matrix &Matrix::get_column(int c){
    if(locked)
        return do_lock_op<Matrix &>("get column (Mat)");
    c--;
    assert(c>=0 && c<N);
    Matrix &X=*new Matrix(M,1);
    for ( int i=0;i<M;i++)
        X(i+1,1)=data[c*M+i];
    return X;
}
Matrix &Matrix::get_row(int r){
    if(locked)
        return do_lock_op<Matrix &>("get row(Mat)");
    r--;
    assert(r>=0 && r<M);
    Matrix &X=*new Matrix(1,N);
    for ( int j=0;j<N;j++)
        X(1,j+1)=data[j*M+r];
    return X;
}
Matrix &hcat(Matrix &A, Matrix &B){
    assert(A.rows() == B.rows());
    int M = A.rows();
    int N = A.cols()+B.cols();
    Matrix &X=*new Matrix(M,N);
    for ( int i=1;i<=A.rows();i++){
        for ( int j=1;j<=A.cols();j++){
            X(i,j)=A(i,j);
        }
    }
    for ( int i=1;i<=A.rows();i++){
        for ( int j=1;j<=B.cols();j++){
            X(i,A.cols()+j)=B(i,j);
        }
    }
    return X;
}
Matrix &Matrix::set_column(int c, Matrix &rhs,bool build){
    if(locked)
        return do_lock_op<Matrix &>("set column(Mat)");
    if (!build){
        c--;
        assert(c<N && c>=0);
        assert((rhs.rows() * rhs.cols()) == M);
    }
    else {
        M = rhs.rows()*rhs.cols();
        N=1;
        c=0;
        if( data)
            delete [] data;
        data = new ElementType [M];
    }
    for ( int i=0;i<M;i++){
        data[c*M+i] = rhs(i+1);
    }
    return *this;
}
Matrix &Matrix::append(ElementType rhs){
    if(locked)
        return do_lock_op<Matrix &>("append to Mat");
    assert(M==1 || N==1);
    if ( M==1){
        Matrix &X=*new Matrix(1,N+1);
        for(int j=0;j<N;j++)
            X(1,j+1) = data[j*M];
        X(1,N+1)=rhs;
        return X;
    }
    if ( N==1){
        Matrix &X=*new Matrix(M+1,1);
        for(int i=0;i<M;i++)
            X(i+1,1) = data[i];
        X(M+1,1)=rhs;
        return X;
    }

    return *this;
}
Matrix &Matrix::diff(){
    if(locked)
        return do_lock_op<Matrix &>("diff(Mat)");
    assert(M==1 || N==1);
    if ( N==1){
        Matrix &X=*new Matrix(M-1,1);
        for ( int i=0;i<M-1;i++)
            X(i+1,1) = data[i]-data[i+1];
        return X;
    }
    if ( M==1){
        Matrix &X=*new Matrix(1,N-1);
        for ( int j=0;j<N-1;j++)
            X(1,j+1) = data[j]-data[j+1];
        return X;
    }
    assert(M!=1 && N!=1);
    return *this;
}
Matrix &Matrix::find(){
    if(locked)
        return do_lock_op<Matrix &>("find(Mat)");
    std::vector<int> nz;
    for(int j=0;j<N;j++)
        for(int i=0;i<M;i++)
            if (data[j*M+i] !=0.0)
                nz.push_back((j+0)*M+(i+1));//(i+1,j+1)=bool(data[j*M+i]);
    Matrix &X=*new Matrix(nz.size(),1);
    for ( unsigned int i=1;i<=nz.size(); i++)
        X(i)=nz[i-1];
    return X;
}
void Matrix::print(const char *s){
    cout << s << endl;
    for ( int i=0; i< M; i++){
        for (int j=0;j<N;j++){
            //printf("%2.2lf,",data[i+j*N]);
            cout << setprecision (8) << setw(10) << data[i+j*M] << ',' ;
        }
        //printf("\n");
        cout << endl;
    }
}
void Matrix::export_data(const char *fn){
    std::ofstream f;
    f.open(fn);
    for ( int i=0; i< M; i++){
        for (int j=0;j<N;j++){
            f << setprecision (8) << setw(10) << data[i+j*M] << ',' ;
        }
        f<< endl;
    }
    f.close();

}
Matrix &max(Matrix &M, int dim){
    return M.max(dim);
}
Matrix &min(Matrix &M, int dim ){
    return M.min(dim);
}
Matrix &floor(Matrix &M){
    return M.floor();
}
Matrix &diff(Matrix &M){
    return M.diff();
}
Matrix &find(Matrix &M){
    return M.find();
}

int length(Matrix &M){return M.rows()*M.cols();}
