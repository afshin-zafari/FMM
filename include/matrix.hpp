#ifndef MATRIX_HPP_INCLUDED
#define MATRIX_HPP_INCLUDED

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

//#define FOR_ALL_ij         for(int i=0;i<M;i++)   for(int j=0;j<N;j++)
#define ROW_RANGE(i)  ((i>=0)&&(i<M))
#define COL_RANGE(j)   ((j>=0)&& (j<N))
#include <vector>
typedef double ElementType;
template <typename T> int sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}
class Matrix
{
private:
    ElementType *data;
    int M,N;
    bool zero_indexing;
    bool locked,temp;
    std::vector<Matrix*> pcols;
    Matrix *parent;
public:
    void lock(){locked = true;if(parent)parent->lock();}
    void unlock(){locked =false;}
    bool is_locked(){return locked;}
    bool is_temp(){return temp;}
    void set_temporary(bool f= true){temp = f;}
    void clear(){
        for ( int i = 0; i<M*N;i++){
            data[i]=0;
        }
    }
    template <typename T>
    Matrix & do_lock_op(const char *s){
        std::cout << "locked op: " << s << std::endl;
        return *this;
    }
    void do_lock_op2(const char *s)const {
        std::cout << "locked op: " << s << std::endl;
    }
    void set_parent(Matrix *p){parent = p;}
    Matrix (int m, int n , double * mem ){
        M = m;
        N = n;
        data = mem;
    }
    Matrix (int m, int n, double val){
        M = m;
        N = n;
        assert(N>0 && M>0);
        data = new ElementType[M*N];
        zero_indexing=false;
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                data[j*M+i]=val;
        locked=false;
        temp = false;
        parent = NULL;
    }
    Matrix( int M_ , int N_ , int p =0):M(M_),N(N_)
    {
        if ( M==0 || N==0)
            M=0;

        assert(N>0 && M>0);
        data = new ElementType[M*N];
        zero_indexing=false;
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                data[j*M+i]=0.0;
        locked=false;
        temp = false;
        parent = NULL;
    }
    ~Matrix ()
    {
        if ( data)
            delete [] data;
    }
    ElementType& operator () (int i, int j=-1 )
    {
        if ( j==-1){
            assert( N ==1 || M== 1);
            if ( i==-1 ){
                if(N==1)i=M,j=1;
                if(M==1)i=1,j=M;
            }
            else{
                if ( N==1){j=1;}
                if(M==1){j=i;i=1;}
            }
        }
        if(!zero_indexing)
        {
            i--;
            j--;
        }
        assert(ROW_RANGE(i));
        assert(COL_RANGE(j));
        assert(data);
        if ( locked){
            do_lock_op2("read access i,j");
            return *data;
        }
        return data[j*M+i];
    }
    ElementType operator () (int i, int j ) const
    {
        if(!zero_indexing)
        {
            i--;
            j--;
        }
        assert(ROW_RANGE(i));
        assert(COL_RANGE(j));
        assert(data);
        if ( locked){
            do_lock_op2("write access i,j");
            return 0;
        }
        return data[j*M+i];
    }

    Matrix &operator *(double rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat * scalar");
        Matrix &  X = * new Matrix(M,N);
        for(int i=0; i<M; i++)
        {
            for(int j=0; j<N; j++)
            {
                X(i+1,j+1) = data[j*M+i] * rhs;
            }
        }

        return X;
    }
    friend Matrix & operator *(double d, Matrix & M )
    {
        return M.operator *(d);
    }
    Matrix &operator +(double rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat + scalar");
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1)=data[j*M+i] + rhs;
        return X;
    }
    Matrix &operator -(double rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat - scalar");
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1)=data[j*M+i] - rhs;
        return X;
    }
    friend Matrix & operator +(double d, Matrix & M )
    {
        return M.operator +(d);
    }
    Matrix & operator +(Matrix &rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat + Mat");
        Matrix &X = * new Matrix(M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X.data[j*M+i] = data[j*M+i] + rhs.data[j*M+i];
        return X;
    }
    Matrix & operator *(Matrix &rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat * Mat ");
        /*ToDo: blas gemm*/
        int K = rhs.cols();
        Matrix *X= new Matrix(M,K);
        for(int i=0; i<M; i++)
        {
            for(int j=0; j<N; j++)
            {
                X->data[j*M+i]=0;
                for(int k=0; k<K; k++)
                {
                    X->data[j*M+i] += data[k*M+i] * rhs.data[j*M+k];
                }
            }
        }
        return *X;

    }
    Matrix & operator -(Matrix &rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat - Mat");
        Matrix *X = new Matrix(M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X->data[j*M+i] = data[j*M+i] - rhs.data[j*M+i];
        return *X;
    }
    Matrix & operator ^(double rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat ^ scalar");
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1)= pow(data[j*M+i],rhs);
        return X;
    }
    Matrix & sqrt()
    {
        if ( locked)
            return do_lock_op<Matrix&>("sqrt ( Mat )");
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1) = std::sqrt((double)data[j*M+i]);
        return X;
    }
    Matrix & log()
    {
        if ( locked)
            return do_lock_op<Matrix&>(" log(Mat)");
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1) = std::log((double)data[j*M+i]);
        return X;
    }

    Matrix & transpose(bool in_place = true)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat^T");
        if ( in_place)
        {
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
            {
                ElementType t = data[j*M+i];
                data[j*M+i] = data[i*M+j];
                data[i*M+j]=t;
            }
            return *this;
        }
        Matrix &X= * new Matrix(N,M);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X.data[i*M+j]=data[j*M+i];
        return X;
    }
    Matrix &column(int c)
    {
        Matrix &X = * new Matrix(&data[c*M],M,1);
        if ( locked)
            X.lock();
        X.set_parent(this);
        return X;
    }
    int rows(){return M;}
    int cols()
    {
        return N;
    }
    Matrix &diag()
    {
        Matrix &X=*new Matrix(M,1);
        if(locked)
            X.lock();
        X.set_parent(this);
        for ( int i=0; i<M; i++)//TOdo : do not create a new copy,use the same memory with proper LDA
            X.data[i]=data[i*M+i];
        return X;
    }
    Matrix &operator >(ElementType rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat > scalar");
        Matrix &X= * new Matrix(M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
            X.data[j*M+i] = (data[j*M+i] > rhs)?1:0;
        return X;
    }
    int numel()
    {
        return M*N;
    }
    Matrix &operator /(double rhs)
    {
        if ( locked)
            return do_lock_op<Matrix&>("Mat / scalar");
        assert ( fabs(rhs) > 1e-100);
        Matrix &X=*new Matrix (M,N);
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                X(i+1,j+1)=data[j*M+i] / rhs;
        return X;
    }
    friend Matrix &operator /(ElementType lhs, Matrix & rhs)
    {
        if ( rhs.is_locked()){
            rhs.do_lock_op2("scalar / Mat");
            return rhs;
        }
        Matrix &X=*new Matrix (rhs.rows(),rhs.cols());
        for ( int i=0; i<rhs.rows(); i++)
            for ( int j=0; j<rhs.cols(); j++)
                X(i+1,j+1) = lhs / rhs(i+1,j+1);
        return X;
    }
    Matrix &col_slice(int fr, int to)
    {
        int len=to-fr+1;
        Matrix &X = * new Matrix(M,len);
        if(locked)
            X.lock();
        X.set_parent(this);
        for(int i=0; i<M; i++)
        {
            for ( int j=fr; j < to; j++)
            {
                int c = j%(len);
                X.data[c*M+i]=data[j*M+i];
            }

        }
        return X;
    }
    Matrix & row_slice(int fr, int to)
    {
        int len = to - fr + 1;
        Matrix &X= * new Matrix(len,N);
        if(locked)
            X.lock();
        X.set_parent(this);
        for ( int i=fr; i<to; i++)
        {
            int r = i%(len);
            for ( int j =0; j<N; j++)
            {
                X.data[j*M+r] = data[j*M+i];
            }
        }
        return X;
    }
    Matrix() {data = nullptr;}
    Matrix & pick_rows(Matrix &P)
    {
        int m = 0;
        if ( P.cols() ==1)
        {
            m = P.rows();
        }
        if ( P.rows() ==1)
        {
            m = P.cols();
        }
        if ( m ==0)
        {
            fprintf(stderr,"row_pick is called with non  vector input\n");
            exit(-1);
        }
        if ( M==1 && N>1)
        {
            Matrix &X = *new Matrix (m,1);
            if(locked)
                X.lock();
            for ( int i=1; i<=m; i++)
            {
                int c=P(i);
                X(i,1)=data[c];
            }
            return X;
        }
        Matrix  &X= *new Matrix (m,N);
        if(locked)
            X.lock();
        for ( int i=1; i<=m; i++)
        {
            int r=P(i);
            for (int j=0; j<N; j++)
                X(i,j+1)=data[j*M+r];
        }
        return X;
    }
    Matrix &round()
    {
        if ( locked)
            return do_lock_op<Matrix&>("round(Mat)");
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                data[j*M+i] = std::round(data[j*M+i]);
        return *this;
    }
    Matrix(ElementType *d,int M_,int N_)
    {
        data=d;
        M=M_;
        N=N_;
    }
    void do_col_partition( int c)
    {
        pcols.clear();
        int len = N/c;
        for ( int j= 0; j<c; j++)
        {
            Matrix *m = new Matrix(&data[j*len*M],M,1);
            pcols.push_back(m);
        }
    }
    Matrix *get_part_column(int c)
    {
        return pcols[c-1];
    }
    void rebuild(int m, int n){
        assert(m*n !=0);
        if (data)
            delete [] data;
        M=m;
        N=n;
        data = new ElementType[M*N];
    }
    Matrix & set_column( int c, Matrix &rhs,bool build=false);
    Matrix & append(double v);
    Matrix &get_row(int r);
    Matrix &get_column(int c);
    Matrix &max(int dim);
    Matrix &min(int dim);
    Matrix &floor();
    Matrix &diff();
    Matrix &find();
    ElementType *get_data_memory(){return data;}
    void print(const char *s="");
    void export_data(const char *fn);
};

Matrix &repmat ( Matrix &,int,int);
Matrix &sqrt(Matrix &);
Matrix &round(Matrix  &);
Matrix &floor(Matrix  &);
Matrix &log_kernel(Matrix &);
int     size (Matrix &M, int rc);
Matrix &diag (Matrix &M);
int     numel(Matrix &M);
int length(Matrix &);
void    SVD  (Matrix &A, Matrix &U, Matrix &S, Matrix &V);
Matrix & zeros( int r,int col=1);
Matrix &max(Matrix &M, int dim=1);
Matrix &min(Matrix &M, int dim=1);
Matrix &diff(Matrix &);
Matrix &find(Matrix &);
Matrix &hcat(Matrix &, Matrix &);

#endif // MATRIX_HPP_INCLUDED
