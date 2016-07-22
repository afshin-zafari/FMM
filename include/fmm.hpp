#ifndef FMM_HPP_INCLUDED
#define FMM_HPP_INCLUDED

#include <iostream>
#include <vector>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <cassert>
#include <sys/time.h>
#include "acml.h"
#include "matrix.hpp"
#include "sgmatrix.hpp"
#include "sg/superglue.hpp"
#include "sg/option/instr_trace.hpp"

#define M_PI		3.14159265358979323846

#undef FMM_3D
#define svd_threshold 1.0e-2

using namespace std;
const int Q = 10;
typedef Matrix& (*KernelFcn)(Matrix &);


template <typename T>
class elastic_vect{
public:
    vector<T> v;
    int size(){return v.size();}
    void clear(){v.clear();}
    T& operator()(unsigned int i){
        i--;
        while(v.size()<=i)
            v.push_back(* new T);
        return v[i];
    }
    operator Matrix&()const{
        int N = v.size();
        Matrix &X= * new Matrix(1,N);
        for ( int i=1;i<=N;i++){
            X(1,i) = v[i-1];
        }
        return X;
    }
    void append(T e){
        v.push_back(*new T(e));

    }
    elastic_vect<T> &operator -(elastic_vect<T> &rhs){
        elastic_vect<T> &X = * new elastic_vect<T>;
        for ( unsigned int i=1;i<=v.size();i++)
            X(i) = v[i-1] - rhs(i);
        return X;
    }
    elastic_vect<T> &operator /(double rhs){
        elastic_vect<T> &X = * new elastic_vect<T>;
        for ( unsigned int i=1;i<=v.size();i++)
            X(i) = v[i-1] / rhs;
        return X;
    }
    elastic_vect<T> &operator +(double rhs){
        elastic_vect<T> &X = * new elastic_vect<T>;
        for ( unsigned int i=1;i<=v.size();i++)
            X(i) = v[i-1] + rhs;
        return X;
    }
    elastic_vect<T> & sign(){
        elastic_vect<T> &X= * new elastic_vect<T>;
        for ( int i=1;i<=size();i++)
            X(i) = sgn(v[i-1]);
        return X;
    }
    elastic_vect<T> & abs(){
        elastic_vect<T> &X= * new elastic_vect<T>;
        for ( int i=1;i<=size();i++)
            X(i) = fabs(v[i-1]);
        return X;
    }
    double max(){
        double m = -1e200;
        for ( int i=0;i<size();i++)
            if ( m<v[i])
                m=v[i];
        return m;
    }
    elastic_vect<T> &operator =(elastic_vect<T> & rhs)
    {
        v.clear();
        for ( int j=1;j<=rhs.size();j++)
            v.push_back(rhs(j));
        return *this;
    }
    elastic_vect<T> &operator =(Matrix &M){
        int m = M.rows();
        int n = M.cols();
        assert(m==1 || n==1);
        v.clear();
        if( m==1 ){
            for ( int j=1;j<=n;j++)
                v.push_back(M(1,j));
            return *this;
        }
        if( n==1 ){
            for ( int i=1;i<=m;i++)
                v.push_back(M(i,1));
        }
        return *this;

    }
};
typedef elastic_vect<double> TempDoubleVector;
typedef elastic_vect<int> TempIntVector;
TempDoubleVector &sign(TempDoubleVector &t);
TempDoubleVector &abs(TempDoubleVector &t);
double max(TempDoubleVector &t);

Matrix &round(TempDoubleVector &t);
/*---------------------------------------*/
template<typename T>
class cell{
private:
    vector<T*> content;
    int M,N;
    bool zero_indexing;
public:
    cell(int a, int b):M(a),N(b){
        for(int i=0;i<M;i++)
            for(int j=0;j<N;j++)
                content.push_back( new T);
        zero_indexing = false;
    }

    T &operator ()(int i,int j){
        if ( !zero_indexing){i--;j--;}
        T *a = static_cast<T*>(NULL);;
        if ( i>=M)
            return *a;
        if (j>=N)
            return *a;
        return *content[i*N+j];
    }
};
typedef int GroupId;
typedef double PointType;
class GroupType{
public:
    elastic_vect<GroupId> child,neighboursList,intList;
    elastic_vect<double> groupcenter;
    SGMatrix V;
    PointType *x,*y,*z;
    int num_point;
    double cubelength;
    void print_lists(){

        cout << "\t\t children(" <<child.size()<<"):\t";
        for ( auto c : child.v)
            cout << c << ",";
        cout << "\n\t\t interact:\t";
        for ( auto i : intList.v)
            cout << i << ",";
        cout << "\n\t\t ngbr:\t\t";
        for ( auto n : neighboursList.v)
            cout << n << ",";
        cout << "\n\t\t CubeLength:\t" << cubelength ;
        cout << "\n\t\t GroupCenter:\t"
            << groupcenter(1) << ' '
            << groupcenter(2) << ' '
            << groupcenter(3) << '\n';
        cout << endl;

    }
};
class Level{
public:
    elastic_vect<GroupType> group;
    cell<SGMatrix> M2M,T;
    SGMatrix AuxSources,AuxObs;

    Level():M2M(2,2),T(7,7){
    }

    void after_group_built(){
        AuxSources.build(Q,group.size(),SGMatrix::Column);
        AuxObs    .build(Q,group.size(),SGMatrix::Column);
    };
    void print_list(){
        for ( int i=1; i<=group.size();i++){
            cout << "\t grp: " << i << endl;
            group(i).print_lists();
        }
    }
};
struct NearFieldBlock{
public:
    int from_row,to_row,from_col,to_col;
    int i,j,ni,nj,p;
    Matrix *M;
    SGMatrix *sgMat;
    SGHandle handle;

    NearFieldBlock(int ii,int jj,int nni,int nnj,int pp,Matrix *MM):
        i(ii),j(jj),ni(nni),nj(nnj),p(pp),M(MM){
            assert(M);
            sgMat = new SGMatrix(*M);
    }
};
struct Partition {
public:
    Partition(int i, int _n, double *m):index(i),n(_n){        
        memory = m;
    }
    int index,n;
    double *memory;
    SGHandle handle;
};
typedef elastic_vect<Partition*> PartitionedVector;
class Tree{
public:
    int Q;
    char method;
    elastic_vect<Level> Levels;
    vector<NearFieldBlock*> NearField;

    Tree(){}
    int num_levels(){return Levels.size();}
    Level &operator ()(int i){
        return Levels(i);
    }
    void print_lists(){
        for ( int l=1;l<=Levels.size();l++){
            cout <<"Level: " << l << endl;
            Levels(l).print_list();
        }
        cout << endl;

    }
	~Tree(){
		fprintf(stdout,"~Tree\n");
	}
};
class MvTask{
public:
    MvTask(SGMatrix &M,SGMatrix &S,int column,SGVector &y){}
};

void cblas_dgemv(const int layout,
                 const bool TransA,
                 const int M, const int N,
                 const double alpha, const double *A, const int lda,
                 const double *X, const int incX,
                 const double beta, double *Y, const int incY);

class SGTaskGemv : public Task<Options, 3> {
private:
    SGMatrix *A,*v,*y;
	string name;
public:
	Time::TimeUnit  s,r,f,d;
	int type;
    bool transA,near_field;
    enum{
		Read=ReadWriteAdd::read,
		Write=ReadWriteAdd::write,
		Add=ReadWriteAdd::add};
    enum{COL_MAJOR,ROW_MAJOR};	
    SGTaskGemv(SGMatrix &A_, SGMatrix &v_, int i1,SGMatrix &Y_, int i2){
        A = &A_;
        v = &v_.get_part(i1);
        y = &Y_.get_part(i2);
		near_field = false;
        register_args();
    }
    SGTaskGemv(SGMatrix &A_, SGMatrix &v_, SGMatrix &Y_){
        A = &A_;
        v = &v_;
        y = &Y_;
		near_field = true;
        register_args();
    }
    void register_args(){
        SGHandle &hA = A->get_handle();
        SGHandle &hv = v->get_handle();
        SGHandle &hy = y->get_handle();
		if (config.h){
			register_access(ReadWriteAdd::read, hA);
			register_access(ReadWriteAdd::read, hv);
		}
		else
			if (!near_field)
				register_access(ReadWriteAdd::read, hv);
        if ( config.w )
            register_access(ReadWriteAdd::write, hy);
        else
            register_access(ReadWriteAdd::add, hy);
        transA = false;

    }
    //void register_access(int axs, SGHandle &h){}
    void run(){
		r = Time::getTime();
		assert(A->get_matrix());
		assert(v->get_matrix());
		assert(y->get_matrix());
        const int M = A->get_matrix()->rows();
        const int N = A->get_matrix()->cols();
        const int lda = M;
        double *Mat = A->get_matrix()->get_data_memory();
        double *X   = v->get_matrix()->get_data_memory();
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
		if ( !config.x)
	        cblas_dgemv(COL_MAJOR,transA,M, N, 1.0, Mat, lda, X, 1, 1.0, Y, 1);
		///if ( !config.x)			dgemv(transA?'T':'N',M, N, 1.0, Mat, lda, X, 1, 1.0, Y, 1);
		f = Time::getTime();
    }
	~SGTaskGemv(){
		d = Time::getTime();
		fprintf(stdout,"%ld, %ld, %ld, %ld, %d\n",s,r,f,d,type);
	}
	std::string get_name(){return name;}
	void set_name(const char *s){name=s;}
};
void submit_all(void);


extern timeval start,finish;
void tic();
double toc();


struct Statistics{
	int t,l,g;
	double dur[6];
};
extern Statistics stats;
template<typename T>
int length(elastic_vect<T> v);
int length(Tree &T);

void unique(Matrix &,Matrix &, TempIntVector &);
void sort(Matrix &, Matrix&);

Tree & BuildOctree(Matrix &pts, double MinCubeSide,Tree &Octree);
void BuildInteractionLists(Tree &OT,int iLev,int iGroup,int jGroup);
void BuildInteractionLists(Tree &OT);
void compute_translator(Tree &OT,KernelFcn kernel);
void compute_radiation( Tree &OT, Matrix &pts, char method, KernelFcn kernel);
void compute_near_field(Tree &OT,Matrix &pts);
void mv_near_field(Tree &OT,SGMatrix &C, SGMatrix &Q);
void MatVec(  Tree & OT,SGMatrix &x , SGMatrix &y);

class EventLog{
public:
	Time::TimeUnit start;
	string s;
	EventLog(const char * s_):s(s_){
		start=Time::getTime();
	}
	void End(){
		Time::TimeUnit end = Time::getTime();
		double dur = end - start;
		fprintf(stdout,"%s: %lf (s).\n",s.c_str(), dur /3e9);
		start = 0;
	}

	~EventLog(){
		if ( start)
			End();
	}
};
extern double submit_time;
#endif // FMM_HPP_INCLUDED
