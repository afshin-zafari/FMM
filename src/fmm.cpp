#include <algorithm>
#include "fmm.hpp"

template<typename T>
int length(elastic_vect<T> v){
    return (int)v.size();
}
template<>
int length(elastic_vect<GroupType> v){
    return (int)v.size();
}
template<>
int length(elastic_vect<GroupId> v){
    return (int)v.size();
}

int length(Tree &T){
    return T.num_levels();
}
TempDoubleVector &sign(TempDoubleVector &t){
    return t.sign();
}
TempDoubleVector &abs(TempDoubleVector &t){
    return t.abs();
}
double max(TempDoubleVector &t){
    return t.max();
}
Matrix &round(TempDoubleVector &t){
    Matrix &X = * new Matrix(1,t.size());
    for ( int i=1;i<=t.size();i++)
        X(1,i)=round(t(i));
    return X;
}
void unique(Matrix &in, Matrix &out, TempIntVector &indx){
    int m = in.rows();
    int n = in.cols();
    assert(m==1 || n==1);
    int len=m*n;
    indx.clear();
    int j=0;
    vector<ElementType> res;
        vector<ElementType>::iterator v(in.get_data_memory());
        for( int i=0;i<len;i++,v++){
            ElementType e = *v;
            vector<ElementType>::iterator p=std::find(res.begin(),res.end(),e);
            bool found = (p !=res.end());
            if (!found ){
                res.push_back(e);
                indx(++j)=i+1;
            }
        }
    if ( n==1){
        out = *new Matrix(res.size(),1);
        for ( unsigned int i=1;i<=res.size();i++)
            out(i,1)=res[i-1];
    }
    if ( m==1){
        out = *new Matrix(1,res.size());
        for ( unsigned int i=1;i<=res.size();i++)
            out(1,i)=res[i-1];
    }
    return;

}
void sort(Matrix &in, Matrix &indx){
    int m= in.rows();
    int n = in.cols();
    assert ( n ==1 || m == 1);
    indx.clear();
    int k = n*m;
    for ( int i =1;i<=k;i++)
        indx(i)=i;
    ElementType t;
    #define SWAP_Val(i,j) t=in(i);in(i)=in(j);in(j)=t;
    #define SWAP_Idx(i,j) t=indx(i);indx(i)=indx(j);indx(j)=t;
    for ( int i=1;i<=k;i++){
        for ( int j=i+1;j<=k;j++){
            if (in(i)>in(j)){
                SWAP_Val(i,j);
                SWAP_Idx(i,j);
                if ( indx(i) ==0){
                    t=0;
                }
            }
        }
    }
}
timeval start,finish;
int XType::LastHandle=0;

void tic()
{
    gettimeofday(&start,NULL);
}
double toc()
{
    gettimeofday(&finish,NULL);
    return ((finish.tv_sec-start.tv_sec)*1000000+finish.tv_usec-start.tv_usec)/1000000.0;
}
