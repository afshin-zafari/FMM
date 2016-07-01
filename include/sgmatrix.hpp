#ifndef SGMATRIX_HPP_INCLUDED
#define SGMATRIX_HPP_INCLUDED

#include "xtype.hpp"
typedef unsigned long Handle;
typedef unsigned int uint;

class SGMatrix: public XType
{
    Handle *sg_handle;
    Matrix *M;
    vector<SGMatrix*> pcols,prows,parts;
    bool trans;
    int pM,pN;
public:
    //using XType::operator=;
    SGMatrix(const char *p):XType(p,(void *)this)
    {
        sg_handle  = new Handle(++LastHandle);
        trans= false;
        M=NULL;
        //cout << "last-h: "<< LastHandle << endl;
        //cout << "  sg-h: "<< *sg_handle << endl;
    }
    SGMatrix(Matrix &m):XType("",(void *)this){
        sg_handle = new Handle(++LastHandle);
        trans=false;
        M = &m;
    }
    SGMatrix & operator =(XType & rhs)
    {
        const char *op="asgn";
        XType &x=* new XType(*this,(const char *)op,rhs);
        //cout << name << " ," << op<<"," << rhs.namex() << endl;
        cout << "================================\n";
        cout << x << endl;

        cout << "================================\n";
        SGMatrix &S = static_cast<SGMatrix&>(x);
        S.interpret();
        return S;
    }

    void interpret()
    {
        if ( op != '=')
            return;
        if ( !left->is_plain())
            return;
        if ( right->opx() != '+')
            return;
        XType *rr = right->get_right_node();
        if ( !rr->is_plain() && rr->opx()!='*')
            return;
        XType *rl = right->get_left_node();
        if ( !rl->is_plain() && rl->opx()!='*')
            return;
        if(rr->is_plain() && rl->opx() !='*')
            return;
        if(rl->is_plain() && rr->opx() !='*')
            return;
        if (!rl->is_plain() && !rr->is_plain())
            return;
        XType *mult;
        if ( rr->opx()=='*')
            mult=rr;
        else
            mult=rl;
        if (!mult->get_left_node()->is_plain())
            return;
        if (!mult->get_right_node()->is_plain())
            return;
        SGMatrix *lfm=(SGMatrix *)left->hx();
        Handle lfh = lfm->get_handle();
        if (rr->is_plain())
        {
            //cout << "a=bc+a" << endl;
            SGMatrix *rrm=(SGMatrix *)rr->hx();
            Handle rrh = rrm->get_handle();
            if ( rrh != lfh)
                return;

        }
        else
        {
            //cout << "a=a+bc" << endl;
            //cout << rl->namex() << endl;
            SGMatrix *rlm=(SGMatrix *)rl->hx();
            Handle rlh = rlm->get_handle();
            if ( rlh != lfh)
                return;
        }
        cout << "gemm/gemv comopatible\n";
        Matrix *lhs = lfm->get_matrix();
        if ( !lhs)
            return;
        if ( lhs->rows() ==1 || lhs->cols() ==1)
            cout << "gemv task\n";
        else
            cout << "gemm task\n";
    }
    void setMatrix(Matrix &A){
        M=&A;
    }
    Matrix *get_matrix()
    {
        return M;
    }
    Handle get_handle()
    {
        return *sg_handle;
    }
    SGMatrix():XType("X",(void *)this)
    {
        sg_handle = new Handle(++LastHandle);
        trans=false;
    }
    typedef enum Partition {Column,Row} Partition;
    void build(int rows, int cols, Partition p)
    {
        M = new Matrix(rows,cols,p);
        if ( p == SGMatrix::Column)
            do_col_partition(cols);
    }
    void operator =(Matrix &rhs)
    {
        if ( !M){
            M = &rhs;
            rhs.lock();
        }
        else
            //ToDo What to do, if Matrix is overwritten?
            ;
    }

    void do_col_partition ( int c)
    {
        pM = 1;
        pN = c;
        pcols.clear();
        assert(M!=NULL);
        int cols=M->cols();
        assert(c>0 && c<=cols);
        M->do_col_partition(c);
        while(c--)
        {
            SGMatrix *m=new SGMatrix;
            m->M = M->get_part_column(c+1);
            pcols.push_back(m);
        }
    }

    SGMatrix &column(int c)
    {
        return *pcols[c-1];
    }

    SGMatrix &transpose()
    {
        SGMatrix &X= *new SGMatrix;
        X.trans= !trans;
        return X;
    }
    SGMatrix &get_part(uint i){
        i--;
        assert ( i>=0);
        assert ( i< parts.size());
        return *parts[i];
    }
    void set_part_size(int n){
        do{
            parts.push_back((SGMatrix *)nullptr);
        }while(n--);

    }
    void set_part(int i, int m ,  int n, double *mem){
        SGMatrix *MM =  new SGMatrix;
        Matrix *M = new Matrix(m,n,mem);
        MM->set_mat(M);
        parts[i-1]=MM;
    }
    void set_mat(Matrix *M_){
            M = M_;
    }
};
class SGVector : public SGMatrix
{
public:
    SGVector():SGMatrix() {}
};



#endif // SGMATRIX_HPP_INCLUDED
