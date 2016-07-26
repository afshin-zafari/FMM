#ifndef XTYPE_HPP_INCLUDED
#define XTYPE_HPP_INCLUDED


#include <iostream>
#include <typeinfo>
#include <type_traits>
using namespace std;
//#define TL fprintf(stdout,"%20.20s, %20.20s, %d\n",__FILE__,__FUNCTION__,__LINE__);
//#define TL cout << __FILE__<< ' ' << __FUNCTION__<< " " << __LINE__ << endl 
#define TL 

class XType
{
protected:
    bool   plain,non_xtype;
    XType  *left,*right;
    char   op;
    string name;
    int    level,h;
    double val;
    void*ex;
public:
    XType(string s,int hd=-1):name(s),h(hd)
    {
        left = right = NULL;
        op=' ';
        plain=true;
        level=0;
        non_xtype = false;
        ex=NULL;
    }
    XType(const char *p, void *pp)
    {
        name.assign(p);
        ex=pp;
        //cout << *(long *)ex << endl;
        left = right = NULL;
        op=' ';
        plain=true;
        level=0;
        non_xtype = false;
    }
    XType(const char * n,double v)
    {
        name.assign(n);
        val=(double)v;
        non_xtype = true;
        plain=true;
        h=-2;
        ex=NULL;
    }
    XType(XType &l,const char *opr, XType &r)
    {
        left = &l;
        right = &r;
        op = opr[0];
        op= (op =='a')?'+':'*';
        op= (opr[1]=='s')?'=':op;
        name.assign(opr);
        plain = false;
        level = l.levelx() +1;
        non_xtype = false;
        h=-3;
        ex=NULL;
    }
    XType & operator *(XType & rhs)
    {
        const char *mult="mult";
        XType &x=* new XType(*this,(const char *)mult,rhs);
        //cout << name << " ," << mult <<"," << rhs.namex() << endl;;
        return x;
    }
    XType & operator *(double  rhs)
    {
        const char *op="mult";
        XType *y = new XType("non-xtype",rhs);
        XType &x=* new XType(*this,(const char *)op,*y);
        //cout << name << ":: ," << op<<"," << y->namex() << endl;;
        return x;
    }
    friend XType & operator *(double  lhs,XType &x )
    {
        return x.operator *(lhs);
    }

    XType & operator +(XType & rhs)
    {
        const char *op="add";
        XType &x=* new XType(*this,(const char *)op,rhs);
        //cout << name << " ," << op<<"," << rhs.namex() << endl;;
        return x;
    }
    XType & operator =(XType & rhs)
    {
        const char *op="asgn";
        XType &x=* new XType(*this,(const char *)op,rhs);
        cout << name << " ," << op<<"," << rhs.namex() << endl;
        cout << "+++++++++++++++++++++++++++++++\n";
        cout << x << endl;
        cout << "+++++++++++++++++++++++++++++++\n";
        return x;
    }

    friend XType & operator +(double  lhs,XType &x )
    {
        return x.operator +(lhs);
    }

    XType & operator +(double   rhs)
    {
        const char *op="add";
        XType *y = new XType("non-xtype",rhs);
        XType &x=* new XType(*this,(const char *)op,*y);
        //cout << name << " ," << op<<"," << y->namex() << endl;;
        return x;
    }
    XType &leftx()
    {
        return *left;
    }
    XType &rightx()
    {
        return *right;
    }
    int levelx()
    {
        return level;
    }
    string namex()
    {
        return name;
    }
    char opx()
    {
        return op;
    }
    double value()
    {
        return val;
    }
    void *hx()
    {
        return ex;
    }

    friend ostream & operator << (ostream &o,XType &x)
    {
        o << "{";
        o << x.namex() ;
        if ( !x.plain)
        {
            o << ',' << x.opx() << ',';
            o << x.leftx() << ',';
            o <<x.rightx() ;
        }
        else if ( x.non_xtype)
            o << ',' << x.value();
        o << "}";
        return o;

    }

    bool is_plain()
    {
        return plain;
    }
    XType*get_right_node()
    {
        return right;
    }
    XType *get_left_node()
    {
        return left;
    }
	static int LastHandle;
};
class YType: public XType
{
private:
    long h;
public:
    YType(const char *p):XType(p,(void *)&h)
    {
        h = ++LastHandle;
        cout << &h << endl;
    }
    //using XType::operator=;
    YType & operator =(XType & rhs)
    {
        const char *op="asgn";
        XType &x=* new XType(*this,(const char *)op,rhs);
        cout << name << " ," << op<<"," << rhs.namex() << endl;
        cout << "================================\n";
        cout << x << endl;
        cout << "================================\n";
        return static_cast<YType &>(x);
    }


};

#endif // XTYPE_HPP_INCLUDED
