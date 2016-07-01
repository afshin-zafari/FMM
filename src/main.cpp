#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>

#include "fmm.hpp"
using namespace std;

void test2(){
    SGMatrix A("A"),B("B"),C("C"),D("D"),X;
    Matrix a(10,10),v(10,1);
    C=v;
    D=a;
    C= A*B+2.*X;
    C=D+A*B;
    C=C+A*B;

}
void test1(){
    Matrix A(10,10),B(10,20),C(10,20),D(10,20);
    D=A*B+C.get_column(1);

}
int N=100;
int L = 3;
void init(Matrix &pts){
    int M = N;
    double theta = 0.0, theta_mod =0.0;
    cout << "Points" << endl;
    for ( int i=1; i<=M;i++){
            pts(i,1) = 2* cos(theta)*(1.0 + .1*cos(theta_mod));
            pts(i,2) = 1* sin(theta)*(1.0 + .1*cos(theta_mod));
            pts(i,3)=0.0;
            theta += 2.*M_PI /N;
            theta_mod += 32*M_PI /N;
//            cout << pts(i,1) << ',' << pts(i,2) << endl;
    }
}
void test_svd(){
    const int M = 6,_N=5;
    double a[M*_N] = {
        8.79,  9.93,  9.83, 5.45,  3.16,
        6.11,  6.91,  5.04, -0.27,  7.98,
       -9.15, -7.93,  4.86, 4.85,  3.01,
        9.57,  1.64,  8.83, 0.74,  5.80,
       -3.49,  4.02,  9.80, 10.00,  4.27,
        9.84,  0.15, -8.99, -6.02, -5.31
    };
    Matrix A(M,_N);
    double *d=A.get_data_memory();
    for ( int i=0;i<M;i++){
        for (int j=0;j<_N;j++){
            d[i+j*M]=a[i*_N+j];
//            printf("%2.2lf ",d[i+j*M]);
            //printf("%2.2lf ",a[i*N+j]);
            cout << setprecision (5) << setw(7) << a[i*_N+j] << ',';
        }
        printf("\n");
    }
    Matrix U,S,V;
    A.print();
    SVD(A,U,S,V);
    S.print();
    U.print();
    V.print();
}
void test_main(){
    Matrix &pts = *new Matrix (N,3);
    Tree &OT=*new Tree;
    OT.Q = 10;
    init(pts);
    double max_dim=max(max(pts) - min(pts))(1,1);
    double min_cube_size=max_dim *pow(2,1-L);
    BuildOctree(pts,min_cube_size,OT);
    //BuildInteractionLists(OT);
    OT.print_lists();
    //compute_translator(OT,log_kernel);
    //compute_radiation (OT, pts, 'T', log_kernel);
    int i=0;
    printf("%d\n",i);
}

class Reader{
public:
    string center_token,
        cubelength_token,
        level_token,
        group_token,
        intact_token,
        nbr_token,
        m2m_cell_token,
        m2m_token,
        xlat_cell_token,
        xlat_token;
    std::ifstream f;
    enum {
        level_type,      //0
        group_type,
        intact_type,
        nbr_type,
        center_type,
        cubelength_type, //5
        m2m_cell_type,
        m2m_type,
        xlat_cell_type,
        xlat_type};
    std::string::size_type pos;
    string f1,f2;
    int level,group;
    Tree *t;
    Reader(string s,string s2,Tree &T){
        f.open(s);
        f2 = s2;
        level = group = -1;
        t = &T;
        center_token    = string("Center");
        cubelength_token= string("CubeLength");
        level_token     = string("Level");
        group_token     = string("Group");
        intact_token    = string("Int Act");
        nbr_token       = string("Neighbor");
        m2m_cell_token  = string("M2M-Cell");
        m2m_token       = string("M2M");
        xlat_cell_token = string("Translator-Cell");
        xlat_token      = string("Translator");
    }
    int get_type(string s){
        pos = s.find(level_token);
        if ( pos != std::string::npos )
            return level_type;
        pos = s.find(group_token);
        if ( pos != std::string::npos)
            return group_type;
        pos = s.find(intact_token);
        if ( pos != std::string::npos )
            return intact_type;
        pos = s.find(nbr_token);
        if ( pos != std::string::npos)
            return nbr_type;
        pos = s.find(center_token);
        if ( pos != std::string::npos)
            return center_type;
        pos = s.find(cubelength_token);
        if ( pos != std::string::npos)
            return cubelength_type;
        pos = s.find(m2m_cell_token);
        if ( pos != std::string::npos)
            return m2m_cell_type;
        pos = s.find(m2m_token);
        if ( pos != std::string::npos)
            return m2m_type;
        pos = s.find(xlat_cell_token);
        if ( pos != std::string::npos)
            return xlat_cell_type;
        pos = s.find(xlat_token);
        if ( pos != std::string::npos)
            return xlat_type;

        return -1;
    }
    void parse(string s,int type)
    {
        string v;
        Tree &T=*t;
        std::string::size_type p1=0,p2;
        while ( true ) {
            p2 = s.find(",",p1);
            if ( p2 == string::npos)
                return;
            v =  s.substr(p1,p2-p1) ;
            cout << v << " , " ;
            switch(type)
            {
                case group_type:
                    T.Levels(level).group(group).child.append(atoi(v.c_str()));
                    break;
                case intact_type:
                    T.Levels(level).group(group).intList.append(atoi(v.c_str()));
                    break;
                case nbr_type:
                    T.Levels(level).group(group).neighboursList.append(atoi(v.c_str()));
                    break;
                default:
                    cout << "Bad Type: " << type << endl;
                    exit(-1);
                    break;
            }
            p1 = p2+1;
        }
    }
    void parse_group_info(string line)
    {
        double v;
        Tree &T=*t;
        string::size_type c;
        string s;
        for ( int i=0; i< 3; i++)
        {
            c=line.find(center_token);
            s=line.substr(c+center_token.size());
            v = atof(s.c_str());
            T.Levels(level).group(group).groupcenter.append(v);
            getline(f,line);
        }
        c=line.find(cubelength_token);
        s=line.substr(c+cubelength_token.size());
        v = atof(s.c_str());
        T.Levels(level).group(group).cubelength=v;


    }
    void parse_v(int g,int r, int c){
        Tree &T=*t;
        Matrix *M=parse_mat(r,c);
        T.Levels(level).group(g).V.setMatrix(*M);
        M->print();
    }
    double get_double(string &s){
        string::size_type p1=0,p2;
        p2 = s.find(",",p1);
        double d=atof(s.substr(p1,p2-p1).c_str());
        s=s.substr(p2+1);
        return d;
    }
    Matrix *parse_mat(int r , int c)
    {
        string s;
        if ( r*c ==0 )
            r=1;
        Matrix &M=*new Matrix(r,c);
        for ( int k=1;k<=r;k++){
            getline(f,s);
            for ( int l=1;l<=c;l++){
                M(k,l)=get_double(s);
            }
        }
//        M.print();
        return &M;
    }
    void parse_m2m(int i,int j,int r,int c)
    {
        Tree &T=*t;
        Matrix *M=parse_mat(r,c);
        T.Levels(level).M2M(i,j).setMatrix(*M);

    }
    void parse_xlat(int i,int j,int r,int c){
        Tree &T=*t;
        Matrix *M=parse_mat(r,c);
        T.Levels(level).T(i,j).setMatrix(*M);
    }
    void read(){
        string s,g;
        string::size_type p ;
        while ( !f.eof())
        {
            getline(f,s);
            int type = get_type(s);
            switch(type)
            {
            case level_type: // Level line
                level = atoi(s.c_str() + level_token.size());
                break;
            case group_type:
                p = s.find(group_token);
                g = s.substr(p+group_token.size());
                group = atoi(g.c_str());
                getline(f,s);
                if (get_type(s) == center_type){
                    parse_group_info(s);
                    getline(f,s);
                }
                parse(s,type);
                break;
            case intact_type:
                getline(f,s);
                parse(s,type);
                break;
            case nbr_type:
                getline(f,s);
                parse(s,type);
                break;

            }
            cout << "Level: " << level <<" Group: " << group << endl;
        }
        cout << "-----------------------\n";
        t->print_lists();
    }
    void read_op(){
        string s,g;
        string::size_type p ;
        int i,j,r,c;
        char temp[10];
        f.close();
        f.open(f2);
        while ( !f.eof())
        {
            getline(f,s);
            int type = get_type(s);
            switch(type)
            {
            case level_type: // Level line
                level = atoi(s.c_str() + level_token.size());
                break;
            case m2m_cell_type:
                r=c=-1;
                sscanf(s.c_str(),"%s %d %d %d %d",temp,&i,&j,&r,&c);
                if (i==0)
                    break;
                break;
            case m2m_type:
//                if(s.find("Cell")!= string::npos)
//                    getline(f,s);
                sscanf(s.c_str(),"%s %d %d %d %d",temp,&i,&j,&r,&c);
                parse_m2m(i,j,r,c);
                break;
            case xlat_cell_type:
                r=c=-1;
                sscanf(s.c_str(),"%s %d %d %d %d",temp,&i,&j,&r,&c);
                break;
            case xlat_type:
//                if(s.find("Cell")!= string::npos)
//                    getline(f,s);
                sscanf(s.c_str(),"%s %d %d %d %d",temp,&i,&j,&r,&c);
                if (r*c==0)
                    break;
                parse_xlat(i,j,r,c);
                break;
            case group_type:
                sscanf(s.c_str(),"%s %d %c %s %d %d",temp,&i,temp,temp,&r,&c);
                parse_v(i,r,c);
                break;
            }
        }
    }
};

Config config;
void parse_args(int argc , char *argv[]){
    std::string::size_type pos;
    string s;
    s = argv[3];
    pos = s.find("n");
    if ( pos != std::string::npos )
        config.n = true;
    pos = s.find("f");
    if ( pos != std::string::npos )
        config.f = true;
    pos = s.find("t");
    if ( pos != std::string::npos )
        config.t = true;
    pos = s.find("s");
    if ( pos != std::string::npos )
        config.s = true;
    pos = s.find("O");
    if ( pos != std::string::npos )
        config.O = true;
    pos = s.find("S");
    if ( pos != std::string::npos )
        config.S = true;
    pos = s.find("a");
    if ( pos != std::string::npos )
        config.a = true;
    pos = s.find("w");
    if ( pos != std::string::npos )
        config.w = true;
    config.n = !config.f;
    config.O = !config.S;
    config.w = !config.a;
    config.s = !config.t;

}
void initialize(){
}
void nbody_solver(){
}
void fmm_solver(){
    Tree &OT=*new Tree;
    char f1[100],f2[100];
    sprintf(f1,"tree_%d_%d.txt",N,L);
    sprintf(f2,"tree_op_%d_%d.txt",N,L);
    Matrix &pts = *new Matrix (N,3);
    Reader rdr(f1,f2,OT);
    rdr.read();
    rdr.read_op();
    OT.Q = 10;
    init(pts);
    compute_near_field(OT,pts);
    Matrix &c = * new Matrix (N,1,1.0);
    Matrix &q = * new Matrix (N,1,0.0);

    SGMatrix &C = *new SGMatrix(c);
    SGMatrix &Q = *new SGMatrix(q);
    mv_near_field(OT,C,Q);
    MatVec(OT,C,Q);
    char res[25];
    sprintf(res,"results_%c_%c_%c_%c.txt",
            config.n?'n':'f',
            config.s?'s':'t',
            config.S?'S':'O',
            config.w?'w':'a'
            );
    Q.get_matrix()->export_data(res);
}
int main(int argc , char *argv[])
{
    if ( argc <2){
        fprintf(stderr,"Usage %s N L nfstSOwa\n",argv[0]);
        exit(-1);
    }
    N=atoi(argv[1]);
    L=atoi(argv[2]);
    parse_args(argc,argv);

    tic();
    fmm_solver();

    cout << " Finished. Time(s): " << toc() << endl;

}


