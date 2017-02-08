#include "util.hpp"
#include <cstring>
#include <fstream>

namespace FMM{
  double submit_time;
  Reader::Reader(string s,string s2,Tree &T){
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
  int Reader::get_type(string s){
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
  void Reader::parse(string s,int type)
  {
    string v;
    Tree &T=*t;
    std::string::size_type p1=0,p2;
    while ( true ) {
      p2 = s.find(",",p1);
      if ( p2 == string::npos)
	return;
      v =  s.substr(p1,p2-p1) ;
      /*cout << v << " , " ;*/
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
  void Reader::parse_group_info(string line)
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
  void Reader::parse_v(int g,int r, int c){
    Tree &T=*t;
    Matrix *M=parse_mat(r,c);
    assert(M);
    T.Levels(level).group(g).V.setMatrix(*M);
  }
  double Reader::get_double(string &s){
    string::size_type p1=0,p2;
    p2 = s.find(",",p1);
    double d=atof(s.substr(p1,p2-p1).c_str());
    s=s.substr(p2+1);
    return d;
  }
  Matrix *Reader::parse_mat(int r , int c)
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
    assert(&M);
    return &M;
  }
  void Reader::parse_m2m(int i,int j,int r,int c)
  {
    Tree &T=*t;
    Matrix *M=parse_mat(r,c);
    T.Levels(level).M2M(i,j).setMatrix(*M);

  }
  void Reader::parse_xlat(int i,int j,int r,int c){
    Tree &T=*t;
    Matrix *M=parse_mat(r,c);
    T.Levels(level).T(i,j).setMatrix(*M);
  }
  void Reader::read(){
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
	//       cout << "Level: " << level <<" Group: " << group << endl;
      }
    // cout << "-----------------------\n";
    //t->print_lists();
  }
  void Reader::read_op(){
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


  void parse_args(int argc , char *argv[]){
    std::string::size_type pos;
    config.omp = false;
    config.N = N = atoi(argv[1]);
    config.L = L = atoi(argv[2]);
    config.NF=config.FF=false;
    config.Q     = atoi(argv[3]);
    config.P     = atoi(argv[4]);
    config.cores = atoi (argv[5]);
    strcpy(config.tree,argv[6]);
    strcpy(config.ops ,argv[7]);

    string s;
    config.l = config.m = config.x = false;
    s = argv[8];
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
    pos = s.find("l");
    if ( pos != std::string::npos )
      config.l = true;
    pos = s.find("m");
    if ( pos != std::string::npos )
      config.m = true;
    pos = s.find("x");
    if ( pos != std::string::npos )
      config.x = true;
    pos = s.find("h");
    if ( pos != std::string::npos )
      config.h = true;
    pos = s.find("F");
    if ( pos != std::string::npos )
      config.FF = true;
    pos = s.find("N");
    if ( pos != std::string::npos )
      config.NF = true;
    pos = s.find("M");
    if ( pos != std::string::npos )
      config.omp = true;
    config.n = !config.f;
    config.O = !config.S;
    config.w = !config.a;
    config.s = !config.t;
    fprintf(stdout,"config pars:  Q:%d, N: %d , P: %d \n",config.Q, config.N, config.P);
  }
}//namespace FMM
