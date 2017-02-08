#ifndef UTIL_HPP
#define UTIL_HPP

#include <string>
#include <fstream>
#include "types.hpp"

extern int N,L;
extern double submit_time;
using namespace std;

void parse_args(int argc , char *argv[]);


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
  Reader(string s,string s2,Tree &T);
  int get_type(string s);
  void parse(string s,int type);
  void parse_group_info(string line);
  void parse_v(int g,int r, int c);
  double get_double(string &s);
  Matrix *parse_mat(int r , int c);
  void parse_m2m(int i,int j,int r,int c);
  void parse_xlat(int i,int j,int r,int c);
  void read();
  void read_op();
};


#endif// UTIL_HPP
