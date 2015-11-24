#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <stack>
#include <cassert>
#include <cerrno>
#include <cstring>
#include "rna.h"
#include "typedef.h"

// static
//inline
uint
RNA::
max_intraction(char x)
{
  switch (x)
  {
    case 'E': return 7; break;
    case 'H': return 5; break;
  }
  return 4;
}

int
RNA::
read(const std::string& fa_name, const std::string& ss_name)
{
  read_fa(fa_name);
  read_ss(ss_name);

  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(this->seq.size()==this->ss.size());
  return this->seq.size();
}

int
RNA::
read_fa(const std::string& fa_name) 
{
  std::string line;
  std::ifstream ifs(fa_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + fa_name);
  while (std::getline(ifs, line)) {
    if (line[0]=='>')
      this->name+=line.substr(1);
    else
      this->seq+=line;
  }

  return this->seq.size();
}

int
RNA::
read_ss(const std::string& ss_name) 
{
  std::string line;
  std::ifstream ifs(ss_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + ss_name);

#if 0
  FILE* fp = fopen((basename+".ss").c_str(), "r");
  if (fp==nullptr) {
    const char* prog=getenv("CENTROID_FOLD");
    if (!prog) prog="centroid_fold";
    char cmd[1000];
    sprintf(cmd, "%s %s", prog, filename.c_str());
    fp = popen(cmd, "r");
  }

  std::string ss_str;
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]=='>' || line[0]=='f') continue;
    if (strchr(".()", line[0])) {
      for (uint i=0; line[i]!='\0'; ++i) {
        switch (line[i]) {
          case '.':
          case '(':
          case ')':
            break;
          case '\n':
          case ' ':
            line[i]='\0';
            break;
        }
      }
      ss_str+=line;
    }
  }
  fclose(fp);
#endif

  std::getline(ifs, line);      // this is the header;
  std::getline(ifs, line);      // this is the sequence.
  ifs >> line;                  // this is the secondary structure
  structural_profile(line, ss);

  return ss.size();
}

//static
char
RNA::
group2(char b)
{
  switch (b) {
    case 'A': case 'G':
      return 'R'; break;
    case 'C': case 'T': case 'U':
      return 'Y'; break;
  }
  return b;
}

// static
void
RNA::
structural_profile(const std::string& ss, std::string& profile)
{
  VU p(ss.size(), -1u);
  std::string ss2(ss);
  std::stack<uint> st;
  for (uint i=0; i!=ss2.size(); ++i) {
    switch (ss2[i]) {
      case '(':
        st.push(i);
        break;
      case ')':
        p[st.top()]=i;
        p[i]=st.top();
        st.pop();
        if (ss2[i+1]!=')' || st.top()+1!=p[i]) {
          ss2[p[i]]='[';
          ss2[i]=']';
        }
        break;
      default:
        break;
    }
  }

  profile.resize(ss.size());
  std::fill(profile.begin(), profile.end(), 'E');
  std::stack<uint> loop_degree;
  std::stack<bool> bulge;
  for (uint i=0; i!=ss2.size(); ++i) {
    switch (ss2[i]) {
      case '(':
        profile[i]='S';
        break;
      case ')':
        profile[i]='S';
        if (ss2[i-1]==']') bulge.top()=true;
        break;
      case '[':
        profile[i]='S';
        if (i>0 && (ss2[i-1]=='(' || ss2[i-1]=='[')) bulge.top()=true;
        loop_degree.push(1);
        bulge.push(false);
        break;
      case ']':
        char ps;
        profile[i]='S';
        if (ss2[i-1]==']') bulge.top()=true;
        switch (loop_degree.top()) {
          case 1: ps='H'; break;
          case 2: ps=bulge.top() ? 'B' : 'I'; break;
          default: ps='M'; break;
        }
        loop_degree.pop();
        if (!loop_degree.empty()) loop_degree.top()++;
        bulge.pop();
        for (uint j=p[i]+1; j!=i; ++j)
          if (profile[j]=='E') profile[j]=ps;
        break;
    }
  }
}
