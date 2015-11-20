#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cassert>
#include "aa.h"
#include "typedef.h"

int
AA::
read(const std::string& fa_name, const std::string& ss_name)
{
  read_fa(fa_name);
  read_ss(ss_name);

  g10.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g10.begin(), group10);
  g8.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g8.begin(), group8);
  g4.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g4.begin(), group4);
  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(seq.size()==ss.size());
  return seq.size();
}

int
AA::
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
  ifs.close();

  return this->seq.size();
}
  

int
AA::
read_ss(const std::string& ss_name)
{
  std::string line;
  std::ifstream ifs(ss_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + ss_name);

  std::getline(ifs, line);
  
  if (line[0]=='#') {           // suppose PSIPRED format
    while (std::getline(ifs, line)) {
      if (line[0]=='\n' || line[0]=='#') continue;
      switch (line[7]) {
        case 'E': case 'H':
          ss.push_back(line[7]);
          break;
        default:
          //assert(!"unreachable"); break;
        case 'C':         
          ss.push_back('C');
          break;
      }
    }
  } else {                      // suppose VIENNA-like format
    std::getline(ifs, line);    // this is the sequence.
    std::getline(ifs, line);    // this is the secondary structure
    for (uint i=0; i!=line.size(); ++i) {
      switch (line[i]) {
        case 'E': case 'H':
          ss.push_back(line[i]);
          break;
        default:
          //assert(!"unreachable"); break;
        case 'C':         
          ss.push_back('C');
          break;
      }
    }    
  }

  return ss.size();
}

// static
//inline
uint
AA::
max_intraction(char x)
{
  return 3;
}

// groups defined by Murphy et al., Protein Eng., 2000, 13(3), pp.149-152
// static
char
AA::
group10(char a)
{
  switch (a) {
    case 'L': case 'V': case 'I': case 'M':
      return 'I'; break; // the smallest alphabet is used as the representative
    case 'C':
      return 'C'; break;
    case 'A':
      return 'A'; break;
    case 'G':
      return 'G'; break;
    case 'S': case 'T':
      return 'S'; break;    
    case 'P':
      return 'P'; break;
    case 'F': case 'Y': case 'W':
      return 'F'; break;
    case 'E': case 'D': case 'N': case 'Q':
      return 'D'; break;
    case 'K': case 'R':
      return 'K'; break;
    case 'H':
      return 'H'; break;
  }
  return a;
}

// static
char
AA::
group8(char a)
{
  switch (a) {
    case 'L': case 'V': case 'I': case 'M': case 'C':
      return 'C'; break;
    case 'A': case 'G':
      return 'A'; break;
    case 'S': case 'T':
      return 'S'; break;    
    case 'P':
      return 'P'; break;
    case 'F': case 'Y': case 'W':
      return 'F'; break;
    case 'E': case 'D': case 'N': case 'Q':
      return 'D'; break;
    case 'K': case 'R':
      return 'K'; break;
    case 'H':
      return 'H'; break;
  }
  return a;
}

// static
char
AA::
group4(char a)
{
  switch (a) {
    case 'L': case 'V': case 'I': case 'M': case 'C':
      return 'C'; break;
    case 'A': case 'G': case 'S': case 'T': case 'P':
      return 'A'; break;
    case 'F': case 'Y': case 'W':
      return 'F'; break;
    case 'E': case 'D': case 'N': case 'Q': case 'K': case 'R': case 'H':
      return 'E'; break;
  }
  return a;
}

// static
char
AA::
group2(char a)
{
  switch (a) {
    case 'L': case 'V': case 'I': case 'M': case 'C':
    case 'A': case 'G': case 'S': case 'T': case 'P':
    case 'F': case 'Y': case 'W':
      return 'A'; break;
    case 'E': case 'D': case 'N': case 'Q': case 'K': case 'R': case 'H':
      return 'D'; break;
  }
  return a;
}

