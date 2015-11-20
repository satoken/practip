#include <iostream>
#include <fstream>
#include <stdexcept>
#include <cerrno>
#include <cstring>
#include "alignment.h"
#include "typedef.h"

template <class T>
Alignment<T>::
Alignment(const char* align, const char* score)
  : q_score_(), idx_()
{
  if (!read_alignment(align))
    throw std::runtime_error(std::string(strerror(errno)) + ": " + align);
  if (!read_scores(score))
    throw std::runtime_error(std::string(strerror(errno)) + ": " + score);
}

template <class T>
bool
Alignment<T>::
read_alignment(const char* filename)
{
  std::ifstream is(filename);
  if (is.fail()) return false;
  std::string line;
  std::string seq;
  while (std::getline(is, line))
  {
    if (line[0]=='>')
    {
      if (!seq.empty())
      {
        idx_.emplace_back(seq.size(), -1U);
        VU& x=idx_.back();
        for (uint i=0, p=0; i!=seq.size(); ++i)
          if (seq[i]!='-') x[i]=p++;
      }
      seq.clear();
    }
    else
    {
      seq+=line;
    }
  }
  if (!seq.empty())
  {
    idx_.emplace_back(seq.size(), -1U);
    VU& x=idx_.back();
    for (uint i=0, p=0; i!=seq.size(); ++i)
      if (seq[i]!='-') x[i]=p++;
  }
  return true;
}

template <class T>
bool
Alignment<T>::
read_scores(const char* filename)
{
  std::ifstream is(filename);
  if (is.fail()) return false;
  uint val;
  while (is >> val)
  {
    q_score_.push_back(val / 100.);
  }
  return true;
}

// instantiation
#include "aa.h"
template class Alignment<AA>;

#include "rna.h"
template class Alignment<RNA>;
