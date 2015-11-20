#ifndef _INC_ALIGNMENT_H__
#define _INC_ALIGNMENT_H__

#include <vector>
#include "typedef.h"

template < class T >
class Alignment
{
public:
Alignment() : q_score_(), idx_() {}
  Alignment(const char* algn, const char* score);

public:
  uint num_sequences() const { return idx_.size(); }
  uint num_columns() const { return q_score_.size(); }
  bool read_alignment(const char* filename);
  bool read_scores(const char* filename);
  const VF& q_score() const { return q_score_; }
  const VVU& idx() const { return idx_; }
  void add_seq(const T& seq) { seq_.push_back(seq); }
  template <class... Args>
    void emplace_add_seq(Args&&... args) { seq_.emplace_back(args...); }
  const T& seq(uint i) const { return seq_[i]; }

private:
  VF q_score_;
  VVU idx_;
  std::vector<T> seq_;
};
  
#endif

// Local Variables:
// mode: C++
// End:
