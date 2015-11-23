#ifndef _INC_FEATURE_H_
#define _INC_FEATURE_H_

#include <vector>
#include <unordered_map>
#include <string>
#include "aa.h"
#include "rna.h"
#include "typedef.h"

class FeatureManager
{
private:
  struct FeatureWeight 
  {
    float weight;
    float sum_squared_grad;
    uint last_updated;

    FeatureWeight() : weight(0.0), sum_squared_grad(0.0), last_updated(0) {}
  };

public:
  FeatureManager(float lambda, float eta0);

  void store_parameters(const char* filename) const;
  void restore_parameters(const char* filename);
  void default_parameters();

  void calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight) const; 
  void update_feature_weight(const AA& aa, const RNA& rna, 
                             const VVU& predicted_int, const VVU& correct_int, float w=1.0);

  float regularization() const;

private:
  float lazy_update(FeatureWeight& f) const;
  float lazy_update(uint fgroup, const char* fname) const;

  template < class Func > 
  void extract_int_feature(const AA& aa, const RNA& rna, uint i, uint j, Func func) const;
  template < class Func > 
  void extract_aa_feature(const AA& aa, uint i, Func func) const;
  template < class Func > 
  void extract_rna_feature(const RNA& rna, uint j, Func func) const;

private:
  typedef std::unordered_map<std::string,FeatureWeight> FM;

  mutable std::vector<FM> feature_weight_; // mutable for lazy updates
  std::vector<bool> use_feature_;
  uint epoch_;
  float lambda_;
  float eta0_;
};

#endif

// Local Variables:
// mode: C++
// End:
