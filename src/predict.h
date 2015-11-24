#ifndef _INC_PREDICT_H_
#define _INC_PREDICT_H_

#include "feature.h"
#include "typedef.h"

class Predictor
{
public:
  Predictor(float exceed_penalty, uint aa_int_max, uint rna_int_max, float mu, float nu, uint n_th)
    : exceed_penalty_(exceed_penalty), 
      aa_int_max_(aa_int_max),
      rna_int_max_(rna_int_max),
      n_th_(n_th)
  {
  }

  float predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna, VVU& predicted_int) const;
  float predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna, 
                            const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p, float mu) const;
  void predict_interaction_object(const AA& aa, const RNA& rna,
                                  const VVF& int_weight, const VF& aa_weight, const VF& rna_weight,
                                  VI& x, VI& y, VVI& z, VI& sl_x, VI& sl_y, IP& ip, float mu) const;
  void predict_interaction_constraints(const AA& aa, const RNA& rna,
                                       const VI& x, const VI& y, const VVI& z, const VI& sl_x, const VI& sl_y, 
                                       IP& ip) const;

  float predict_common_interaction(const FeatureManager& fm, const Alignment<AA>& aa, const Alignment<RNA>& rna, VVVU& predicted_int);
  float predict_common_interaction(const FeatureManager& fm, const Alignment<AA>& aa, const Alignment<RNA>& rna, 
                                   const VVVF& int_weight, const VVF& aa_weight, const VVF& rna_weight, 
                                   VVVU& predicted_int) const;

  float calculate_score(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, const VVU& interactions) const;

private:
  float exceed_penalty_;
  uint aa_int_max_;
  uint rna_int_max_;
  float mu_;
  float nu_;
  uint n_th_;
};

#endif

// Local Variables:
// mode: C++
// End:
