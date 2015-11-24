#ifndef _INC_PRACTIP_H_
#define _INC_PRACTIP_H_

#include <vector>
#include <string>
#include <unordered_map>
#include "aa.h"
#include "rna.h"
#include "alignment.h"
#include "feature.h"
#include "typedef.h"

class PRactIP
{
public:
  PRactIP()
    : train_mode_(false),
      n_th_(1),
      aa_int_max_(-1u),
      rna_int_max_(-1u)
  { }
  
public:
  PRactIP& parse_options(int& argc, char**& argv);
  int run();

private:  
  uint load_labeled_data(const std::string& filename);
  uint load_unlabeled_data(const std::string& filename);

  void supervised_training(FeatureManager& fm) { semisupervised_training(fm); }
  void supervised_training(FeatureManager& fm, const VU& use_idx) { semisupervised_training(fm, use_idx); }
  float supervised_training(FeatureManager& fm, const AA& aa, const RNA& rna, const VVU& correct_int, 
                            bool max_margin=true, float w=1.0);
  void semisupervised_training(FeatureManager& fm);
  void semisupervised_training(FeatureManager& fm, const VU& use_idx);

  void cross_validation(uint n);
  void read_correct_interaction(const std::string& filename, VVU& correct_int) const; 

  void penalize_correct_interaction(VVF& int_weight, VF& aa_weight, VF& rna_weight, const VVU& correct_int) const;

  void show_result(const AA& aa, const RNA& rna, const VVU& predicted_int, float score) const;

private:
  std::vector<AA> labeled_aa_;
  std::vector<RNA> labeled_rna_;
  std::vector<VVU> labeled_int_;
  std::vector<Alignment<AA>> unlabeled_aa_;
  std::vector<Alignment<RNA>> unlabeled_rna_;
  
  bool train_mode_;
  std::string out_file_;
  std::string param_file_;
  float pos_w_;
  float neg_w_;
  float lambda_;
  float mu_;
  float nu_;
  float eta0_;
  uint d_max_;
  uint g_max_;
  uint cv_fold_;

  float exceed_penalty_;
  uint n_th_;
  uint aa_int_max_;
  uint rna_int_max_;
  std::vector<std::string> args_;
};

#endif

// Local Variables:
// mode: C++
// End:
