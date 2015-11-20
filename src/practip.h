#ifndef _INC_PRACTIP_H_
#define _INC_PRACTIP_H_

#include <vector>
#include <string>
#include <unordered_map>
#include "aa.h"
#include "rna.h"
#include "alignment.h"
#include "typedef.h"

class PRactIP
{
public:
  struct FeatureWeight {
    float weight;
    float sum_of_grad2;
    uint last_updated;

    FeatureWeight() : weight(0.0), sum_of_grad2(0.0), last_updated(0) {}
  };

  typedef std::unordered_map<std::string,FeatureWeight> FM;

  // Feature Groups
  enum {
    FG_P_3 = 0,
    FG_P_5,
    FG_R_3,
    FG_R_5,
    FG_P_3_R_3,
    FG_P_5_R_5,
    FG_Pss_3,
    FG_Pss_5,
    FG_Rss_3,
    FG_Rss_5,
    FG_Pss_3_Rss_3,
    FG_Pss_5_Rss_5,
    FG_Pg10_5,
    FG_Pg10_7,
    FG_Pg10_3_R_3,
    FG_Pg10_3_Rss_3,
    FG_Pg10_5_R_5,
    FG_Pg10_5_Rss_5,
    FG_Pg4_5,
    FG_Pg4_7,
    FG_Pg4_3_R_3,
    FG_Pg4_3_Rss_3,
    FG_Pg4_5_R_5,
    FG_Pg4_5_Rss_5,
    FG_NUM
  };

public:
  PRactIP()
    : feature_weight_(FG_NUM),
      use_feature_(FG_NUM, true),
      train_mode_(false),
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
  void store_parameters(const char* filename) const;
  void restore_parameters(const char* filename);
  void default_parameters();

  void supervised_training() { semisupervised_training(); }
  void supervised_training(const VU& use_idx) { semisupervised_training(use_idx); }
  float supervised_training(const AA& aa, const RNA& rna, const VVU& correct_int, bool max_margin=true, float w=1.0);
  void semisupervised_training();
  void semisupervised_training(const VU& use_idx);

  void cross_validation(uint n);


  template < class Func > void extract_int_feature(const AA& aa, const RNA& rna, uint i, uint j, Func func) const;
  template < class Func > void extract_aa_feature(const AA& aa, uint i, Func func) const;
  template < class Func > void extract_rna_feature(const RNA& rna, uint j, Func func) const;

  void read_correct_interaction(const std::string& filename, VVU& correct_int) const; 
  void calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight);
  void penalize_correct_interaction(VVF& int_weight, VF& aa_weight, VF& rna_weight, const VVU& correct_int) const;
  void update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_int, const VVU& correct_int, float w=1.0);

  float update_fobos(uint fgroup, const char* fname);
  float regularization_fobos();

  float predict_interaction(const AA& aa, const RNA& rna, VVU& predicted_int, float w=1.0);
  float predict_interaction(const AA& aa, const RNA& rna, const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p, float w=1.0) const;
  void predict_interaction_object(const AA& aa, const RNA& rna,
                                  const VVF& int_weight, const VF& aa_weight, const VF& rna_weight,
                                  VI& x, VI& y, VVI& z, VI& sl_x, VI& sl_y, IP& ip, float w=1.0) const;
  void predict_interaction_constraints(const AA& aa, const RNA& rna,
                                       const VI& x, const VI& y, const VVI& z, const VI& sl_x, const VI& sl_y, 
                                       IP& ip) const;

  float predict_common_interaction(const Alignment<AA>& aa, const Alignment<RNA>& rna, VVVU& predicted_int);
  float predict_common_interaction(const Alignment<AA>& aa, const Alignment<RNA>& rna, 
                                   const VVVF& int_weight, const VVF& aa_weight, const VVF& rna_weight, 
                                   VVVU& predicted_int) const;

  float calculate_score(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, const VVU& interactions) const;

  void show_result(const AA& aa, const RNA& rna, const VVU& predicted_int, float score) const;

private:
  std::vector<AA> labeled_aa_;
  std::vector<RNA> labeled_rna_;
  std::vector<VVU> labeled_int_;
  std::vector<Alignment<AA>> unlabeled_aa_;
  std::vector<Alignment<RNA>> unlabeled_rna_;
  
  std::vector<FM> feature_weight_;
  std::vector<bool> use_feature_;

  bool train_mode_;
  std::string param_file_;
  float pos_w_;
  float neg_w_;
  float lambda_;
  float mu_;
  float nu_;
  float eta0_;
  float exceed_penalty_;
  uint n_th_;
  uint d_max_;
  uint g_max_;
  uint cv_fold_;
  uint aa_int_max_;
  uint rna_int_max_;
  std::vector<std::string> args_;

private:
  static uint epoch;
};

#endif

// Local Variables:
// mode: C++
// End:
