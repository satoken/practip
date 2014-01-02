#ifndef _INC_PRACTIP_H_
#define _INC_PRACTIP_H_

#include <vector>
#include <string>

#define USE_SPARSEHASH
#ifdef USE_SPARSEHASH
#include <google/sparse_hash_map>
#else
#include <map>
#endif

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;
#ifdef USE_SPARSEHASH
typedef google::sparse_hash_map<std::string,float> FM;
typedef google::sparse_hash_map<std::string,uint> FC;
#else
typedef std::map<std::string,float> FM;
typedef std::map<std::string,uint> FC;
#endif

class PRactIP
{
public:
  struct RNA {
    std::string seq;
    std::string ss;

    int read(const std::string& filename);
    static char group2(char);
    static void structural_profile(const std::string& ss, std::string& profile);
  };

  struct AA {
    std::string seq;
    std::string ss;

    int read(const std::string& filename);

    static char group10(char a);
    static char group8(char a);
    static char group4(char a);
    static char group2(char a);
  };

  // Feature Groups
  enum {
    FG_Rss5 = 0,
    FG_Pss5,
    FG_PsRs,
    FG_Ps3Rs3,
    FG_AAp,
    FG_AAab,
    FG_AAh,
    FG_Rpp,
    FG_P1,
    FG_P1R1,
    FG_P3R1,
    FG_P3R3,
    FG_P5R1,
    FG_P1R5,
    FG_P5R5,
    FG_NUM
  };
public:
  PRactIP()
    : feature_weight_(FG_NUM),
      feature_count_(FG_NUM),
      feature_group_weight_(FG_NUM, 0.0),
      feature_group_count_(FG_NUM, 0),
      use_feature_(FG_NUM, true),
      n_th_(1)
  { }
  
public:
  PRactIP& parse_options(int& argc, char**& argv);
  int run();

private:  
  uint load_labeled_data(const std::string& filename);
  uint load_unlabeled_data(const std::string& filename);
  void supervised_training();
  void supervised_training(const VU& use_idx);
  float supervised_training(const AA& aa, const RNA& rna, const VVU& correct_int, float eta);
  float unsupervised_training(const AA& aa, const RNA& rna,
                              std::vector<FC>& fc, VU& tc) const;
  void semisupervised_training();
  void semisupervised_training(const VU& use_idx);
  void cross_validation(uint n, void (PRactIP::*train)(const VU&));
  void supervised_cross_validation(uint n);
  void semisupervised_cross_validation(uint n);
  float predict_interaction(const AA& aa, const RNA& rna, VVU& predicted_int) const;

  void read_correct_interaction(const std::string& filename, VVU& correct_int) const; 
  template < class Func > void extract_int_feature(const AA& aa, const RNA& rna, uint i, uint j, Func& func) const;
  template < class Func > void extract_aa_feature(const AA& aa, uint i, Func& func) const;
  template < class Func > void extract_rna_feature(const RNA& rna, uint j, Func& func) const;
  void calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight) const;
  void penalize_correct_interaction(VVF& int_weight, VF& aa_weight, VF& rna_weight, const VVU& correct_int) const;
  void update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_int, const VVU& correct_int, float eta);
  void count_feature(const AA& aa, const RNA& rna, const VVU& predicted_int, std::vector<FC>& fc, VU& tc) const;
  float regularization_fobos(float eta);
  float predict_interaction(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p) const;
  float calculate_score(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, const VVU& interactions) const;

private:
  std::vector<AA> labeled_aa_;
  std::vector<RNA> labeled_rna_;
  std::vector<VVU> labeled_int_;
  std::vector<AA> unlabeled_aa_;
  std::vector<RNA> unlabeled_rna_;
  
  std::vector<FM> feature_weight_;
  std::vector<FC> feature_count_;
  VF feature_group_weight_;
  VU feature_group_count_;
  std::vector<bool> use_feature_;

  float pos_w_;
  float neg_w_;
  float lambda_;
  float mu_;
  float eta0_;
  uint n_th_;
  uint d_max_;
  uint g_max_;
  uint cv_fold_;

};

#endif
