#ifndef _TEST_SSVM_SEMI_R_H_
#define _TEST_SSVM_SEMI_R_H_

#include <vector>
#include <string>

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;
typedef std::map<std::string,float> FM;
typedef std::map<std::string,uint> FC;

class PRactIP
{
public:
  struct RNA {
    std::string seq;
    std::string ss;

    int read(const std::string& filename);
    static int pp2(char);
  };

  struct AA {
    std::string seq;
    std::string ss;

    int read(const std::string& filename);
    static int group7(char);
    static int ab3(char);
    static int hydro3(char);
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
    : labeled_aa_(), labeled_rna_(), labeled_matching_(),
      unlabeled_aa_(), unlabeled_rna_(), 
      feature_weight_(FG_NUM),
      feature_count_(FG_NUM),
      feature_group_weight_(FG_NUM, 0.0),
      feature_group_count_(FG_NUM, 0),
      use_feature_(FG_NUM, true),
      positive_penalty_(1.0),
      negative_penalty_(1.0),
      lambda_(1.0),
      mu_(1.0),
      eta0_(0.5),
      n_th_(1),
      t_max_(1000),
      u_max_(1000)
  { }
  
private:
  std::vector<AA> labeled_aa_;
  std::vector<RNA> labeled_rna_;
  std::vector<VVU> labeled_matching_;
  std::vector<AA> unlabeled_aa_;
  std::vector<RNA> unlabeled_rna_;
  
  std::vector<FM> feature_weight_;
  std::vector<FC> feature_count_;
  VF feature_group_weight_;
  VU feature_group_count_;
  std::vector<bool> use_feature_;

  float positive_penalty_;
  float negative_penalty_;
  float lambda_;
  float mu_;
  float eta0_;
  uint n_th_;
  uint t_max_;
  uint u_max_;

  uint load_labeled_data(const std::string& filename);
  uint load_unlabeled_data(const std::string& filename);
  void supervised_training();
  void supervised_training(const VU& use_idx);
  float supervised_training(const AA& aa, const RNA& rna, const VVU& correct_edges, float eta);
  float unsupervised_training(const AA& aa, const RNA& rna,
                              std::vector<FC>& fc, VU& tc) const;
  void semisupervised_training();
  void semisupervised_training(const VU& use_idx);
  void cross_validation(uint n, void (PRactIP::*train)(const VU&));
  void supervised_cross_validation(uint n);
  void semisupervised_cross_validation(uint n);
  float predict_matching(const AA& aa, const RNA& rna, VVU& predicted_edges) const;

  void read_correct_matching(const std::string& filename, VVU& correct_edges) const; 
  template < class Func > void extract_feature(const AA& aa, const RNA& rna, uint i, uint j, Func& func) const;
  void calculate_edge_weight(const AA& aa, const RNA& rna, VVF& edge_weight) const;
  void penalize_correct_matching(VVF& edge_weight, const VVU& correct_edges) const;
  void update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_edges, const VVU& correct_edges, float eta);
  void count_feature(const AA& aa, const RNA& rna, const VVU& edge, std::vector<FC>& fc, VU& tc) const;
  float regularization_fobos(float eta);
  float predict_matching(const VVF& edge_weight, VVU& p) const;
};

#endif
