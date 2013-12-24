#ifndef _TEST_SSVM_SEMI_R_H_
#define _TEST_SSVM_SEMI_R_H_

#include <string>

typedef unsigned int uint;
typedef std::vector<float> VF;
typedef std::vector<VF> VVF;
typedef std::vector<int> VI;
typedef std::vector<VI> VVI;
typedef std::vector<uint> VU;
typedef std::vector<VU> VVU;

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
    : feature_weight_(),
      use_feature_(FG_NUM, true)
  { }
  
private:
  std::vector<std::map<std::string,float> > feature_weight_;
  std::vector<bool> use_feature_;

  float positive_penalty_;
  float negative_penalty_;
  double f_clip_;          // fobos clip
  double eta_;         // weight bonus
  uint n_th_;

  float label_learn(const AA& aa, const RNA& rna, const VVU& correct_edges);
  void read_correct_matching(const std::string& filename, VVU& correct_edges) const; 
  template < class Func > void extract_feature(const AA& aa, const RNA& rna, uint i, uint j, Func& func) const;
  void calculate_edge_weight(const AA& aa, const RNA& rna, VVF& edge_weight) const;
  void penalize_correct_matching(VVF& edge_weight, const VVU& correct_edges) const;
  void update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_edges, const VVU& correct_edges);
  void count_feature(const AA& aa, const RNA& rna, const VVU& edge, std::vector<std::map<std::string,uint> >& fc, VU& tc) const;
  void regularization_fobos();
  float predict_matching(const VVF& edge_weight, VVU& p) const;
};

#endif
