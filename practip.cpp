#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stack>
#include <functional>
#include <getopt.h>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cassert>


#include "ip.h"
#include "practip.h"
#include "cmdline.h"

#define FOREACH(itr, i, v) for (itr i=(v).begin(); i!=(v).end(); ++i)

struct Accuracy
{
  uint tp, tn, fp, fn;
  float ppv, sen, fval;
  void calculate(const std::vector<bool>& c, const std::vector<bool>& p)
  {
    tp=tn=fp=fn=0;
    assert(c.size()==p.size());
    for (uint i=0; i!=c.size(); ++i)
    {
      if (c[i] && p[i]) tp++;
      else if (!c[i] && !p[i]) tn++;
      else if (!c[i] && p[i]) fp++;
      else /*if (c[i] && !p[i])*/ fn++;
    }
    ppv = tp+fp!=0 ? (float)tp/(tp+fp) : 0;
    sen = tp+fn!=0 ? (float)tp/(tp+fn) : 0;
    fval = ppv+sen!=0 ? 2*ppv*sen/(ppv+sen) : 0;
  }
};

struct AccuracySummary
{
  uint tp, tn, fp, fn;
  float ppv, sen, fval;
  float ppv2, sen2, fval2;
  float ppv_avg, sen_avg, fval_avg;
  float ppv_sd, sen_sd, fval_sd;
  uint n;

  AccuracySummary()
    : tp(0), tn(0), fp(0), fn(0),
      ppv(0.0), sen(0.0), fval(0.0),
      ppv2(0.0), sen2(0.0), fval2(0.0),
      n(0)
  { }

  void add(const Accuracy& acc)
  {
    tp += acc.tp;
    tn += acc.tn;
    fp += acc.fp;
    fn += acc.fn;
    ppv += acc.ppv;
    ppv2 += acc.ppv*acc.ppv;
    sen += acc.sen;
    sen2 += acc.sen*acc.sen;
    fval2 += acc.fval*acc.fval;
    n++;
  }

  float avg(float sum) const
  {
    return sum/n;
  }

  float sd(float sum, float sum2) const
  {
    float a=avg(sum);
    float b=avg(sum2);
    return b-a*a>=0.0 ? std::sqrt(b-a*a) : 0.0;
  }

  void summary()
  {
    ppv_avg = avg(ppv); ppv_sd = sd(ppv, ppv2);
    sen_avg = avg(sen); sen_sd = sd(sen, sen2);
    fval_avg = avg(fval); fval_sd = sd(fval, fval2);
  }
};

static
void
calculate_accuracy(const PRactIP::AA& aa, const PRactIP::RNA& rna,
                   const VVU& predicted_int, const VVU& correct_int,
                   Accuracy& int_acc, Accuracy& aa_acc, Accuracy& rna_acc)
{
  std::vector<bool> c_int(aa.seq.size()*rna.seq.size(), false);
  std::vector<bool> p_int(aa.seq.size()*rna.seq.size(), false);
  std::vector<bool> c_aa(aa.seq.size(), false);
  std::vector<bool> p_aa(aa.seq.size(), false);
  std::vector<bool> c_rna(rna.seq.size(), false);
  std::vector<bool> p_rna(rna.seq.size(), false);
  for (uint i=0; i!=correct_int.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_int[i])
      c_int[rna.seq.size()*i+*j]=c_aa[i]=c_rna[*j]=true;
  for (uint i=0; i!=predicted_int.size(); ++i)
    FOREACH (VU::const_iterator, j, predicted_int[i])
      p_int[rna.seq.size()*i+*j]=p_aa[i]=p_rna[*j]=true;
  int_acc.calculate(c_int, p_int);
  aa_acc.calculate(c_aa, p_aa);
  rna_acc.calculate(c_rna, p_rna);
}

uint
PRactIP::
load_labeled_data(const std::string& filename)
{
  std::string aa_seq, rna_seq, matching;
  std::ifstream is(filename.c_str());
  std::cout << "loading labeled data" << std::endl;
  while (is >> aa_seq >> rna_seq >> matching) {
    std::cout << aa_seq << " " << rna_seq << " " << matching << std::endl;
    labeled_aa_.push_back(AA());
    uint aa_len=labeled_aa_.back().read(aa_seq);
    labeled_rna_.push_back(RNA());
    labeled_rna_.back().read(rna_seq);
    labeled_int_.push_back(VVU(aa_len));
    read_correct_interaction(matching, labeled_int_.back());
  }
  return labeled_aa_.size();
}

uint
PRactIP::
load_unlabeled_data(const std::string& filename)
{
  std::string aa_seq, rna_seq;
  std::ifstream is(filename.c_str());
  std::cout << "loading unlabeled data" << std::endl;
  while (is >> aa_seq >> rna_seq) {
    std::cout << aa_seq << " " << rna_seq << std::endl;
    unlabeled_aa_.push_back(AA());
    unlabeled_aa_.back().read(aa_seq);
    unlabeled_rna_.push_back(RNA());
    unlabeled_rna_.back().read(rna_seq);
  }
  return unlabeled_aa_.size();
}

void
PRactIP::
supervised_training()
{
  // use all labeled data
  VU idx(labeled_aa_.size());
  for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
  supervised_training(idx);
}

void
PRactIP::
supervised_training(const VU& use_idx)
{
  VU idx(use_idx);
  for (uint t=0; t!=d_max_; ++t) {
    float eta=eta0_/std::sqrt(t+1); // update learning rate
    float total_loss=0.0;
    std::random_shuffle(idx.begin(), idx.end()); // shuffle the order of training data
    FOREACH (VU::const_iterator, it, idx) {
      total_loss += supervised_training(labeled_aa_[*it], labeled_rna_[*it], labeled_int_[*it], eta);
    }
  }
}

float
PRactIP::
calculate_score(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, const VVU& interactions) const
{
  float score = 0.0;
  for (uint i=0; i!=interactions.size(); ++i)
    FOREACH (VU::const_iterator, j, interactions[i])
      score += int_weight[i][*j];
  for (uint i=0; i!=interactions.size(); ++i)
    score += aa_weight[i];
  std::vector<bool> rna_has_int(rna_weight.size(), false);
  for (uint i=0; i!=interactions.size(); ++i)
    FOREACH (VU::const_iterator, j, interactions[i])
      rna_has_int[*j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      score += rna_weight[j];
  return score;
}


float
PRactIP::
supervised_training(const AA& aa, const RNA& rna, const VVU& correct_int, float eta)
{
  float loss=0.0;
  VVF int_weight;
  VF aa_weight, rna_weight;
  calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  loss -= calculate_score(int_weight, aa_weight, rna_weight, correct_int);
  penalize_correct_interaction(int_weight, aa_weight, rna_weight, correct_int);
  VVU predicted_int;
  predict_interaction(int_weight, aa_weight, rna_weight, predicted_int);
  loss += calculate_score(int_weight, aa_weight, rna_weight, predicted_int);
#if 0
  for (uint i=0; i!=correct_int.size(); ++i) {
    if (!correct_int[i].empty() || !predicted_int[i].empty())
    {
      std::cout << i << ": [ ";
      FOREACH (VU::const_iterator, j, correct_int[i])
        std::cout << *j << "(" << int_weight[i][*j] << ") ";
      std::cout << "], [ ";
      FOREACH (VU::const_iterator, j, predicted_int[i])
        std::cout << *j << "(" << int_weight[i][*j] << ") ";
      std::cout << "]" << std::endl;
    }
  }
  std::cout << std::endl;
#endif  
  update_feature_weight(aa, rna, predicted_int, correct_int, eta);
  loss += regularization_fobos(eta);
  return loss;
}

void
PRactIP::
cross_validation(uint n, void (PRactIP::*training)(const VU&))
{
  AccuracySummary int_summary, aa_summary, rna_summary;
  for (uint i=0; i!=n; ++i) {
    VU train;
    VU test;
    for (uint j=0; j!=labeled_aa_.size(); ++j) {
      if (j%n==i)
        test.push_back(j);
      else
        train.push_back(j);
    }

    std::cout << "[Set " << i << "]" << std::endl;
    (this->*training)(train);
    
    FOREACH (VU::const_iterator, j, test) {
      VVU predicted_int;
      predict_interaction(labeled_aa_[*j], labeled_rna_[*j], predicted_int);
      Accuracy int_acc, aa_acc, rna_acc;
      calculate_accuracy(labeled_aa_[*j], labeled_rna_[*j],
                         predicted_int, labeled_int_[*j],
                         int_acc, aa_acc, rna_acc);
      int_summary.add(int_acc);
      aa_summary.add(aa_acc);
      rna_summary.add(rna_acc);
    }
  }
  int_summary.summary();
  aa_summary.summary();
  rna_summary.summary();
  std::cout << "Interaction: "
            << "PPV=" << int_summary.ppv_avg << "(" << int_summary.ppv_sd << ") "
            << "SEN=" << int_summary.sen_avg << "(" << int_summary.sen_sd << ") "
            << "F=" << int_summary.fval_avg << "(" << int_summary.fval_sd << ") ["
            << int_summary.tp << "," << int_summary.tn << ","
            << int_summary.fp << "," << int_summary.fn << "]"
            << std::endl;
  std::cout << "AA: "
            << "PPV=" << aa_summary.ppv_avg << "(" << aa_summary.ppv_sd << ") "
            << "SEN=" << aa_summary.sen_avg << "(" << aa_summary.sen_sd << ") "
            << "F=" << aa_summary.fval_avg << "(" << aa_summary.fval_sd << ") ["
            << aa_summary.tp << "," << aa_summary.tn << ","
            << aa_summary.fp << "," << aa_summary.fn << "]"
            << std::endl;
  std::cout << "RNA: "
            << "PPV=" << rna_summary.ppv_avg << "(" << rna_summary.ppv_sd << ") "
            << "SEN=" << rna_summary.sen_avg << "(" << rna_summary.sen_sd << ") "
            << "F=" << rna_summary.fval_avg << "(" << rna_summary.fval_sd << ") ["
            << rna_summary.tp << "," << rna_summary.tn << ","
            << rna_summary.fp << "," << rna_summary.fn << "]"
            << std::endl;
}

void
PRactIP::
supervised_cross_validation(uint n)
{
  void (PRactIP::*train)(const VU&) = &PRactIP::supervised_training;
  cross_validation(n, train);
}

void
PRactIP::
semisupervised_cross_validation(uint n)
{
  void (PRactIP::*train)(const VU&) = &PRactIP::semisupervised_training;
  cross_validation(n, train);
}

void
PRactIP::
semisupervised_training()
{
  VU idx(labeled_aa_.size());
  for (uint i=0; i!=idx.size(); ++i) idx[i]=i;
  semisupervised_training(idx);
}

void
PRactIP::
semisupervised_training(const VU& use_idx)
{
  supervised_training(use_idx);

  VU idx(use_idx);
  for (uint u=0; u!=g_max_; ++u) {
    float total_score=0.0;
    std::vector<FC> fc(FG_NUM);
    VU tc(FG_NUM, 0);
    for (uint i=0; i!=unlabeled_aa_.size(); ++i)
      total_score += unsupervised_training(unlabeled_aa_[i], unlabeled_rna_[i], fc, tc);
    std::swap(feature_count_, fc);
    std::swap(feature_group_count_, tc);

    float eta=eta0_/std::sqrt(u+1);
    float total_loss=0.0;
    std::random_shuffle(idx.begin(), idx.end());
    FOREACH (VU::const_iterator, it, idx) {
      total_loss += supervised_training(labeled_aa_[*it], labeled_rna_[*it], labeled_int_[*it], eta);
    }
  }    
}

float
PRactIP::
unsupervised_training(const AA& aa, const RNA& rna,
                      std::vector<FC>& fc, VU& tc) const
{
  VVU p;
  float score = predict_interaction(aa, rna, p);
  count_feature(aa, rna, p, fc, tc);
  return score;
}

float
PRactIP::
predict_interaction(const AA& aa, const RNA& rna, VVU& predicted_int) const
{
  VVF int_weight;
  VF aa_weight, rna_weight;
  calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  return predict_interaction(int_weight, aa_weight, rna_weight, predicted_int);
}

void
PRactIP::
read_correct_interaction(const std::string& filename, VVU& correct_int) const
{
  uint r, a;
  std::ifstream is(filename.c_str());
  while (is >> a >> r)
    correct_int[a].push_back(r);
}

inline
static
char*
feature_string(const char* str, uint str_len, uint p, uint w, char* fname)
{
  char* f=fname;
  for (uint i=w; i>=1; --i)
    *f++ = p>=i ? str[p-i] : '-';
  *f++ = str[p];
  for (uint i=1; i<=w; ++i)
    *f++ = p+i<str_len ? str[p+i] : '-';
  *f++='\0';
  return fname;
}

template < class Func >
void
PRactIP::
extract_int_feature(const AA& aa, const RNA& rna, uint i, uint j, Func& func) const
{
  struct {
    uint id;
    const char* aa_str;
    uint aa_w;
    const char* rna_str;
    uint rna_w;
  } features[] = {
    { FG_P_3_R_3,      aa.seq.c_str(), 1, rna.seq.c_str(), 1 },
    { FG_P_5_R_5,      aa.seq.c_str(), 2, rna.seq.c_str(), 2 },
    { FG_Pss_3_Rss_3,  aa.ss.c_str(),  1, rna.ss.c_str(),  1 },
    { FG_Pss_5_Rss_5,  aa.ss.c_str(),  2, rna.ss.c_str(),  2 },
    { FG_Pg10_3_R_3,   aa.g10.c_str(), 1, rna.seq.c_str(), 1 },
    { FG_Pg10_3_Rss_3, aa.g10.c_str(), 1, rna.ss.c_str(),  1 },
    { FG_Pg10_5_R_5,   aa.g10.c_str(), 2, rna.seq.c_str(), 2 },
    { FG_Pg10_5_Rss_5, aa.g10.c_str(), 2, rna.ss.c_str(),  2 },
    { FG_Pg4_3_R_3,    aa.g4.c_str(),  1, rna.seq.c_str(), 1 },
    { FG_Pg4_3_Rss_3,  aa.g4.c_str(),  1, rna.ss.c_str(),  1 },
    { FG_Pg4_5_R_5,    aa.g4.c_str(),  2, rna.seq.c_str(), 2 },
    { FG_Pg4_5_Rss_5,  aa.g4.c_str(),  2, rna.ss.c_str(),  2 },
    { -1u, NULL, 0, NULL, 0 }
  };
  
  const uint rna_len=rna.seq.size();
  const uint aa_len=aa.seq.size();
  char buf[20];
  for (uint k=0; features[k].id!=-1u; ++k) {
    if (use_feature_[features[k].id]) {
      feature_string(features[k].aa_str, aa_len, i, features[k].aa_w, buf);
      buf[features[k].aa_w*2+1]=',';
      feature_string(features[k].rna_str, rna_len, j, features[k].rna_w, buf+features[k].aa_w*2+2);
      func(features[k].id, buf, i, j);
      //std::cout << features[k].id << " " << buf << std::endl;
    }
  }
}

template < class Func >
void
PRactIP::
extract_aa_feature(const AA& aa, uint i, Func& func) const
{
  struct {
    uint id;
    const char* aa_str;
    uint aa_w;
  } features[] = {
    { FG_P_3,    aa.seq.c_str(), 1 },
    { FG_P_5,    aa.seq.c_str(), 2 },
    { FG_Pss_3,  aa.ss.c_str(),  1 },
    { FG_Pss_5,  aa.ss.c_str(),  2 },
    { FG_Pg10_5, aa.g10.c_str(), 2 },
    { FG_Pg10_7, aa.g10.c_str(), 3 },
    { FG_Pg4_5,  aa.g4.c_str(),  2 },
    { FG_Pg4_7,  aa.g4.c_str(),  3 },
    { -1u, NULL, 0 }
  };
  
  const uint aa_len=aa.seq.size();
  char buf[20];

  for (uint k=0; features[k].id!=-1u; ++k) {
    if (use_feature_[features[k].id]) {
      feature_string(features[k].aa_str, aa_len, i, features[k].aa_w, buf);
      func(features[k].id, buf, i);
      //std::cout << features[k].id << " " << buf << std::endl;
    }
  }
}

template < class Func >
void
PRactIP::
extract_rna_feature(const RNA& rna, uint j, Func& func) const
{
  struct {
    uint id;
    const char* rna_str;
    uint rna_w;
  } features[] = {
    { FG_R_3,   rna.seq.c_str(), 1 },
    { FG_R_5,   rna.seq.c_str(), 2 },
    { FG_Rss_3, rna.ss.c_str(),  1 },
    { FG_Rss_5, rna.ss.c_str(),  2 },
    { -1u, NULL, 0 }
  };
  
  const uint rna_len=rna.seq.size();
  char buf[20];

  for (uint k=0; features[k].id!=-1u; ++k) {
    if (use_feature_[features[k].id]) {
      feature_string(features[k].rna_str, rna_len, j, features[k].rna_w, buf);
      func(features[k].id, buf, j);
      //std::cout << features[k].id << " " << buf << std::endl;
    }
  }
}

struct EdgeWeightCalculator
{
  EdgeWeightCalculator(const std::vector<FM>& feature_weight,
                       const std::vector<FC>& feature_count,
                       const VF& feature_group_weight,
                       const VU& feature_group_count,
                       VVF& edge_weight)
    : feature_weight_(feature_weight),
      feature_count_(feature_count),
      feature_group_weight_(feature_group_weight),
      feature_group_count_(feature_group_count),
      edge_weight_(edge_weight)
  { }

  inline void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    FM::const_iterator m;
    m=feature_weight_[fgroup].find(fname);
    if (m!=feature_weight_[fgroup].end())
      edge_weight_[i][j] += m->second;

    if (feature_group_count_[fgroup]) {
      FC::const_iterator m;
      m=feature_count_[fgroup].find(fname);
      if (m!=feature_count_[fgroup].end())
        edge_weight_[i][j] += feature_group_weight_[fgroup]*m->second/feature_group_count_[fgroup];
    }
  }

  const std::vector<FM>& feature_weight_;
  const std::vector<FC>& feature_count_;
  const VF& feature_group_weight_;
  const VU& feature_group_count_;
  VVF& edge_weight_;
};

struct NodeWeightCalculator
{
  NodeWeightCalculator(const std::vector<FM>& feature_weight,
                       const std::vector<FC>& feature_count,
                       const VF& feature_group_weight,
                       const VU& feature_group_count,
                       VF& node_weight)
    : feature_weight_(feature_weight),
      feature_count_(feature_count),
      feature_group_weight_(feature_group_weight),
      feature_group_count_(feature_group_count),
      node_weight_(node_weight)
  { }

  inline void operator()(uint fgroup, const char* fname, uint i)
  {
    FM::const_iterator m;
    m=feature_weight_[fgroup].find(fname);
    if (m!=feature_weight_[fgroup].end())
      node_weight_[i] += m->second;

    if (feature_group_count_[fgroup]) {
      FC::const_iterator m;
      m=feature_count_[fgroup].find(fname);
      if (m!=feature_count_[fgroup].end())
        node_weight_[i] += feature_group_weight_[fgroup]*m->second/feature_group_count_[fgroup];
    }
  }

  const std::vector<FM>& feature_weight_;
  const std::vector<FC>& feature_count_;
  const VF& feature_group_weight_;
  const VU& feature_group_count_;
  VF& node_weight_;
};

void
PRactIP::
calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight) const
{
  EdgeWeightCalculator int_weight_calculator(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, int_weight);
  int_weight.resize(aa.seq.size());
  for (uint i=0; i!=int_weight.size(); ++i) {
    int_weight[i].resize(rna.seq.size());
    for (uint j=0; j!=int_weight[i].size(); ++j) {
      int_weight[i][j] = 0.0;
      extract_int_feature(aa, rna, i, j, int_weight_calculator);
    }
  }

  NodeWeightCalculator aa_weight_calculator(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, aa_weight);
  aa_weight.resize(aa.seq.size());
  for (uint i=0; i!=aa_weight.size(); ++i) {
    aa_weight[i] = 0.0;
    extract_aa_feature(aa, i, aa_weight_calculator);
  }

  NodeWeightCalculator rna_weight_calculator(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, rna_weight);
  rna_weight.resize(rna.seq.size());
  for (uint j=0; j!=rna_weight.size(); ++j) {
    rna_weight[j] = 0.0;
    extract_rna_feature(rna, j, rna_weight_calculator);
  }
}

void
PRactIP::
penalize_correct_interaction(VVF& int_weight, VF& aa_weight, VF& rna_weight, const VVU& correct_int) const
{
  assert(int_weight.size()==correct_int.size());
  assert(int_weight.size()==aa_weight.size());

  for (uint i=0; i!=int_weight.size(); ++i) 
    for (uint j=0; j!=int_weight[i].size(); ++j)
      int_weight[i][j] += neg_w_;

  for (uint i=0; i!=correct_int.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_int[i])
      int_weight[i][*j] -= pos_w_+neg_w_;

  for (uint i=0; i!=aa_weight.size(); ++i) 
    aa_weight[i] += neg_w_;

  for (uint i=0; i!=correct_int.size(); ++i)
    if (!correct_int[i].empty())
      aa_weight[i] -= pos_w_+neg_w_;

  for (uint j=0; j!=rna_weight.size(); ++j)
    rna_weight[j] += neg_w_;

  std::vector<bool> rna_has_int(rna_weight.size(), false);
  for (uint i=0; i!=correct_int.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_int[i])
      rna_has_int[*j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      rna_weight[j] -= pos_w_+neg_w_;
}

struct FeatureWeightUpdater
{
  FeatureWeightUpdater(std::vector<FM >& feature_weight,
                       const std::vector<FC >& feature_count,
                       VF& feature_group_weight,
                       const VU& feature_group_count,
                       float eta)
    : feature_weight_(feature_weight),
      feature_count_(feature_count),
      feature_group_weight_(feature_group_weight),
      feature_group_count_(feature_group_count),
      eta_(eta)
  { }

  inline void operator()(uint fgroup, const char* fname, uint i=-1u, uint j=-1u)
  {
    //std::cout << fgroup << ":" << fname << "=" << feature_weight_[fgroup][fname]; 
    feature_weight_[fgroup].insert(std::make_pair(std::string(fname),0.0f)).first->second += eta_;
    //std::cout << " -> " << feature_weight_[fgroup][fname] << std::endl;
    if (feature_group_count_[fgroup]) {
      FC::const_iterator m;
      m=feature_count_[fgroup].find(fname);
      if (m!=feature_count_[fgroup].end())
        feature_group_weight_[fgroup] += eta_*m->second/feature_group_count_[fgroup];
    }
  }

  std::vector<FM >& feature_weight_;
  const std::vector<FC >& feature_count_;
  VF& feature_group_weight_;
  const VU& feature_group_count_;
  float eta_;
};

void
PRactIP::
update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_int, const VVU& correct_int, float eta)
{
  assert(predicted_int.size()==correct_int.size());

  // update feature weights in correct interactions, i.e., add eta
  //   interactions
  FeatureWeightUpdater f(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, eta);
  for (uint i=0; i!=correct_int.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_int[i])
      extract_int_feature(aa, rna, i, *j, f);
  //   amino acids
  for (uint i=0; i!=correct_int.size(); ++i)
    if (!correct_int[i].empty())
      extract_aa_feature(aa, i, f);
  //   RNAs
  std::vector<bool> rna_has_int(correct_int.size(), false);
  for (uint i=0; i!=correct_int.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_int[i])
      rna_has_int[*j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      extract_rna_feature(rna, j, f);

  // update feature weights in predicted interactions, i.e., add -eta
  //   interactions
  FeatureWeightUpdater g(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, -eta);
  for (uint i=0; i!=predicted_int.size(); ++i) 
    FOREACH (VU::const_iterator, j, predicted_int[i]) 
      extract_int_feature(aa, rna, i, *j, g);
  //   amino acids
  for (uint i=0; i!=predicted_int.size(); ++i)
    if (!predicted_int[i].empty())
      extract_aa_feature(aa, i, g);
  //   RNAs
  std::fill(rna_has_int.begin(), rna_has_int.end(), false);
  for (uint i=0; i!=predicted_int.size(); ++i)
    FOREACH (VU::const_iterator, j, predicted_int[i])
      rna_has_int[*j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      extract_rna_feature(rna, j, g);
}

struct FeatureCounter
{
  FeatureCounter(std::vector<FC >& fc, VU& tc)
    : fc_(fc), tc_(tc)
  { }

  inline void operator()(uint fgroup, const char* fname, uint i=-1u, uint j=-1u)
  {
    fc_[fgroup].insert(std::make_pair(std::string(fname),0u)).first->second++;
    tc_[fgroup]++;
  }

  std::vector<FC >& fc_;
  VU& tc_;
};

void
PRactIP::
count_feature(const AA& aa, const RNA& rna, const VVU& predicted_int, std::vector<FC>& fc, VU& tc) const
{
  FeatureCounter c(fc, tc);
  for (uint i=0; i!=predicted_int.size(); ++i)
    FOREACH (VU::const_iterator, j, predicted_int[i])
      extract_int_feature(aa, rna, i, *j, c);
  for (uint i=0; i!=predicted_int.size(); ++i)
    if (!predicted_int[i].empty())
      extract_aa_feature(aa, i, c);
  std::vector<bool> rna_has_int(predicted_int.size(), false);
  for (uint i=0; i!=predicted_int.size(); ++i)
    FOREACH (VU::const_iterator, j, predicted_int[i])
      rna_has_int[*j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      extract_rna_feature(rna, j, c);
}

static inline
float
clip(float w, float c)
{
  if(w >= 0)
    return w>c ? w-c : 0.0;
  else
    return -clip(-w, c);
}

float
PRactIP::
regularization_fobos(float eta)
{
  float sum1=0.0;
  for (uint i=0; i!=feature_weight_.size(); ++i)
  {
    FM::iterator j=feature_weight_[i].begin();
    while (j!=feature_weight_[i].end()) {
      //std::cout << i << ":" << j->first << "=" << j->second;
      j->second = clip(j->second, eta*lambda_);
      //std::cout << " => " << j->second << std::endl;
      sum1 += std::abs(j->second);
      if (j->second==0.0) {
        feature_weight_[i].erase(j++);
      } else
        ++j;
    }
  }

  float sum2=0.0;
  FOREACH (VF::iterator, fg, feature_group_weight_) {
    *fg /= (1.0+2.0*eta*mu_);
    sum2 += *fg * *fg;
  }
  
  return lambda_ * sum1 + mu_ * std::sqrt(sum2);
}

float
PRactIP::
predict_interaction(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p) const
{
  const uint aa_len = int_weight.size();
  const uint rna_len = int_weight[0].size();
  
  IP ip(IP::MAX, n_th_);

  VI x(aa_len, -1);             // binding site in AA
  VI y(rna_len, -1);            // binding site in RNA
  VVI z(aa_len, VI(rna_len, -1)); // interactions
  for (uint i=0; i!=aa_len; ++i)
    if (aa_weight[i]>0.0)
      x[i] = ip.make_variable(aa_weight[i]);
  for (uint j=0; j!=rna_len; ++j)
    if (rna_weight[j]>0.0)
      y[j] = ip.make_variable(rna_weight[j]);
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (int_weight[i][j]>0.0)
        z[i][j] = ip.make_variable(int_weight[i][j]);
  
  ip.update();

  for (uint i=0; i!=aa_len; ++i) {
    if (x[i]>=0) {
      //x[i] >= z[i][j]
      for (uint j=0; j!=rna_len; ++j) {
        if (z[i][j]>=0) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, x[i], 1);
          ip.add_constraint(row, z[i][j], -1);
        }
      }

      //sum_j z[i][j] >= x[i]
      int row = ip.make_constraint(IP::LO, 0, 0);
      ip.add_constraint(row, x[i], -1);
      for (uint j=0; j!=rna_len; ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }
  }

  for (uint j=0; j!=rna_len; ++j) {
    if (y[j]>=0) {
      //y[j] >= z[i][j]
      for (uint i=0; i!=aa_len; ++i) {
        if (z[i][j]>=0) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, y[j], 1);
          ip.add_constraint(row, z[i][j], -1);
        }
      }

      //sum_i z[i][j] >= y[j]
      int row = ip.make_constraint(IP::LO, 0, 0);
      ip.add_constraint(row, y[j], -1);
      for (uint i=0; i!=aa_len; ++i)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }
  }

  for (uint i=0; i!=aa_len; ++i) {
    //sum_j z[i][j] <= 1
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j=0; j!=rna_len; ++j) 
      if (z[i][j]>=0)
        ip.add_constraint(row, z[i][j], 1);
  }

  for (uint j=0; j!=rna_len; ++j) {
    //sum_i z[i][j] <= 1
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint i=0; i!=aa_len; ++i)
      if (z[i][j]>=0)
        ip.add_constraint(row, z[i][j], 1);
  }

  float s = ip.solve();
  
  p.resize(aa_len);
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (z[i][j]>=0 && ip.get_value(z[i][j])>0.5)
        p[i].push_back(j);

  return s;
}

int
PRactIP::AA::
read(const std::string& filename)
{
  char line[BUFSIZ];
  FILE* fp = fopen((filename+".seq").c_str(), "r");
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]!='>') {
      char* nr = strchr(line, '\n');
      if (nr) *nr='\0';
      this->seq+=line;
    }
  }
  fclose(fp);

#if 0
  fp = fopen((filename+".2nd").c_str(), "r");
  while (fgets(line, sizeof(line), fp)) {
    for (uint i=0; line[i]!='\0'; ++i) {
      switch (line[i]) {
        case 'C': line[i]='.'; break;
        case 'E': line[i]='>'; break;
        case 'H': line[i]='='; break;
        case '\n': line[i]='\0'; break;
      }
    }
    this->ss += line;
  }
#else
  fp = fopen((filename+".ss2").c_str(), "r");
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]=='\n' || line[0]=='#') continue;
    switch (line[7]) {
      case 'C': case 'E': case 'H':
        ss.push_back(line[7]);
        break;
      default: assert(!"unreachable"); break;
    }
  }
#endif
  fclose(fp);

  g10.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g10.begin(), group10);
  g8.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g8.begin(), group8);
  g4.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g4.begin(), group4);
  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(this->seq.size()==this->ss.size());
  return this->seq.size();
}

int
PRactIP::RNA::
read(const std::string& filename) 
{
  char line[BUFSIZ];
  FILE* fp = fopen(filename.c_str(), "r");
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]!='>') {
      char* nr = strchr(line, '\n');
      if (nr) *nr='\0';
      this->seq+=line;
    }
  }
  fclose(fp);

  std::string ss_file = filename.substr(0, filename.find_last_of('.'))+".ss";
  fp = fopen(ss_file.c_str(), "r");
  if (fp==NULL) {
    const char* prog=getenv("CENTROID_FOLD");
    if (!prog) prog="centroid_fold";
    char cmd[1000];
    sprintf(cmd, "%s %s", prog, filename.c_str());
    fp = popen(cmd, "r");
  }
  
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]=='>' || line[0]=='f') continue;
    if (strchr(".()", line[0])) {
      for (uint i=0; line[i]!='\0'; ++i) {
        switch (line[i]) {
          case '.':
          case '(':
          case ')':
            break;
          case '\n':
          case ' ':
            line[i]='\0';
            break;
        }
      }
      this->ss+=line;
    }
  }
  fclose(fp);

  structural_profile(ss, ss);

  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(this->seq.size()==this->ss.size());
  return this->seq.size();
}

// groups defined by Murphy et al., Protein Eng., 2000, 13(3), pp.149-152
// static
char
PRactIP::AA::
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
PRactIP::AA::
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
PRactIP::AA::
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
PRactIP::AA::
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

//static
char
PRactIP::RNA::
group2(char b)
{
  switch (b) {
    case 'A': case 'G':
      return 'R'; break;
    case 'C': case 'T': case 'U':
      return 'Y'; break;
  }
  return b;
}

// static
void
PRactIP::RNA::
structural_profile(const std::string& ss, std::string& profile)
{
  VU p(ss.size(), -1u);
  std::string ss2(ss);
  std::stack<uint> st;
  for (uint i=0; i!=ss2.size(); ++i) {
    switch (ss2[i]) {
      case '(':
        st.push(i);
        break;
      case ')':
        p[st.top()]=i;
        p[i]=st.top();
        st.pop();
        if (ss2[i+1]!=')' || st.top()+1!=p[i]) {
          ss2[p[i]]='[';
          ss2[i]=']';
        }
        break;
      default:
        break;
    }
  }

  profile.resize(ss.size());
  std::fill(profile.begin(), profile.end(), 'E');
  std::stack<uint> loop_degree;
  std::stack<bool> bulge;
  for (uint i=0; i!=ss2.size(); ++i) {
    switch (ss2[i]) {
      case '(':
        profile[i]='S';
        break;
      case ')':
        profile[i]='S';
        if (ss2[i-1]==']') bulge.top()=true;
        break;
      case '[':
        profile[i]='S';
        if (i>0 && (ss2[i-1]=='(' || ss2[i-1]=='[')) bulge.top()=true;
        loop_degree.push(1);
        bulge.push(false);
        break;
      case ']':
        char ps;
        profile[i]='S';
        if (ss2[i-1]==']') bulge.top()=true;
        switch (loop_degree.top()) {
          case 1: ps='H'; break;
          case 2: ps=bulge.top() ? 'B' : 'I'; break;
          default: ps='M'; break;
        }
        loop_degree.pop();
        loop_degree.top();
        if (!loop_degree.empty()) loop_degree.top()++;
        bulge.pop();
        for (uint j=p[i]+1; j!=i; ++j)
          if (profile[j]=='E') profile[j]=ps;
        break;
    }
  }
}

PRactIP&
PRactIP::
parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);

  pos_w_ = args_info.pos_w_arg;
  neg_w_ = args_info.neg_w_arg;
  lambda_ = args_info.discriminative_arg;
  mu_ = args_info.generative_arg;
  eta0_ = args_info.eta_arg;
  d_max_ = args_info.d_max_arg;
  g_max_ = args_info.g_max_arg;
  cv_fold_ = args_info.cross_validation_arg;
  
  if (args_info.inputs_num==0)
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }

  if (args_info.inputs_num>0)
    load_labeled_data(args_info.inputs[0]);
  if (args_info.inputs_num>1)
    load_unlabeled_data(args_info.inputs[1]);

  cmdline_parser_free(&args_info);
  return *this;
}

int
PRactIP::
run()
{
  if (cv_fold_>0) {
    if (unlabeled_aa_.size()>0)
      semisupervised_cross_validation(cv_fold_);
    else
      supervised_cross_validation(cv_fold_);
  } else {
    assert(!"not implemented yet");
  }
  return 0;
}

int
main(int argc, char* argv[])
{
  try {
    PRactIP practip;
    return practip.parse_options(argc, argv).run();
  } catch (const char* str) {
    std::cerr << str << std::endl;
  }
}
