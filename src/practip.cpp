#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <stack>
#include <functional>
#include <numeric>
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

//static
uint PRactIP::epoch = 0;

static
const char *groupname[] = {
  "Pro_3",                      // FG_P_3,
  "Pro_5",                      // FG_P_5,
  "RNA_3",                      // FG_R_3,
  "RNA_5",                      // FG_R_5,
  "Pro_3-RNA_3",                // FG_P_3_R_3,
  "Pro_5-RNA_5",                // FG_P_5_R_5,
  "ProSS_3",                    // FG_Pss_3,
  "ProSS_5",                    // FG_Pss_5,
  "RNASS_3",                    // FG_Rss_3,
  "RNASS_5",                    // FG_Rss_5,
  "ProSS_3-RNASS_3",            // FG_Pss_3_Rss_3,
  "ProSS_5-RNASS_5",            // FG_Pss_5_Rss_5,
  "ProG10_5",                   // FG_Pg10_5,
  "ProG10_7",                   // FG_Pg10_7,
  "ProG10_3-RNA_3",             // FG_Pg10_3_R_3,
  "ProG10_3-RNASS_3",           // FG_Pg10_3_Rss_3,
  "ProG10_5-RNA_5",             // FG_Pg10_5_R_5,
  "ProG10_5-RNASS_5",           // FG_Pg10_5_Rss_5,
  "ProG4_5",                    // FG_Pg4_5,
  "ProG4_7",                    // FG_Pg4_7,
  "ProG4_3-RNA_3",              // FG_Pg4_3_R_3,
  "ProG4_3-RNASS_3",            // FG_Pg4_3_Rss_3,
  "ProG4_5-RNA_5",              // FG_Pg4_5_R_5,
  "ProG4_5-RNASS_5",            // FG_Pg4_5_Rss_5,
};

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
    fval += acc.fval;
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

  void print(std::ostream& os) const
  {
    os << "PPV=" << ppv_avg << "(" << ppv_sd << ") "
       << "SEN=" << sen_avg << "(" << sen_sd << ") "
       << "F=" << fval_avg << "(" << fval_sd << ") ["
       << tp << "," << tn << "," << fp << "," << fn << "]"
       << std::endl;
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
    for (auto j : correct_int[i])
      c_int[rna.seq.size()*i+j]=c_aa[i]=c_rna[j]=true;
  for (uint i=0; i!=predicted_int.size(); ++i)
    for (auto j : predicted_int[i])
      p_int[rna.seq.size()*i+j]=p_aa[i]=p_rna[j]=true;
  int_acc.calculate(c_int, p_int);
  aa_acc.calculate(c_aa, p_aa);
  rna_acc.calculate(c_rna, p_rna);
}

uint
PRactIP::
load_labeled_data(const std::string& filename)
{
  std::string aa_seq, aa_ss, rna_seq, rna_ss, matching;
  std::ifstream is(filename.c_str());
  std::cout << "loading labeled data" << std::endl;
  while (is >> aa_seq >> aa_ss >> rna_seq >> rna_ss >> matching) {
    std::cout << aa_seq << " " << aa_ss << " " << rna_seq << " " << rna_ss << " " << matching << std::endl;
    labeled_aa_.push_back(AA());
    uint aa_len=labeled_aa_.back().read(aa_seq, aa_ss);
    labeled_rna_.push_back(RNA());
    labeled_rna_.back().read(rna_seq, rna_ss);
    labeled_int_.push_back(VVU(aa_len));
    read_correct_interaction(matching, labeled_int_.back());
  }
  return labeled_aa_.size();
}

uint
PRactIP::
load_unlabeled_data(const std::string& filename)
{
  std::string aa_seq1, aa_ss1, rna_seq1, rna_ss1;
  std::string aa_seq2, aa_ss2, rna_seq2, rna_ss2;
  std::string aa_al, aa_sc;
  std::string rna_al, rna_sc;
  std::ifstream is(filename.c_str());
  std::cout << "loading unlabeled data" << std::endl;
  while (is >> aa_seq1 >> aa_ss1 >> rna_seq1 >> rna_ss1 >>
         aa_seq2 >> aa_ss2 >> rna_seq2 >> rna_ss2 >> 
         aa_al >> aa_sc >> rna_al >> rna_sc) {
    std::cout << aa_seq1 << " " << aa_ss1 << " " << rna_seq1 << " " << rna_ss2 << " "
              << aa_seq2 << " " << aa_ss2 << " " << rna_seq2 << " " << rna_ss2 << " "
              << aa_al << " " << aa_sc << " " 
              << rna_al << " " << rna_sc << std::endl;

    unlabeled_aa_.push_back(Alignment<AA>(aa_al.c_str(), aa_sc.c_str()));
    unlabeled_aa_.back().add_seq(AA(aa_seq1, aa_ss1));
    unlabeled_aa_.back().add_seq(AA(aa_seq2, aa_ss2));
    unlabeled_rna_.push_back(Alignment<RNA>(rna_al.c_str(), rna_sc.c_str()));
    unlabeled_rna_.back().add_seq(RNA(rna_seq1, rna_ss1));
    unlabeled_rna_.back().add_seq(RNA(rna_seq2, rna_ss2));
  }
  return unlabeled_aa_.size();
}

float
PRactIP::
calculate_score(const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, const VVU& interactions) const
{
  float score = 0.0;
  for (uint i=0; i!=interactions.size(); ++i)
    for (auto j : interactions[i])
      score += int_weight[i][j];
  for (uint i=0; i!=interactions.size(); ++i)
    score += aa_weight[i];
  std::vector<bool> rna_has_int(rna_weight.size(), false);
  for (uint i=0; i!=interactions.size(); ++i)
    for (auto j : interactions[i])
      rna_has_int[j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      score += rna_weight[j];
  return score;
}


float
PRactIP::
supervised_training(const AA& aa, const RNA& rna, const VVU& correct_int, bool max_margin /*=true*/, float w /*=1.0*/)
{
  float loss=0.0;
  VVF int_weight;
  VF aa_weight, rna_weight;
  calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  loss -= calculate_score(int_weight, aa_weight, rna_weight, correct_int);
  if (max_margin)
    penalize_correct_interaction(int_weight, aa_weight, rna_weight, correct_int);
  VVU predicted_int;
  predict_interaction(aa, rna, int_weight, aa_weight, rna_weight, predicted_int);
  loss += calculate_score(int_weight, aa_weight, rna_weight, predicted_int);
#if 0
  for (uint i=0; i!=correct_int.size(); ++i) {
    if (!correct_int[i].empty() || !predicted_int[i].empty())
    {
      std::cout << i << ": [ ";
      for (auto j : correct_int[i])
        std::cout << j << "(" << int_weight[i][j] << ") ";
      std::cout << "], [ ";
      for (auto j : predicted_int[i])
        std::cout << j << "(" << int_weight[i][j] << ") ";
      std::cout << "]" << std::endl;
    }
  }
  std::cout << std::endl;
#endif  
  update_feature_weight(aa, rna, predicted_int, correct_int, w);

  return loss;
}

void
PRactIP::
cross_validation(uint n)
{
  AccuracySummary int_summary, aa_summary, rna_summary;
  AccuracySummary int_summary_tr, aa_summary_tr, rna_summary_tr;
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
    semisupervised_training(train);

    for (auto j : test) {
      VVU predicted_int;
      predict_interaction(labeled_aa_[j], labeled_rna_[j], predicted_int);
      Accuracy int_acc, aa_acc, rna_acc;
      calculate_accuracy(labeled_aa_[j], labeled_rna_[j],
                         predicted_int, labeled_int_[j],
                         int_acc, aa_acc, rna_acc);
      int_summary.add(int_acc);
      aa_summary.add(aa_acc);
      rna_summary.add(rna_acc);
    }
    for (auto j : train) {
      VVU predicted_int;
      predict_interaction(labeled_aa_[j], labeled_rna_[j], predicted_int);
      Accuracy int_acc, aa_acc, rna_acc;
      calculate_accuracy(labeled_aa_[j], labeled_rna_[j],
                         predicted_int, labeled_int_[j],
                         int_acc, aa_acc, rna_acc);
      int_summary_tr.add(int_acc);
      aa_summary_tr.add(aa_acc);
      rna_summary_tr.add(rna_acc);
    }
  }

  std::cout << "[training data]" << std::endl;
  std::cout << "Interaction: ";
  int_summary_tr.summary(); int_summary_tr.print(std::cout);
  std::cout << "AA: ";
  aa_summary_tr.summary(); aa_summary_tr.print(std::cout);
  std::cout << "RNA: ";
  rna_summary_tr.summary(); rna_summary_tr.print(std::cout);

  std::cout << "[test data]" << std::endl;
  std::cout << "Interaction: ";
  int_summary.summary(); int_summary.print(std::cout);
  std::cout << "AA: ";
  aa_summary.summary(); aa_summary.print(std::cout);
  std::cout << "RNA: ";
  rna_summary.summary(); rna_summary.print(std::cout);
}

void
PRactIP::
semisupervised_training()
{
  VU idx(labeled_aa_.size());
  std::iota(std::begin(idx), std::end(idx), 0);
  semisupervised_training(idx);
}

void
PRactIP::
semisupervised_training(const VU& use_idx)
{
  epoch=0;
  // initial supervised learning
  VU idx(use_idx);
  for (uint t=0; t!=d_max_; ++t) {
    float total_loss=0.0;
    std::random_shuffle(std::begin(idx), std::end(idx)); // shuffle the order of training data
    for (auto i : idx) {
      total_loss += supervised_training(labeled_aa_[i], labeled_rna_[i], labeled_int_[i]);
      epoch++;
    }
  }

  // semisupervised learning
  if (unlabeled_aa_.size()>0)
  {
    idx.resize(use_idx.size()+unlabeled_aa_.size());
    for (uint i=0; i!=idx.size(); ++i)
      idx[i] = i<use_idx.size() ? use_idx[i] : i-use_idx.size()+labeled_aa_.size();

    for (uint u=0; u!=g_max_; ++u)
    {
      std::random_shuffle(idx.begin(), idx.end());

      float unlabeled_loss=0.0;
      float total_loss=0.0;
      for (auto i : idx)
      {
        if (i<labeled_aa_.size()) // labeled data
        {
          //std::cout << ">> supervised" << std::endl;
          total_loss += supervised_training(labeled_aa_[i], labeled_rna_[i], labeled_int_[i]);
          epoch++;
        }
        else                    // unlabeled data
        {
          const auto iu=i-labeled_aa_.size();
          // common structure prediction
          VVVU common_int;
          predict_common_interaction(unlabeled_aa_[iu], unlabeled_rna_[iu], common_int);
          assert(common_int.size()==unlabeled_aa_[iu].num_sequences());
          assert(common_int.size()==unlabeled_rna_[iu].num_sequences());
          // train from predicted common structures
          for (uint k=0; k!=common_int.size(); ++k)
          {
            //std::cout << ">> semi-supervised" << std::endl;
            unlabeled_loss += supervised_training(unlabeled_aa_[iu].seq(k), unlabeled_rna_[iu].seq(k), common_int[k], false, mu_);
            epoch++;
          }
        }
      }
    }
  }

  regularization_fobos();
}

float
PRactIP::
predict_interaction(const AA& aa, const RNA& rna, VVU& predicted_int, float w /*=1.0*/)
{
  VVF int_weight;
  VF aa_weight, rna_weight;
  calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  return predict_interaction(aa, rna, int_weight, aa_weight, rna_weight, predicted_int, w);
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
extract_int_feature(const AA& aa, const RNA& rna, uint i, uint j, Func func) const
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
extract_aa_feature(const AA& aa, uint i, Func func) const
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
extract_rna_feature(const RNA& rna, uint j, Func func) const
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

void
PRactIP::
store_parameters(const char* filename) const
{
  std::ofstream os(filename);
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  for (uint k=0; k!=feature_weight_.size(); ++k)
  {
    os << "[ " << groupname[k] << " weight ]" << std::endl;
    for (const auto& e : feature_weight_[k])
    {
      if (e.second.weight!=0.0)
        os << e.first << " " << e.second.weight << std::endl;
    }
    os << std::endl;
  }
}

void
PRactIP::
restore_parameters(const char* filename)
{
  std::ifstream is(filename);
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  for (auto f : feature_weight_)
    f.clear();

  uint k=-1u;
  std::string first, second;
  while (is >> first >> second)
  {
    if (first=="[")             // group name
    {
      //std::cout << second << std::endl;
      std::string temp1, temp2;
      is >> temp1 >> temp1;
      auto p = std::find(std::begin(groupname), std::end(groupname), second);
      if (p!=std::end(groupname))
        k = p - std::begin(groupname);
      else
        throw std::runtime_error(std::string("unknown groupname: ") + second);
    }
    else if (k!=-1u)
    {
      feature_weight_[k][first].weight = std::atof(second.c_str());
    }
  }
}

void
PRactIP::
default_parameters()
{
  struct { const char* name; float value; } params[] = {
#include "defparam.dat"
  };

  for (auto f : feature_weight_)
    f.clear();

  uint k=-1u;
  for (const auto& x : params) {
    if (x.value==0.0)
    {
      auto p = std::find(std::begin(groupname), std::end(groupname), x.name);
      if (p!=std::end(groupname))
        k = p - std::begin(groupname);
      else
        throw std::runtime_error(std::string("unknown groupname: ") + x.name);
    }
    else if (k!=-1u)
    {
      feature_weight_[k][x.name].weight = x.value;
    }
  }
}

void
PRactIP::
calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight)
{
  int_weight.resize(aa.seq.size());
  for (uint i=0; i!=int_weight.size(); ++i) 
  {
    int_weight[i].resize(rna.seq.size());
    for (uint j=0; j!=int_weight[i].size(); ++j) 
    {
      int_weight[i][j] = 0.0;
      extract_int_feature(aa, rna, i, j, 
                          [&](uint fgroup, const char* fname, uint i, uint j) 
                          {
                            int_weight[i][j] += update_fobos(fgroup, fname);
                          }
        );
    }
  }

  aa_weight.resize(aa.seq.size());
  for (uint i=0; i!=aa_weight.size(); ++i) 
  {
    aa_weight[i] = 0.0;
    extract_aa_feature(aa, i,
                       [&](uint fgroup, const char* fname, uint i) 
                       {
                         aa_weight[i] += update_fobos(fgroup, fname);
                       }
      );
  }

  rna_weight.resize(rna.seq.size());
  for (uint j=0; j!=rna_weight.size(); ++j) 
  {
    rna_weight[j] = 0.0;
    extract_rna_feature(rna, j, 
                       [&](uint fgroup, const char* fname, uint j) 
                        {
                          rna_weight[j] += update_fobos(fgroup, fname);
                        }
      );
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
    for (auto j : correct_int[i])
      int_weight[i][j] -= pos_w_+neg_w_;

  for (uint i=0; i!=aa_weight.size(); ++i) 
    aa_weight[i] += neg_w_;

  for (uint i=0; i!=correct_int.size(); ++i)
    if (!correct_int[i].empty())
      aa_weight[i] -= pos_w_+neg_w_;

  for (uint j=0; j!=rna_weight.size(); ++j)
    rna_weight[j] += neg_w_;

  std::vector<bool> rna_has_int(rna_weight.size(), false);
  for (uint i=0; i!=correct_int.size(); ++i)
    for (auto j : correct_int[i])
      rna_has_int[j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
      rna_weight[j] -= pos_w_+neg_w_;
}

void
PRactIP::
update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_int, const VVU& correct_int, float w /*=1.0*/)
{
  assert(predicted_int.size()==correct_int.size());

  // calculate gradients in correct interactions
  //   interactions
  typedef std::unordered_map<std::string,int> GM;
  std::vector<GM> gr(FG_NUM);
  for (uint i=0; i!=correct_int.size(); ++i)
    for (auto j : correct_int[i])
    {
      extract_int_feature(aa, rna, i, j,
                          [&] (uint fgroup, const char* fname, uint i, uint j)
                          {
                            gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += -1;
                          }
        );
    }
  //   amino acids
  for (uint i=0; i!=correct_int.size(); ++i)
    if (!correct_int[i].empty())
    {
      extract_aa_feature(aa, i, 
                         [&] (uint fgroup, const char* fname, uint i)
                         {
                           gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += -1;
                         }
        );
    }
  //   RNAs
  std::vector<bool> rna_has_int(rna.seq.size(), false);
  for (uint i=0; i!=correct_int.size(); ++i)
    for (auto j : correct_int[i])
      rna_has_int[j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
    {
      extract_rna_feature(rna, j, 
                          [&] (uint fgroup, const char* fname, uint j)
                          {
                            gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += -1;
                          }
        );
    }

  // calculate gradients in predicted interactions
  //   interactions
  for (uint i=0; i!=predicted_int.size(); ++i) 
    for (auto j : predicted_int[i])
    {
      extract_int_feature(aa, rna, i, j, 
                          [&] (uint fgroup, const char* fname, uint i, uint j)
                          {
                            gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += +1;
                          }
        );
    }
  //   amino acids
  for (uint i=0; i!=predicted_int.size(); ++i)
    if (!predicted_int[i].empty())
    {
      extract_aa_feature(aa, i, 
                         [&] (uint fgroup, const char* fname, uint i)
                         {
                           gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += +1;
                         }
        );
    }
  //   RNAs
  std::fill(rna_has_int.begin(), rna_has_int.end(), false);
  for (uint i=0; i!=predicted_int.size(); ++i)
    for (auto j : predicted_int[i])
      rna_has_int[j] = true;
  for (uint j=0; j!=rna_has_int.size(); ++j)
    if (rna_has_int[j])
    {
      extract_rna_feature(rna, j, 
                          [&] (uint fgroup, const char* fname, uint j)
                          {
                            gr[fgroup].insert(std::make_pair(std::string(fname), 0)).first->second += +1;
                          }
        );
    }

  // update feature weights by AdaGrad
  for (uint k=0; k!=gr.size(); ++k)
  {
    for (const auto& e : gr[k])
    {
      if (e.second!=0)
      {
        auto fe = feature_weight_[k].insert(std::make_pair(e.first, PRactIP::FeatureWeight())).first;
        fe->second.weight -= w*e.second * eta0_/std::sqrt(1.0+fe->second.sum_of_grad2);
        fe->second.sum_of_grad2 += w*e.second * w*e.second;
#if 0
        std::cerr << epoch << ", "
                  << k << ", "
                  << it->first << ", "
                  << it->second << ", "
                  << fe->second.weight << ", "
                  << eta0_/std::sqrt(1.0+fe->second.sum_of_grad2) << std::endl;
#endif
        //assert(fe->second.last_updated == epoch);
      }
    }
  }
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
update_fobos(uint fgroup, const char* fname)
{
  auto m=feature_weight_[fgroup].find(fname);
  if (m!=feature_weight_[fgroup].end())
  {
    // lazy update for FOBOS
    if (m->second.last_updated<epoch)
    {
      const float eta = eta0_/std::sqrt(1.0+m->second.sum_of_grad2);
      const uint t = epoch - m->second.last_updated;
      m->second.weight = clip(m->second.weight, lambda_*eta*t);
      m->second.last_updated = epoch;
    }
    return m->second.weight;
  }
  return 0.0;
}

float
PRactIP::
regularization_fobos()
{
  // L1-norm
  float sum1=0.0;
  for (uint k=0; k!=feature_weight_.size(); ++k)
  {
    for (auto& e : feature_weight_[k])
    {
      if (e.second.last_updated<epoch)
      {
        const float eta = eta0_/std::sqrt(1.0+e.second.sum_of_grad2);
        const uint t = epoch - e.second.last_updated;
        e.second.weight = clip(e.second.weight, lambda_*eta*t);
        e.second.last_updated = epoch;
      }
      sum1 += std::abs(e.second.weight);
    }
  }
  return lambda_ * sum1;
}

void
PRactIP::
predict_interaction_object(const AA& aa, const RNA& rna,
                           const VVF& int_weight, const VF& aa_weight, const VF& rna_weight,
                           VI& x, VI& y, VVI& z, VI& sl_x, VI& sl_y, IP& ip, float w /*=1.0*/) const
{
  const uint aa_len = x.size();
  const uint rna_len = y.size();
  const int MAX_INTERACTION = 20;

  for (uint i=0; i!=aa_len; ++i)
    if (aa_weight[i]>0.0) {
      x[i] = ip.make_variable(w*aa_weight[i]);
      sl_x[i] = ip.make_variable(-exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint j=0; j!=rna_len; ++j)
    if (rna_weight[j]>0.0) {
      y[j] = ip.make_variable(w*rna_weight[j]);
      sl_y[j] = ip.make_variable(-exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (int_weight[i][j]>0.0)
        z[i][j] = ip.make_variable(w*int_weight[i][j]);
}

void
PRactIP::
predict_interaction_constraints(const AA& aa, const RNA& rna,
                                const VI& x, const VI& y, const VVI& z, const VI& sl_x, const VI& sl_y, 
                                IP& ip) const
{
  const uint aa_len = x.size();
  const uint rna_len = y.size();

  for (uint i=0; i!=aa_len; ++i) {
    if (x[i]>=0) {
      //x_i >= z_ij
      for (uint j=0; j!=rna_len; ++j) {
        if (z[i][j]>=0) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, x[i], 1);
          ip.add_constraint(row, z[i][j], -1);
        }
      }

      //sum_j z_ij >= x_i
      int row = ip.make_constraint(IP::LO, 0, 0);
      ip.add_constraint(row, x[i], -1);
      for (uint j=0; j!=rna_len; ++j)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }
  }

  for (uint j=0; j!=rna_len; ++j) {
    if (y[j]>=0) {
      //y_j >= z_ij
      for (uint i=0; i!=aa_len; ++i) {
        if (z[i][j]>=0) {
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, y[j], 1);
          ip.add_constraint(row, z[i][j], -1);
        }
      }

      //sum_i z_ij >= y_j
      int row = ip.make_constraint(IP::LO, 0, 0);
      ip.add_constraint(row, y[j], -1);
      for (uint i=0; i!=aa_len; ++i)
        if (z[i][j]>=0)
          ip.add_constraint(row, z[i][j], 1);
    }
  }

  for (uint i=0; i!=aa_len; ++i) {
    //sum_j z_ij <= max_int(ss_i)
    int row = ip.make_constraint(IP::UP, 0, aa_int_max_==-1u ? AA::max_intraction(aa.ss[i]) : aa_int_max_);
    if (sl_x[i]>=0) ip.add_constraint(row, sl_x[i], -1);
    for (uint j=0; j!=rna_len; ++j) 
      if (z[i][j]>=0)
        ip.add_constraint(row, z[i][j], 1);
  }

  for (uint j=0; j!=rna_len; ++j) {
    //sum_i z_ij <= max_int(ss_j)
    int row = ip.make_constraint(IP::UP, 0, rna_int_max_==-1u ? RNA::max_intraction(rna.ss[j]) : rna_int_max_);
    if (sl_y[j]>=0) ip.add_constraint(row, sl_y[j], -1);
    for (uint i=0; i!=aa_len; ++i)
      if (z[i][j]>=0)
        ip.add_constraint(row, z[i][j], 1);
  }

  for (uint j=0; j!=rna_len; ++j) {
    if (y[j]>=0) {
      // y_{j-1}+(1-y_j)+y_{j+1}>=1
      int row = ip.make_constraint(IP::LO, 0, 0);
      if (j>=1 && y[j-1]>=0) ip.add_constraint(row, y[j-1], 1);
      ip.add_constraint(row, y[j], -1);
      if (j+1<rna_len && y[j+1]>=0) ip.add_constraint(row, y[j+1], 1);
    }
  }
}

float
PRactIP::
predict_interaction(const AA& aa, const RNA& rna,
                    const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p, float w /*=1.0*/) const
{
  const uint aa_len = int_weight.size();
  const uint rna_len = int_weight[0].size();
  
  IP ip(IP::MAX, n_th_);

  VI x(aa_len, -1);             // binding site in AA
  VI y(rna_len, -1);            // binding site in RNA
  VVI z(aa_len, VI(rna_len, -1)); // interactions
  VI sl_x(aa_len, -1);            // slack variables for AA to relax some constraints
  VI sl_y(rna_len, -1);           // slack variables for RNA to relax some constraints

  predict_interaction_object(aa, rna, int_weight, aa_weight, rna_weight, x, y, z, sl_x, sl_y, ip, w);
  ip.update();
  predict_interaction_constraints(aa, rna, x, y, z, sl_x, sl_y, ip);
  float s = ip.solve();

  p.resize(aa_len);
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (z[i][j]>=0 && ip.get_value(z[i][j])>0.5)
        p[i].push_back(j);

  return s;
}

float
PRactIP::
predict_common_interaction(const Alignment<AA>& aa, const Alignment<RNA>& rna, VVVU& predicted_int)
{
  const uint n_seq=aa.num_sequences();
  assert(n_seq==rna.num_sequences());

  VVVF int_weight(n_seq);
  VVF aa_weight(n_seq), rna_weight(n_seq);
  for (uint i=0; i!=n_seq; ++i)
    calculate_feature_weight(aa.seq(i), rna.seq(i), int_weight[i], aa_weight[i], rna_weight[i]);
  return predict_common_interaction(aa, rna, int_weight, aa_weight, rna_weight, predicted_int);
}

float
PRactIP::
predict_common_interaction(const Alignment<AA>& aa, const Alignment<RNA>& rna, 
                           const VVVF& int_weight, const VVF& aa_weight, const VVF& rna_weight, 
                           VVVU& predicted_int) const
{
  const uint n_seq=aa.num_sequences();

  IP ip(IP::MAX, n_th_);

  VVI x(n_seq);          // binding site in AA
  VVI y(n_seq);          // binding site in RNA
  VVVI z(n_seq);         // interactions
  VVI sl_x(n_seq);       // slack variables for AA to relax some constraints
  VVI sl_y(n_seq);       // slack variables for RNA to relax some constraints
  for (uint k=0; k!=n_seq; ++k)
  {
    const uint aa_len = int_weight[k].size();
    const uint rna_len = int_weight[k][0].size();
    x[k].resize(aa_len, -1);
    y[k].resize(rna_len, -1);
    z[k].resize(aa_len, VI(rna_len, -1));
    sl_x[k].resize(aa_len, -1);
    sl_y[k].resize(rna_len, -1);

    predict_interaction_object(aa.seq(k), rna.seq(k), int_weight[k], aa_weight[k], rna_weight[k], x[k], y[k], z[k], sl_x[k], sl_y[k], ip, mu_);
  }

  VVI x_al(n_seq*(n_seq-1)/2, VI(aa.num_columns(), -1));
  VVI y_al(n_seq*(n_seq-1)/2, VI(rna.num_columns(), -1));
  VVVI z_al(n_seq*(n_seq-1)/2, VVI(aa.num_columns(), VI(rna.num_columns(), -1)));
  for (uint k=0, m=0; k!=n_seq-1; ++k) 
  {
    const auto& aa_k=aa.idx()[k];
    const auto& rna_k=rna.idx()[k];
    for (uint l=k+1; l!=n_seq; ++l, ++m) 
    {
      const auto& aa_l=aa.idx()[l];
      const auto& rna_l=rna.idx()[l];
      for (uint i=0; i!=aa.num_columns(); ++i)
        if (aa.q_score()[i]>0.0 && 
            aa_k[i]!=-1u && x[k][aa_k[i]]>=0 && 
            aa_l[i]!=-1u && x[l][aa_l[i]]>=0)
          x_al[m][i] = ip.make_variable(-nu_*aa.q_score()[i]);
      for (uint j=0; j!=rna.num_columns(); ++j)
        if (rna.q_score()[j]>0.0 &&
            rna_k[j]!=-1u && y[k][rna_k[j]]>=0 &&
            rna_l[j]!=-1u && y[l][rna_l[j]]>=0)
          y_al[m][j] = ip.make_variable(-nu_*rna.q_score()[j]);
      for (uint i=0; i!=aa.num_columns(); ++i)
        if (aa.q_score()[i]>0.0 && aa_k[i]!=-1u && aa_l[i]!=-1u)
          for (uint j=0; j!=rna.num_columns(); ++j)
            if (rna.q_score()[j]>0.0 && 
                rna_k[j]!=-1u && z[k][aa_k[i]][rna_k[j]]>=0 && 
                rna_l[j]!=-1u && z[l][aa_l[i]][rna_l[j]]>=0)
              z_al[m][i][j] = ip.make_variable(-nu_*aa.q_score()[i]*rna.q_score()[j]);
    }
  }
  ip.update();

  for (uint k=0; k!=n_seq; ++k)
    predict_interaction_constraints(aa.seq(k), rna.seq(k), x[k], y[k], z[k], sl_x[k], sl_y[k], ip);

  for (uint k=0, m=0; k!=n_seq-1; ++k) 
  {
    const auto& aa_k=aa.idx()[k];
    const auto& rna_k=rna.idx()[k];
    for (uint l=k+1; l!=n_seq; ++l, ++m) 
    {
      const auto& aa_l=aa.idx()[l];
      const auto& rna_l=rna.idx()[l];
      for (uint i=0; i!=aa.num_columns(); ++i)
      {
        if (x_al[m][i]>=0)
        {
          // x[k][i] - x[l][i] >= -x_al[m][i]
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, x_al[m][i], 1);
          ip.add_constraint(row, x[k][aa_k[i]], 1);
          ip.add_constraint(row, x[l][aa_l[i]], -1);
          // x[k][i] - x[l][i] <= x_al[m][i]
          row = ip.make_constraint(IP::UP, 0, 0);
          ip.add_constraint(row, x_al[m][i], -1);
          ip.add_constraint(row, x[k][aa_k[i]], 1);
          ip.add_constraint(row, x[l][aa_l[i]], -1);
        }
      }
      for (uint j=0; j!=rna.num_columns(); ++j)
      {
        if (y_al[m][j]>=0)
        {
          // y[k][j] - y[l][j] >= -y_al[m][j]
          int row = ip.make_constraint(IP::LO, 0, 0);
          ip.add_constraint(row, y_al[m][j], 1);
          ip.add_constraint(row, y[k][rna_k[j]], 1);
          ip.add_constraint(row, y[l][rna_l[j]], -1);
          // y[k][j] - y[l][j] <= y_al[m][j]
          row = ip.make_constraint(IP::UP, 0, 0);
          ip.add_constraint(row, y_al[m][j], -1);
          ip.add_constraint(row, y[k][rna_k[j]], 1);
          ip.add_constraint(row, y[l][rna_l[j]], -1);
        }
      }
      for (uint i=0; i!=aa.num_columns(); ++i)
      {
        for (uint j=0; j!=rna.num_columns(); ++j)
        {
          if (z_al[m][i][j]>=0)
          {
            // z[k][i][j] - z[l][i][j] >= -z_al[m][i][j]
            int row = ip.make_constraint(IP::LO, 0, 0);
            ip.add_constraint(row, z_al[m][i][j], 1);
            ip.add_constraint(row, z[k][aa_k[i]][rna_k[j]], 1);
            ip.add_constraint(row, z[l][aa_l[i]][rna_l[j]], -1);
            // z[k][i][j] - z[l][i][j] <= z_al[m][i][j]
            row = ip.make_constraint(IP::UP, 0, 0);
            ip.add_constraint(row, z_al[m][i][j], -1);
            ip.add_constraint(row, z[k][aa_k[i]][rna_k[j]], 1);
            ip.add_constraint(row, z[l][aa_l[i]][rna_l[j]], -1);
          }
        }
      }
    }
  }

  float s = ip.solve();

  predicted_int.resize(n_seq);
  for (uint k=0; k!=n_seq; ++k)
  {
    const uint aa_len = int_weight[k].size();
    const uint rna_len = int_weight[k][0].size();
    predicted_int[k].resize(aa_len);
    for (uint i=0; i!=aa_len; ++i)
      for (uint j=0; j!=rna_len; ++j)
        if (z[k][i][j]>=0 && ip.get_value(z[k][i][j])>0.5)
          predicted_int[k][i].push_back(j);
  }

  return s;
}

void
PRactIP::
show_result(const AA& aa, const RNA& rna, const VVU& predicted_int, float score) const
{
#if 0
  std::cout << ">" << aa.name << std::endl
            << aa.seq << std::endl;
  for (uint i=0; i!=aa.seq.size(); ++i)
    std::cout << i%10;
  std::cout << std::endl;
  
  std::cout << ">" << rna.name << std::endl
            << rna.seq << std::endl;
  for (uint i=0; i!=rna.seq.size(); ++i)
    std::cout << i%10;
  std::cout << std::endl;
#endif

  std::cout << ">Score=" << score << std::endl;
  for (uint i=0; i!=predicted_int.size(); ++i)
    for (auto j : predicted_int[i])
      std::cout << i << " " << j << std::endl;
}

int
PRactIP::AA::
read(const std::string& fa_name, const std::string& ss_name)
{
  read_fa(fa_name);
  read_ss(ss_name);

  g10.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g10.begin(), group10);
  g8.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g8.begin(), group8);
  g4.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g4.begin(), group4);
  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(seq.size()==ss.size());
  return seq.size();
}

int
PRactIP::AA::
read_fa(const std::string& fa_name)
{
  std::string line;
  std::ifstream ifs(fa_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + fa_name);
  while (std::getline(ifs, line)) {
    if (line[0]=='>')
      this->name+=line.substr(1);
    else
      this->seq+=line;
  }
  ifs.close();

  return this->seq.size();
}
  

int
PRactIP::AA::
read_ss(const std::string& ss_name)
{
  std::string line;
  std::ifstream ifs(ss_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + ss_name);

  std::getline(ifs, line);
  
  if (line[0]=='#') {           // suppose PSIPRED format
    while (std::getline(ifs, line)) {
      if (line[0]=='\n' || line[0]=='#') continue;
      switch (line[7]) {
        case 'E': case 'H':
          ss.push_back(line[7]);
          break;
        default:
          //assert(!"unreachable"); break;
        case 'C':         
          ss.push_back('C');
          break;
      }
    }
  } else {                      // suppose VIENNA-like format
    std::getline(ifs, line);    // this is the sequence.
    std::getline(ifs, line);    // this is the secondary structure
    for (uint i=0; i!=line.size(); ++i) {
      switch (line[i]) {
        case 'E': case 'H':
          ss.push_back(line[i]);
          break;
        default:
          //assert(!"unreachable"); break;
        case 'C':         
          ss.push_back('C');
          break;
      }
    }    
  }

  return ss.size();
}

// static
inline
uint
PRactIP::AA::
max_intraction(char x)
{
  return 3;
}

// static
inline
uint
PRactIP::RNA::
max_intraction(char x)
{
  switch (x)
  {
    case 'E': return 7; break;
    case 'H': return 5; break;
  }
  return 4;
}

int
PRactIP::RNA::
read(const std::string& fa_name, const std::string& ss_name)
{
  read_fa(fa_name);
  read_ss(ss_name);

  g2.resize(seq.size());
  std::transform(seq.begin(), seq.end(), g2.begin(), group2);

  assert(this->seq.size()==this->ss.size());
  return this->seq.size();
}

int
PRactIP::RNA::
read_fa(const std::string& fa_name) 
{
  std::string line;
  std::ifstream ifs(fa_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + fa_name);
  while (std::getline(ifs, line)) {
    if (line[0]=='>')
      this->name+=line.substr(1);
    else
      this->seq+=line;
  }
  ifs.close();

  return this->seq.size();
}

int
PRactIP::RNA::
read_ss(const std::string& ss_name) 
{
  std::string line;
  std::ifstream ifs(ss_name);
  if (!ifs)
    throw std::runtime_error(std::string(strerror(errno)) + ": " + ss_name);

#if 0
  FILE* fp = fopen((basename+".ss").c_str(), "r");
  if (fp==NULL) {
    const char* prog=getenv("CENTROID_FOLD");
    if (!prog) prog="centroid_fold";
    char cmd[1000];
    sprintf(cmd, "%s %s", prog, filename.c_str());
    fp = popen(cmd, "r");
  }

  std::string ss_str;
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
      ss_str+=line;
    }
  }
  fclose(fp);
#endif

  std::getline(ifs, line);      // this is the header;
  std::getline(ifs, line);      // this is the sequence.
  std::getline(ifs, line);      // this is the secondary structure
  structural_profile(line, ss);

  return ss.size();
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
        if (!loop_degree.empty()) loop_degree.top()++;
        bulge.pop();
        for (uint j=p[i]+1; j!=i; ++j)
          if (profile[j]=='E') profile[j]=ps;
        break;
    }
  }
}

template <class T>
PRactIP::Alignment<T>::
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
PRactIP::Alignment<T>::
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
        idx_.push_back(VU(seq.size(), -1U));
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
    idx_.push_back(VU(seq.size(), -1U));
    VU& x=idx_.back();
    for (uint i=0, p=0; i!=seq.size(); ++i)
      if (seq[i]!='-') x[i]=p++;
  }
  return true;
}

template <class T>
bool
PRactIP::Alignment<T>::
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

PRactIP&
PRactIP::
parse_options(int& argc, char**& argv)
{
  gengetopt_args_info args_info;
  if (cmdline_parser(argc, argv, &args_info)!=0) exit(1);

  if (args_info.train_given)
  {
    train_mode_ = true;
    param_file_ = args_info.train_arg;
  }
  else if (args_info.predict_given)
  {
    train_mode_ = false;
    param_file_ = args_info.predict_arg;
  }
  pos_w_ = args_info.pos_w_arg;
  neg_w_ = args_info.neg_w_arg;
  lambda_ = args_info.lambda_arg;
  mu_ = args_info.mu_arg;
  nu_ = args_info.nu_arg;
  eta0_ = args_info.eta_arg;
  d_max_ = args_info.d_max_arg;
  g_max_ = args_info.g_max_arg;
  cv_fold_ = args_info.cross_validation_arg;
  exceed_penalty_ = args_info.exceeding_penalty_arg;
  if (args_info.aa_int_max_given)
    aa_int_max_ = args_info.aa_int_max_arg;
  if (args_info.rna_int_max_given)
    rna_int_max_ = args_info.rna_int_max_arg;
  if (args_info.threads_given)
    n_th_ = args_info.threads_arg;
  
  if (args_info.inputs_num==0)
  {
    cmdline_parser_print_help();
    cmdline_parser_free(&args_info);
    exit(1);
  }

  args_.resize(args_info.inputs_num);
  for (uint i=0; i!=args_info.inputs_num; ++i)
    args_[i]=args_info.inputs[i];

  cmdline_parser_free(&args_info);

  return *this;
}

int
PRactIP::
run()
{
  if (cv_fold_>0 || train_mode_)
  {
    if (args_.size()>0)
      load_labeled_data(args_[0]);
    if (args_.size()>1)
      load_unlabeled_data(args_[1]);

    if (cv_fold_>0) {
      cross_validation(cv_fold_);
    } else if (train_mode_){
      if (unlabeled_aa_.size()>0)
        semisupervised_training();
      else
        supervised_training();
      store_parameters(param_file_.c_str());
    }
  }
  else
  {
    if (args_.size()<4)
    {
      cmdline_parser_print_help();
      return 0;
    }
    if (param_file_.empty())
      default_parameters();
    else
      restore_parameters(param_file_.c_str());

    AA aa(args_[0], args_[1]); 
    RNA rna(args_[2], args_[3]);
    VVU predicted_int;
    float s=predict_interaction(aa, rna, predicted_int);
    show_result(aa, rna, predicted_int, s);
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
