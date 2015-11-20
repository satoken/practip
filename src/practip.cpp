#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
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
calculate_accuracy(const AA& aa, const RNA& rna,
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
    labeled_aa_.emplace_back();
    uint aa_len=labeled_aa_.back().read(aa_seq, aa_ss);
    labeled_rna_.emplace_back();
    labeled_rna_.back().read(rna_seq, rna_ss);
    labeled_int_.emplace_back(aa_len);
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

    unlabeled_aa_.emplace_back(aa_al.c_str(), aa_sc.c_str());
    unlabeled_aa_.back().emplace_add_seq(aa_seq1, aa_ss1);
    unlabeled_aa_.back().emplace_add_seq(aa_seq2, aa_ss2);
    unlabeled_rna_.emplace_back(rna_al.c_str(), rna_sc.c_str());
    unlabeled_rna_.back().emplace_add_seq(rna_seq1, rna_ss1);
    unlabeled_rna_.back().emplace_add_seq(rna_seq2, rna_ss2);
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
supervised_training(FeatureManager& fm, const AA& aa, const RNA& rna, const VVU& correct_int,
                    bool max_margin /*=true*/, float w /*=1.0*/)
{
  VVF int_weight;
  VF aa_weight, rna_weight;
  fm.calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  float score_for_correct = calculate_score(int_weight, aa_weight, rna_weight, correct_int);

  if (max_margin)
    penalize_correct_interaction(int_weight, aa_weight, rna_weight, correct_int);
  VVU predicted_int;
  predict_interaction(fm, aa, rna, int_weight, aa_weight, rna_weight, predicted_int);
  float score_for_predict = calculate_score(int_weight, aa_weight, rna_weight, predicted_int);

#if 1
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

  fm.update_feature_weight(aa, rna, predicted_int, correct_int, w);

  return score_for_predict - score_for_correct;
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
    FeatureManager fm(lambda_, eta0_);
    semisupervised_training(fm, train);

    for (auto j : test) {
      VVU predicted_int;
      predict_interaction(fm, labeled_aa_[j], labeled_rna_[j], predicted_int);
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
      predict_interaction(fm, labeled_aa_[j], labeled_rna_[j], predicted_int);
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
semisupervised_training(FeatureManager& fm)
{
  VU idx(labeled_aa_.size());
  std::iota(std::begin(idx), std::end(idx), 0);
  semisupervised_training(fm, idx);
}

void
PRactIP::
semisupervised_training(FeatureManager& fm, const VU& use_idx)
{
  // initial supervised learning
  VU idx(use_idx);
  for (uint t=0; t!=d_max_; ++t) {
    float total_loss=0.0;
    std::random_shuffle(std::begin(idx), std::end(idx)); // shuffle the order for each round
    for (auto i : idx) {
      total_loss += supervised_training(fm, labeled_aa_[i], labeled_rna_[i], labeled_int_[i]);
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
      std::random_shuffle(idx.begin(), idx.end()); // shuffle the order for each round

      float unlabeled_loss=0.0;
      float total_loss=0.0;
      for (auto i : idx)
      {
        if (i<labeled_aa_.size()) // labeled data
        {
          //std::cout << ">> supervised" << std::endl;
          total_loss += supervised_training(fm, labeled_aa_[i], labeled_rna_[i], labeled_int_[i]);
        }
        else                    // unlabeled data
        {
          const auto iu=i-labeled_aa_.size();
          // common structure prediction
          VVVU common_int;
          predict_common_interaction(fm, unlabeled_aa_[iu], unlabeled_rna_[iu], common_int);
          assert(common_int.size()==unlabeled_aa_[iu].num_sequences());
          assert(common_int.size()==unlabeled_rna_[iu].num_sequences());
          // train from predicted common structures
          for (uint k=0; k!=common_int.size(); ++k)
          {
            //std::cout << ">> semi-supervised" << std::endl;
            unlabeled_loss += supervised_training(fm, unlabeled_aa_[iu].seq(k), unlabeled_rna_[iu].seq(k), 
                                                  common_int[k], false, mu_);
          }
        }
      }
    }
  }

  fm.regularization_fobos();
}

float
PRactIP::
predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna, VVU& predicted_int) const
{
  VVF int_weight;
  VF aa_weight, rna_weight;
  fm.calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  return predict_interaction(fm, aa, rna, int_weight, aa_weight, rna_weight, predicted_int);
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
predict_interaction_object(const AA& aa, const RNA& rna,
                           const VVF& int_weight, const VF& aa_weight, const VF& rna_weight,
                           VI& x, VI& y, VVI& z, VI& sl_x, VI& sl_y, IP& ip) const
{
  const uint aa_len = x.size();
  const uint rna_len = y.size();
  const int MAX_INTERACTION = 20;

  for (uint i=0; i!=aa_len; ++i)
    if (aa_weight[i]>0.0) {
      x[i] = ip.make_variable(aa_weight[i]);
      sl_x[i] = ip.make_variable(-exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint j=0; j!=rna_len; ++j)
    if (rna_weight[j]>0.0) {
      y[j] = ip.make_variable(rna_weight[j]);
      sl_y[j] = ip.make_variable(-exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (int_weight[i][j]>0.0)
        z[i][j] = ip.make_variable(int_weight[i][j]);
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
predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna,
                    const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p) const
{
  const uint aa_len = int_weight.size();
  const uint rna_len = int_weight[0].size();
  
  IP ip(IP::MAX, n_th_);

  VI x(aa_len, -1);             // binding site in AA
  VI y(rna_len, -1);            // binding site in RNA
  VVI z(aa_len, VI(rna_len, -1)); // interactions
  VI sl_x(aa_len, -1);            // slack variables for AA to relax some constraints
  VI sl_y(rna_len, -1);           // slack variables for RNA to relax some constraints

  predict_interaction_object(aa, rna, int_weight, aa_weight, rna_weight, x, y, z, sl_x, sl_y, ip);
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
predict_common_interaction(const FeatureManager& fm, const Alignment<AA>& aa, const Alignment<RNA>& rna, VVVU& predicted_int)
{
  const uint n_seq=aa.num_sequences();
  assert(n_seq==rna.num_sequences());

  VVVF int_weight(n_seq);
  VVF aa_weight(n_seq), rna_weight(n_seq);
  for (uint i=0; i!=n_seq; ++i)
    fm.calculate_feature_weight(aa.seq(i), rna.seq(i), int_weight[i], aa_weight[i], rna_weight[i]);
  return predict_common_interaction(fm, aa, rna, int_weight, aa_weight, rna_weight, predicted_int);
}

float
PRactIP::
predict_common_interaction(const FeatureManager& fm, const Alignment<AA>& aa, const Alignment<RNA>& rna, 
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

    predict_interaction_object(aa.seq(k), rna.seq(k), int_weight[k], aa_weight[k], rna_weight[k], 
                               x[k], y[k], z[k], sl_x[k], sl_y[k], ip);
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
      FeatureManager fm(lambda_, eta0_);
      if (unlabeled_aa_.size()>0)
        semisupervised_training(fm);
      else
        supervised_training(fm);
      fm.store_parameters(param_file_.c_str());
    }
  }
  else
  {
    if (args_.size()<4)
    {
      cmdline_parser_print_help();
      return 0;
    }
    FeatureManager fm(lambda_, eta0_);
    if (param_file_.empty())
      fm.default_parameters();
    else
      fm.restore_parameters(param_file_.c_str());

    AA aa(args_[0], args_[1]); 
    RNA rna(args_[2], args_[3]);
    VVU predicted_int;
    float s=predict_interaction(fm, aa, rna, predicted_int);
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
