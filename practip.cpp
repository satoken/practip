#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
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
                   const VVU& predicted_edges, const VVU& correct_edges,
                   Accuracy& edge_acc, Accuracy& aa_acc, Accuracy& rna_acc)
{
  std::vector<bool> c_edge(aa.seq.size()*rna.seq.size(), false);
  std::vector<bool> p_edge(aa.seq.size()*rna.seq.size(), false);
  std::vector<bool> c_aa(aa.seq.size(), false);
  std::vector<bool> p_aa(aa.seq.size(), false);
  std::vector<bool> c_rna(rna.seq.size(), false);
  std::vector<bool> p_rna(rna.seq.size(), false);
  for (uint i=0; i!=correct_edges.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_edges[i])
      c_edge[rna.seq.size()*i+*j]=c_aa[i]=c_rna[*j]=true;
  for (uint i=0; i!=predicted_edges.size(); ++i)
    FOREACH (VU::const_iterator, j, predicted_edges[i])
      p_edge[rna.seq.size()*i+*j]=p_aa[i]=p_rna[*j]=true;
  edge_acc.calculate(c_edge, p_edge);
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
    labeled_matching_.push_back(VVU(aa_len));
    read_correct_matching(matching, labeled_matching_.back());
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
    float eta=eta0_/std::sqrt(t+1);
    float total_loss=0.0;
    std::random_shuffle(idx.begin(), idx.end());
    FOREACH (VU::const_iterator, it, idx) {
      total_loss += supervised_training(labeled_aa_[*it], labeled_rna_[*it], labeled_matching_[*it], eta);
    }
  }
}

float
PRactIP::
supervised_training(const AA& aa, const RNA& rna, const VVU& correct_edges, float eta)
{
  float loss=0.0;
  VVF edge_weight;
  calculate_edge_weight(aa, rna, edge_weight);
  for (uint i=0; i!=correct_edges.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_edges[i])
      loss -= edge_weight[i][*j];
  penalize_correct_matching(edge_weight, correct_edges);
  VVU predicted_edges;
  loss += predict_matching(edge_weight, predicted_edges);
#if 0
  for (uint i=0; i!=correct_edges.size(); ++i) {
    if (!correct_edges[i].empty() || !predicted_edges[i].empty())
    {
      std::cout << i << ": [ ";
      FOREACH (VU::const_iterator, j, correct_edges[i])
        std::cout << *j << "(" << edge_weight[i][*j] << ") ";
      std::cout << "], [ ";
      FOREACH (VU::const_iterator, j, predicted_edges[i])
        std::cout << *j << "(" << edge_weight[i][*j] << ") ";
      std::cout << "]" << std::endl;
    }
  }
  std::cout << std::endl;
#endif  
  update_feature_weight(aa, rna, predicted_edges, correct_edges, eta);
  loss += regularization_fobos(eta);
  return loss;
}

void
PRactIP::
cross_validation(uint n, void (PRactIP::*training)(const VU&))
{
  AccuracySummary edge_summary, aa_summary, rna_summary;
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
      VVU predicted_edges;
      predict_matching(labeled_aa_[*j], labeled_rna_[*j], predicted_edges);
      Accuracy edge_acc, aa_acc, rna_acc;
      calculate_accuracy(labeled_aa_[*j], labeled_rna_[*j],
                         predicted_edges, labeled_matching_[*j],
                         edge_acc, aa_acc, rna_acc);
      edge_summary.add(edge_acc);
      aa_summary.add(aa_acc);
      rna_summary.add(rna_acc);
    }
  }
  edge_summary.summary();
  aa_summary.summary();
  rna_summary.summary();
  std::cout << "Edge: "
            << "PPV=" << edge_summary.ppv_avg << "(" << edge_summary.ppv_sd << ") "
            << "SEN=" << edge_summary.sen_avg << "(" << edge_summary.sen_sd << ") "
            << "F=" << edge_summary.fval_avg << "(" << edge_summary.fval_sd << ") "
            << std::endl;
  std::cout << "AA: "
            << "PPV=" << aa_summary.ppv_avg << "(" << aa_summary.ppv_sd << ") "
            << "SEN=" << aa_summary.sen_avg << "(" << aa_summary.sen_sd << ") "
            << "F=" << aa_summary.fval_avg << "(" << aa_summary.fval_sd << ") "
            << std::endl;
  std::cout << "RNA: "
            << "PPV=" << rna_summary.ppv_avg << "(" << rna_summary.ppv_sd << ") "
            << "SEN=" << rna_summary.sen_avg << "(" << rna_summary.sen_sd << ") "
            << "F=" << rna_summary.fval_avg << "(" << rna_summary.fval_sd << ") "
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
      total_loss += supervised_training(labeled_aa_[*it], labeled_rna_[*it], labeled_matching_[*it], eta);
    }
  }    
}

float
PRactIP::
unsupervised_training(const AA& aa, const RNA& rna,
                      std::vector<FC>& fc, VU& tc) const
{
  VVU p;
  float score = predict_matching(aa, rna, p);
  count_feature(aa, rna, p, fc, tc);
  return score;
}

float
PRactIP::
predict_matching(const AA& aa, const RNA& rna, VVU& predicted_edges) const
{
  VVF edge_weight;
  calculate_edge_weight(aa, rna, edge_weight);
  return predict_matching(edge_weight, predicted_edges);
}

void
PRactIP::
read_correct_matching(const std::string& filename, VVU& correct_edges) const
{
  uint r, a;
  std::ifstream is(filename.c_str());
  while (is >> a >> r)
    correct_edges[a].push_back(r);
}

template < class Func >
void
PRactIP::
extract_feature(const AA& aa, const RNA& rna, uint i, uint j, Func& func) const
{
  const uint rna_len=rna.seq.size();
  const uint aa_len=aa.seq.size();
  char buf[20];
  if (use_feature_[FG_Rss5]) {
    sprintf(buf, "%c%c%c%c%c",
            j>=2 ? rna.ss[j-2] : '-',
            j>=1 ? rna.ss[j-1] : '-',
            rna.ss[j],
            j+1<rna_len ? rna.ss[j+1] : '-',
            j+2<rna_len ? rna.ss[j+2] : '-');
    func(FG_Rss5, buf, i, j);
  }
  // Proteinss 
  if (use_feature_[FG_Pss5]) {
    sprintf(buf, "%c%c%c%c%c",
            i>=2 ? aa.ss[i-2] : '-',
            i>=1 ? aa.ss[i-1] : '-',
            aa.ss[i],
            i+1<aa_len ? aa.ss[i+1] : '-',
            i+2<aa_len ? aa.ss[i+2] : '-');
    func(FG_Pss5, buf, i, j);
  }
  // Pro-RNAss 
  if (use_feature_[FG_PsRs]) {
    sprintf(buf, "%c%c", aa.ss[i], rna.ss[j]);
    func(FG_PsRs, buf, i, j);
  }
  // Pro-RNA ss 3-3
  if (use_feature_[FG_Ps3Rs3]) {
    sprintf(buf, "%c%c%c%c%c%c",
            i>=1 ? aa.ss[i-1] : '-',
            aa.ss[i],
            i+1<=aa_len ? aa.ss[i+1] :'-',
            j>=1 ? rna.ss[j-1] : '-',
            rna.ss[j],
            j+1<=rna_len ? rna.ss[j+1] : '-');
    func(FG_Ps3Rs3, buf, i, j);
  }
  // AAp
  if (use_feature_[FG_AAp]) {
    sprintf(buf, "%d", AA::group7(aa.seq[i]));
    func(FG_AAp, buf, i, j);
  }
  // AAab
  if (use_feature_[FG_AAab]) {
    sprintf(buf, "%d", AA::ab3(aa.seq[i]));
    func(FG_AAab, buf, i, j);
  }
  // AAh
  if (use_feature_[FG_AAh]) {
    sprintf(buf, "%d", AA::hydro3(aa.seq[i]));
    func(FG_AAh, buf, i, j);
  }
  // Rpp
  if (use_feature_[FG_Rpp]) {
    sprintf(buf, "%d", RNA::pp2(rna.seq[j]));
    func(FG_Rpp, buf, i, j);
  }
  // Protein 1 -
  if (use_feature_[FG_P1]) {
    sprintf(buf, "%c", aa.seq[i]);
    func(FG_P1, buf, i, j);
  }      
  // Protein 1 - RNA 1
  if (use_feature_[FG_P1R1]) {
    sprintf(buf, "%c%c", aa.seq[i], rna.seq[j]);
    func(FG_P1R1, buf, i, j);
  }      
  // Protein 3 - RNA 1
  if (use_feature_[FG_P3R1]) {
    sprintf(buf, "%c%c%c%c",
            i>=1 ? aa.seq[i-1] : '-',
            aa.seq[i],
            i+1<aa_len ? aa.seq[i+1] : '-',
            rna.seq[j]);
    func(FG_P3R1, buf, i, j);
  }      
  // Protein 3 - RNA 3
  if (use_feature_[FG_P3R3]) {
    sprintf(buf, "%c%c%c%c%c%c",
            i>=1 ? aa.seq[i-1] : '-',
            aa.seq[i],
            i+1<aa_len ? aa.seq[i+1] : '-',
            j>=1 ? rna.seq[j-1] : '-',
            rna.seq[j],
            j+1<rna_len ? rna.seq[j+1] : '-');
    func(FG_P3R3, buf, i, j);
  }      
  // Protein 5 - RNA 1
  if (use_feature_[FG_P5R1]) {
    sprintf(buf, "%c%c%c%c%c%c",
            i>=2 ? aa.seq[i-2] : '-',
            i>=1 ? aa.seq[i-1] : '-',
            aa.seq[i],
            i+1<aa_len ? aa.seq[i+1] : '-',
            i+2<aa_len ? aa.seq[i+2] : '-',
            rna.seq[j]);
    func(FG_P5R1, buf, i, j);
  }      
  // Protein 1 - RNA 5
  if (use_feature_[FG_P1R5]) {
    sprintf(buf, "%c%c%c%c%c%c",
            aa.seq[i],
            j>=2 ? rna.seq[j-2] : '-',
            j>=1 ? rna.seq[j-1] : '-',
            rna.seq[j],
            j+1<rna_len ? rna.seq[j+1] : '-',
            j+2<rna_len ? rna.seq[j+2] : '-');
    func(FG_P1R5, buf, i, j);
  }      
  // Protein 5 - RNA 5
  if (use_feature_[FG_P5R5]) {
    sprintf(buf, "%c%c%c%c%c%c%c%c%c%c",
            i>=2 ? aa.seq[i-2] : '-',
            i>=1 ? aa.seq[i-1] : '-',
            aa.seq[i],
            i+1<aa_len ? aa.seq[i+1] : '-',
            i+2<aa_len ? aa.seq[i+2] : '-',
            j>=2 ? rna.seq[j-2] : '-',
            j>=1 ? rna.seq[j-1] : '-',
            rna.seq[j],
            j+1<rna_len ? rna.seq[j+1] : '-',
            j+2<rna_len ? rna.seq[j+2] : '-');
    func(FG_P5R5, buf, i, j);
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

void
PRactIP::
calculate_edge_weight(const AA& aa, const RNA& rna, VVF& edge_weight) const
{
  EdgeWeightCalculator e(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, edge_weight);
  edge_weight.resize(aa.seq.size());
  for (uint i=0; i!=edge_weight.size(); ++i) {
    edge_weight[i].resize(rna.seq.size());
    for (uint j=0; j!=edge_weight[i].size(); ++j) {
      edge_weight[i][j] = 0.0;
      extract_feature(aa, rna, i, j, e);
    }
  }
}

void
PRactIP::
penalize_correct_matching(VVF& edge_weight, const VVU& correct_edges) const
{
  for (uint i=0; i!=edge_weight.size(); ++i) 
    for (uint j=0; j!=edge_weight[i].size(); ++j)
      edge_weight[i][j] += neg_w_;

  for (uint i=0; i!=correct_edges.size(); ++i)
    FOREACH (VU::const_iterator, j, correct_edges[i])
      edge_weight[i][*j] -= pos_w_+neg_w_;
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

  inline void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    feature_weight_[fgroup].insert(std::make_pair(std::string(fname),0.0f)).first->second += eta_;
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
update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_edges, const VVU& correct_edges, float eta)
{
  FeatureWeightUpdater f(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, eta);
  for (uint i=0; i!=correct_edges.size(); ++i) {
    FOREACH (VU::const_iterator, j, correct_edges[i]) {
      extract_feature(aa, rna, i, *j, f);
    }
  }

  FeatureWeightUpdater g(feature_weight_, feature_count_, feature_group_weight_, feature_group_count_, -eta);
  for (uint i=0; i!=predicted_edges.size(); ++i) {
    FOREACH (VU::const_iterator, j, predicted_edges[i]) {
      extract_feature(aa, rna, i, *j, g);
    }
  }
}

struct FeatureCounter
{
  FeatureCounter(std::vector<FC >& fc, VU& tc)
    : fc_(fc), tc_(tc)
  { }

  inline void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    fc_[fgroup].insert(std::make_pair(std::string(fname),0u)).first->second++;
    tc_[fgroup]++;
  }

  std::vector<FC >& fc_;
  VU& tc_;
};

void
PRactIP::
count_feature(const AA& aa, const RNA& rna, const VVU& edge, std::vector<FC>& fc, VU& tc) const
{
  FeatureCounter c(fc, tc);
  for (uint i=0; i!=aa.seq.size(); ++i) {
    for (uint j=0; j!=rna.seq.size(); ++j) {
      if (edge[i][j]) {
        extract_feature(aa, rna, i, j, c);
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
regularization_fobos(float eta)
{
  float sum1=0.0;
  std::vector<FM >::iterator i;
  for (i=feature_weight_.begin(); i!=feature_weight_.end(); ++i)
  {
    FM::iterator j=i->begin();
    while (j!=i->end()) {
      j->second = clip(j->second, eta*lambda_);
      sum1 += std::abs(j->second);
      if (j->second==0.0)
        i->erase(j++);
      else
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
predict_matching(const VVF& edge_weight, VVU& p) const
{
  const uint aa_len = edge_weight.size();
  const uint rna_len = edge_weight[0].size();
  
  IP ip(IP::MAX, n_th_);

  VVI x(aa_len, VI(rna_len, -1));
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (edge_weight[i][j]>0.0)
        x[i][j] = ip.make_variable(edge_weight[i][j]);
  
  ip.update();

  for (uint i=0; i!=aa_len; ++i) {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint j=0; j!=rna_len; ++j) 
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  for (uint j=0; j!=rna_len; ++j) {
    int row = ip.make_constraint(IP::UP, 0, 1);
    for (uint i=0; i!=aa_len; ++i)
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  float s = ip.solve();
  
  p.resize(aa_len);
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (x[i][j]>=0 && ip.get_value(x[i][j])>0.5)
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
      case 'C': ss.push_back('.'); break;
      case 'E': ss.push_back('>'); break;
      case 'H': ss.push_back('='); break;
      default: assert(!"unreachable"); break;
    }
  }
#endif
  fclose(fp);

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

  const char* prog=getenv("CENTROID_FOLD");
  if (!prog) prog="centroid_fold";
  char cmd[1000];
  sprintf(cmd, "%s %s", prog, filename.c_str());
  fp = popen(cmd, "r");
  while (fgets(line, sizeof(line), fp)) {
    if (line[0]=='>') continue;
    if (strchr(".()[]{}<>", line[0])) {
      for (uint i=0; line[i]!='\0'; ++i) {
        switch (line[i]) {
          case '.':
            line[i]='.'; break;
          case '(': case ')':
          case '[': case ']':
          case '{': case '}':
          case '<': case '>':
            line[i]='|'; break;
          case '\n':
          case ' ':
            line[i]='\0'; break;
        }
      }
      this->ss+=line;
    }
  }

  assert(this->seq.size()==this->ss.size());
  return this->seq.size();
}

// static
int
PRactIP::AA::
group7(char a)
{
  if(a == 'A' || a == 'G' || a == 'V') {return 0;}
  if(a == 'I' || a == 'L' || a == 'F' || a == 'P') {return 1;}
  if(a == 'Y' || a == 'M' || a == 'T' || a == 'S') {return 2;}
  if(a == 'H' || a == 'N' || a == 'Q' || a == 'W') {return 3;}
  if(a == 'R' || a == 'K') {return 4;}
  if(a == 'D' || a == 'E') {return 5;}
  if(a == 'C') {return 6;}
  return -1;
}

// static
int
PRactIP::AA::
ab3(char a)
{
  if(a == 'D' || a == 'E') {return 0;}
  if(a == 'R' || a == 'H' || a == 'K') {return 1;}
  else{return 2;}
}

// static
int
PRactIP::AA::
hydro3(char a)
{
  if(a == 'A' || a == 'M' || a == 'C' || a == 'F' || a == 'L' || a == 'V' || a == 'I'){return 0;}
  if(a == 'R' || a == 'K' || a == 'Q' || a == 'N' || a == 'E' || a == 'D' || a == 'H'){return 1;}
  else{return 2;}
}

//static
int
PRactIP::RNA::
pp2(char b)
{
  if(b == 'A' || b == 'G'){return 0;}
  else{return 1;}
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
