#include <algorithm>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cstring>
#include <cerrno>
#include "feature.h"

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


FeatureManager::
FeatureManager(float lambda, float eta0)
  : feature_weight_(FG_NUM),
    use_feature_(FG_NUM, true),
    epoch_(0),
    lambda_(lambda),
    eta0_(eta0)
{}

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
FeatureManager::
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
    { -1u, nullptr, 0, nullptr, 0 }
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
FeatureManager::
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
    { -1u, nullptr, 0 }
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
FeatureManager::
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
    { -1u, nullptr, 0 }
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
FeatureManager::
store_parameters(const char* filename) const
{
  std::ofstream os(filename);
  if (!os) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  for (uint k=0; k!=feature_weight_.size(); ++k)
  {
    os << "[ " << groupname[k] << " weight ]" << std::endl;
    for (const auto& e : feature_weight_[k])
    {
      assert(e.second.last_updated==epoch_);
#if 0
      if (e.second.weight!=0.0)
        os << e.first << " " << e.second.weight << " " << e.second.sum_squared_grad << std::endl;
#else
      os << e.first << " " << e.second.weight << " " << e.second.sum_squared_grad << std::endl;
#endif
    }
    os << std::endl;
  }
}

void
FeatureManager::
restore_parameters(const char* filename)
{
  std::ifstream is(filename);
  if (!is) throw std::runtime_error(std::string(strerror(errno)) + ": " + filename);

  for (auto f : feature_weight_)
    f.clear();

  uint k=-1u;
  std::string first, second, third;
  while (is >> first >> second >> third)
  {
    if (first=="[")             // group name
    {
      //std::cout << second << std::endl;
      std::string temp;
      is >> temp; // "]"
      auto p = std::find(std::begin(groupname), std::end(groupname), second);
      if (p!=std::end(groupname))
        k = p - std::begin(groupname);
      else
        throw std::runtime_error(std::string("unknown groupname: ") + second);
    }
    else if (k!=-1u)
    {
      auto& e = feature_weight_[k][first];
      e.weight = std::atof(second.c_str());
      e.sum_squared_grad = std::atof(third.c_str());
    }
  }
}

void
FeatureManager::
default_parameters()
{
  static struct { const char* name; float value; } params[] = {
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

//mutable
void
FeatureManager::
calculate_feature_weight(const AA& aa, const RNA& rna, VVF& int_weight, VF& aa_weight, VF& rna_weight) const
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
                            int_weight[i][j] += lazy_update(fgroup, fname);
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
                         aa_weight[i] += lazy_update(fgroup, fname);
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
                          rna_weight[j] += lazy_update(fgroup, fname);
                        }
      );
  }
}

void
FeatureManager::
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
        auto fe = feature_weight_[k].insert(std::make_pair(e.first, FeatureWeight())).first;
        fe->second.weight -= w*e.second * eta0_/std::sqrt(1.0+fe->second.sum_squared_grad);
        fe->second.sum_squared_grad += w*e.second * w*e.second;
#if 0
        std::cerr << epoch << ", "
                  << k << ", "
                  << it->first << ", "
                  << it->second << ", "
                  << fe->second.weight << ", "
                  << eta0_/std::sqrt(1.0+fe->second.sum_squared_grad) << std::endl;
#endif
        //assert(fe->second.last_updated == epoch);
      }
    }
  }
  epoch_++;
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
FeatureManager::
lazy_update(FeatureWeight& f) const
{
  // lazy update for FOBOS
  if (f.last_updated<epoch_) 
  {
    const float eta = eta0_/std::sqrt(1.0+f.sum_squared_grad);
    const uint t = epoch_ - f.last_updated;
    f.weight = clip(f.weight, lambda_*eta*t);
    f.last_updated = epoch_;
  }
  return f.weight;
}

float
FeatureManager::
lazy_update(uint fgroup, const char* fname) const
{
  auto m=feature_weight_[fgroup].find(fname);
  if (m!=feature_weight_[fgroup].end())
    return lazy_update(m->second);
  return 0.0;
}

float
FeatureManager::
regularization() const
{
  // L1-norm
  float sum1=0.0;
  for (uint k=0; k!=feature_weight_.size(); ++k)
  {
    for (auto& e : feature_weight_[k])
    {
      auto w = lazy_update(e.second);
      sum1 += std::abs(w);
    }
  }
  return lambda_ * sum1;
}

