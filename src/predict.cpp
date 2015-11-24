#include <cassert>
#include "ip.h"
#include "feature.h"
#include "aa.h"
#include "rna.h"
#include "alignment.h"
#include "predict.h"

float
Predictor::
predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna, VVU& predicted_int) const
{
  VVF int_weight;
  VF aa_weight, rna_weight;
  fm.calculate_feature_weight(aa, rna, int_weight, aa_weight, rna_weight);
  return predict_interaction(fm, aa, rna, int_weight, aa_weight, rna_weight, predicted_int, 1.0);
}

void
Predictor::
predict_interaction_object(const AA& aa, const RNA& rna,
                           const VVF& int_weight, const VF& aa_weight, const VF& rna_weight,
                           VI& x, VI& y, VVI& z, VI& sl_x, VI& sl_y, IP& ip, float mu) const
{
  const uint aa_len = x.size();
  const uint rna_len = y.size();
  const int MAX_INTERACTION = 20;

  for (uint i=0; i!=aa_len; ++i)
    if (aa_weight[i]>0.0) {
      x[i] = ip.make_variable(mu*aa_weight[i]);
      sl_x[i] = ip.make_variable(-mu*exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint j=0; j!=rna_len; ++j)
    if (rna_weight[j]>0.0) {
      y[j] = ip.make_variable(mu*rna_weight[j]);
      sl_y[j] = ip.make_variable(-mu*exceed_penalty_, 0, MAX_INTERACTION);
    }
  for (uint i=0; i!=aa_len; ++i)
    for (uint j=0; j!=rna_len; ++j)
      if (int_weight[i][j]>0.0)
        z[i][j] = ip.make_variable(mu*int_weight[i][j]);
}

void
Predictor::
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
Predictor::
predict_interaction(const FeatureManager& fm, const AA& aa, const RNA& rna,
                    const VVF& int_weight, const VF& aa_weight, const VF& rna_weight, VVU& p, float mu) const
{
  const uint aa_len = int_weight.size();
  const uint rna_len = int_weight[0].size();
  
  IP ip(IP::MAX, n_th_);

  VI x(aa_len, -1);             // binding site in AA
  VI y(rna_len, -1);            // binding site in RNA
  VVI z(aa_len, VI(rna_len, -1)); // interactions
  VI sl_x(aa_len, -1);            // slack variables for AA to relax some constraints
  VI sl_y(rna_len, -1);           // slack variables for RNA to relax some constraints

  predict_interaction_object(aa, rna, int_weight, aa_weight, rna_weight, x, y, z, sl_x, sl_y, ip, mu);
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
Predictor::
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
Predictor::
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
                               x[k], y[k], z[k], sl_x[k], sl_y[k], ip, mu_);
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

float
Predictor::
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


