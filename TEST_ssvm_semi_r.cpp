#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string.h>
#include <getopt.h>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <ctime>
#include <cassert>


#include "ip.h"
#include "TEST_ssvm_semi_r.h"

float
PRactIP::
label_learn(const AA& aa, const RNA& rna, const VVU& correct_edges)
{
  VVF edge_weight;
  calculate_edge_weight(aa, rna, edge_weight);
  penalize_correct_matching(edge_weight, correct_edges);
  VVU predicted_edges;
  predict_matching(edge_weight, predicted_edges);
  update_feature_weight(aa, rna, predicted_edges, correct_edges);
  regularization_fobos();
  return 0.0;
}

void
PRactIP::
read_correct_matching(const std::string& filename, VVU& correct_edges) const
{
  uint r, a;
  std::ifstream is(filename.c_str());
  while (is >> r >> a)
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
  EdgeWeightCalculator(const std::vector<std::map<std::string,float> >& feature_weight,
                       VVF& edge_weight)
    : feature_weight_(feature_weight),
      edge_weight_(edge_weight)
  { }

  void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    std::map<std::string,float>::const_iterator m;
    m=feature_weight_[fgroup].find(fname);
    if (m!=feature_weight_[fgroup].end())
      edge_weight_[i][j] += m->second;
  }

  const std::vector<std::map<std::string,float> >& feature_weight_;
  VVF& edge_weight_;
};

void
PRactIP::
calculate_edge_weight(const AA& aa, const RNA& rna, VVF& edge_weight) const
{
  EdgeWeightCalculator e(feature_weight_, edge_weight);
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
      edge_weight[i][j] -= negative_penalty_;

  for (uint i=0; i!=correct_edges.size(); ++i) 
    for (VU::const_iterator j=correct_edges[i].begin(); j!=correct_edges[i].end(); ++j)
      edge_weight[i][*j] += positive_penalty_+negative_penalty_;
}

struct FeatureWeightUpdater
{
  FeatureWeightUpdater(std::vector<std::map<std::string,float> >& feature_weight, float eta)
    : feature_weight_(feature_weight), eta_(eta)
  { }

  void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    feature_weight_[fgroup].insert(std::make_pair(std::string(fname),0.0f)).first->second += eta_;
  }

  std::vector<std::map<std::string,float> >& feature_weight_;
  float eta_;
};

void
PRactIP::
update_feature_weight(const AA& aa, const RNA& rna, const VVU& predicted_edges, const VVU& correct_edges)
{
  FeatureWeightUpdater f(feature_weight_, eta_);
  for (uint i=0; i!=correct_edges.size(); ++i) {
    for (VU::const_iterator j=correct_edges[i].begin(); j!=correct_edges[i].end(); ++j) {
      extract_feature(aa, rna, i, *j, f);
    }
  }

  FeatureWeightUpdater g(feature_weight_, -eta_);
  for (uint i=0; i!=predicted_edges.size(); ++i) {
    for (VU::const_iterator j=predicted_edges[i].begin(); j!=predicted_edges[i].end(); ++j) {
      extract_feature(aa, rna, i, *j, g);
    }
  }
}

struct CountFeature
{
  CountFeature(std::vector<std::map<std::string,uint> >& fc, VU& tc)
    : fc_(fc), tc_(tc)
  { }

  void operator()(uint fgroup, const char* fname, uint i, uint j)
  {
    fc_[fgroup].insert(std::make_pair(std::string(fname),0u)).first->second++;
    tc_[fgroup]++;
  }

  std::vector<std::map<std::string,uint> >& fc_;
  VU& tc_;
};

void
PRactIP::
count_feature(const AA& aa, const RNA& rna, const VVU& edge, std::vector<std::map<std::string,uint> >& fc, VU& tc) const
{
  CountFeature c(fc, tc);
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

void
PRactIP::
regularization_fobos()
{
  std::vector<std::map<std::string,float> >::iterator i;
  for (i=feature_weight_.begin(); i!=feature_weight_.end(); ++i)
  {
    std::map<std::string,float>::iterator j;
    for (j=i->begin(); j!=i->end(); ++j)
      j->second = clip(j->second, f_clip_);
  }
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
      if (edge_weight[i][j]>=0.0)
        x[i][j] = ip.make_variable(edge_weight[i][j]);
  
  ip.update();

  for (uint i=0; i!=aa_len; ++i) {
    int row = ip.make_constraint(IP::DB, 0, 1);
    for (uint j=0; j!=rna_len; ++j) 
      if (x[i][j]>=0)
        ip.add_constraint(row, x[i][j], 1);
  }

  for (uint j=0; j!=rna_len; ++j) {
    int row = ip.make_constraint(IP::DB, 0, 1);
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

  fp = fopen((filename+".ss2").c_str(), "r");
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

void
printUsage()
{
  std::cout << "Usage: [-option] labeled_data [unlabeled_data]" << std::endl
            << "\t-b : Hamming bonus \tdefault = 0.01" << std::endl
            << "\t-e : learning eta \tdefault = 0.01" << std::endl
            << "\t-f : FOBOS clip \tdefault= 0.001" << std::endl
            << "\t-r : threshold rate \tdefault = 0.9" << std::endl
            << "\t-c : cross-validation \tdefault = 4" << std::endl
            << "-o : output only result" << std::endl;
}

#if 0
int
main(int argc, char** argv)
{
  int label_num = 0;
  int unlabel_num = 0;
  string Pseq, Rseq, Match;
  vector<string> pseq, rseq, match;
  vector<string> p_tr, r_tr, m_tr, p_te, r_te, m_te;
  vector<string> u_p, u_r;
  double H_bonus_ = 0.01;
  double eta_ = 0.01;
  double f_clip_ = 0.001;
  double rate_ = 0.9;
  int CV = 4;
  int output_ = 1;

  int input = 0;
  while((input = getopt(argc, argv, "hotb:e:f:r:c:")) != -1) {
    switch(input) {
    case 'b':
      H_bonus_ = atof(optarg);
      cout << "H_bonus = " << H_bonus_ << endl;
      break;
    case 'e':
      eta_ = atof(optarg);
      cout << "eta = " << eta_ << endl;
      break;      
    case 'f':
      f_clip_ = atof(optarg);
      cout << "f_clip = " << f_clip_ << endl;
      break;      
    case 'r':
      rate_ = atof(optarg);
      cout << "rate = " << rate_ << endl;
      break;
    case 'c':
      CV = atof(optarg);
      cout << "CV = " << CV << endl;
      break;
    case 'o':
      output_ = 0;
      break;
    case 'h':   
    case ':':
    case '?':
      printUsage();
      return 0;
    }
  }
  if(!argv[optind]) {
    printUsage();
    return 0;
  }

  ifstream labeled_data(argv[optind]);
  ifstream unlabeled_data(argv[optind+1]);

  if(labeled_data.fail()) {
    cerr << "Labeled_Data File do not exist.\n";
    return 0;
  }
  while(!labeled_data.eof()) {
    labeled_data >> Pseq >> Rseq >> Match;
    if(output_>0){ cout << Pseq << " " << Rseq << " " << Match << endl;}
    pseq.push_back(Pseq);
    rseq.push_back(Rseq);
    match.push_back(Match);
    label_num++;
  }
  
  if(argv[optind+1] != NULL){
    if(unlabeled_data.fail()) {
      cerr << "Unlabeled_Data File do not exist.\n";
      return 0;
    }
    if(output_>0){ cout << "--- unlabeled data ---" << endl;}
    while(!unlabeled_data.eof()) {
      unlabeled_data >> Pseq >> Rseq;
      if(output_>0){ cout << " " << Pseq << " " << Rseq << " " << endl;}
      u_p.push_back(Pseq);
      u_r.push_back(Rseq);
      unlabel_num++;
    }
  }

  srand( unsigned (time(0)));

  if(CV == 1){
    ssvm train(pseq, rseq, match, pseq, rseq, match, u_p, u_r, unlabel_num, H_bonus_, f_clip_, eta_, rate_, output_);
  }

  if(CV != 1) {
    for (int j=0; j<CV; j++) {
      for(int i=0; i<label_num; i++) {
	if(i % CV == j) {
	  p_te.push_back(pseq[i]);
	  r_te.push_back(rseq[i]);
	  m_te.push_back(match[i]);
	}
	else {
	  p_tr.push_back(pseq[i]);
	  r_tr.push_back(rseq[i]);
	  m_tr.push_back(match[i]);
	}
      }
      cout << "CV " << j+1 << "-" << CV << endl;
      ssvm train(p_tr, r_tr, m_tr, p_te, r_te, m_te, u_p, u_r, unlabel_num, H_bonus_, f_clip_, eta_, rate_, output_);
      p_te.clear(); r_te.clear(); m_te.clear();
      p_tr.clear(); r_tr.clear(); m_tr.clear();
    }
  }

  cout << endl;
  cout << "PPV_edge = " << PPV_edge/CV << endl;
  cout << "PPV_Pro = " << PPV_Pro/CV << endl;
  cout << "PPV_RNA = " << PPV_RNA/CV << endl;
  cout << "SEN_edge = " << SEN_edge/CV << endl;
  cout << "SEN_Pro = " << SEN_Pro/CV << endl;
  cout << "SEN_RNA = " << SEN_RNA/CV << endl;
  cout << "SPE_edge = " << SPE_edge/CV << endl;
  cout << "SPE_Pro = " << SPE_Pro/CV << endl;
  cout << "SPE_RNA = " << SPE_RNA/CV << endl;

  return 0;
}
#endif
