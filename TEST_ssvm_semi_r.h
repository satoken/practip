#ifndef _TEST_SSVM_SEMI_R_H_
#define _TEST_SSVM_SEMI_R_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string.h>
#include <getopt.h>
#include <stdlib.h>
#include <algorithm>
#include <ctime>

using namespace std;

const int MAX_len = 4000;
const double INF = 5000;

const char base[4] = {'A','C','G','U'};
const char aa[20] = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};


  double PPV_edge = 0, PPV_Pro = 0, PPV_RNA = 0;
  double SEN_edge = 0, SEN_Pro = 0, SEN_RNA = 0;
  double SPE_edge = 0, SPE_Pro = 0, SPE_RNA = 0;

class ssvm
{
 public:
  int Pro_len;
  char* protein;
  int RNA_len;
  char* RNA;
  double** cost;
  int** correct_match;

  double learn_num;

  double H_bonus;       // Hamming bonus
  double f_clip;          // fobos clip
  double eta;         // weight bonus
  double rate;          // weight threshold rate
  int output;          // 0:print only result 2:output_feature

  int result;
  int unlabel; // unlabeled
  double W;
  int predict_edge;
  double un_predict_edge; // unlabeled
  int all_correct_edge;
  int correct_predict;
  int correct_predict_protein;
  int* correct_protein;
  int correct_predict_rna;
  int* correct_rna;
  int all_edge; int all_protein; int all_RNA;
  
  int f_Rss5;/////
  int f_Pss5;
  int f_PsRs;
  int f_Ps3Rs3;
  int f_AAp;
  int f_AAab;
  int f_AAh;
  int f_Rpp;
  int f_P1;
  int f_P1R1;
  int f_P3R1;
  int f_P3R3;
  int f_P5R1;
  int f_P1R5;
  int f_P5R5;
  int out_train_weight;
  int out_cost;

  // RNA secondary structure : centroidfold
  char* rna_ss; // '.':unpair, '|':pair, '-':seqend
  // 5 - mer
  char* cRss5;
  char* cRss5_r;
  map<string, int> mcorrect_Rss5;
  map<string, int> mRss5;
  map<string, int> u_mRss5; // unlabeled
  map<string, double> mRss5_w;
  map<string, double> u_mRss5_w; // unlabeled
  map<string, double> u_rate_mRss5; // unlabeled
 
  // Protein secondary structure : PDB data
  char* pro_ss; //  '=':helix , '>':sheet ,  '.':other
  // 5 - mer
  char* cPss5;
  char* cPss5_r;
  map<string, int> mcorrect_Pss5;
  map<string, int> mPss5;
  map<string, int> u_mPss5; // unlabeled
  map<string, double> mPss5_w;
  map<string, double> u_mPss5_w; // unlabeled
  map<string, double> u_rate_mPss5; // unlabeled  

  // Pro-RNA secondary
  char* cPsRs;
  map<string, int> mcorrect_PsRs;
  map<string, int> mPsRs;
  map<string, int> u_mPsRs; // unlabeled
  map<string, double> mPsRs_w;
  map<string, double> u_mPsRs_w; // unlabeled
  map<string, double> u_rate_mPsRs; // unlabeled 

  // Pro-RNA secondary 3-3
  char* cPs3Rs3;
  char* cPs3Rs3_r;
  map<string, int> mcorrect_Ps3Rs3;
  map<string, int> mPs3Rs3;
  map<string, int> u_mPs3Rs3; // unlabeled
  map<string, double> mPs3Rs3_w;
  map<string, double> u_mPs3Rs3_w; // unlabeled
  map<string, double> u_rate_mPs3Rs3; // unlabeled  

  // AA property  Muppirala,2011  AGV, ILFP, YMTS, HNQW, RK, DE, C
  int* mcorrect_AAp;
  int* mAAp;
  map<int, double> u_mAAp;  // unlabeled
  map<int, double> mAAp_w;
  map<int, double> u_mAAp_w; // ubnlabeled
  map<int, double> u_rate_mAAp; // unlabeled
  
  // AA acid or base or neutral   DE RHK
  int* mcorrect_AAab;
  int* mAAab;
  map<int, double> u_mAAab; // unlabeled
  map<int, double> mAAab_w;
  map<int, double> u_mAAab_w; // unlabeled
  map<int, double> u_rate_mAAab; // unlabeled
  
  // AA hydropathy 3  MFLVIAC:phobic, RKQNEDH:philic
  int* mcorrect_AAh;
  int* mAAh;
  map<int, double> u_mAAh; // unlabeled
  map<int, double> mAAh_w;
  map<int, double> u_mAAh_w; // unlabeled
  map<int, double> u_rate_mAAh; // unlabeled

  // RNA purine(AG) or pyrimidine(CU)
  int* mcorrect_Rpp;
  int* mRpp;
  map<int, double> u_mRpp; // unlabeled
  map<int, double> mRpp_w;
  map<int, double> u_mRpp_w; // unlabeled
  map<int, double> u_rate_mRpp; // unlabeled
  
  // Protein 1 -
  map<char, int> mcorrect_P1;
  map<char, double> mP1;
  map<char, double> u_mP1; // unlabeled
  map<char, double> mP1_w;
  map<char, double> u_mP1_w; // unlabeled
  map<char, double> u_rate_mP1; // unlabeled
  
  // Protein 1 - RNA 1
  char* cP1R1;
  map<string, int> mcorrect_P1R1;
  map<string, int> mP1R1;
  map<string, double> u_mP1R1; // unlabeled
  map<string, double> mP1R1_w;
  map<string, double> u_mP1R1_w; // unlabeled
  map<string, double> u_rate_mP1R1; // unlabeled
  
  // Protein 3 - RNA 1
  char* cP3R1;
  char* cP3R1_r;
  map<string, int> mcorrect_P3R1;
  map<string, int> mP3R1;
  map<string, double> u_mP3R1; // unlabeled
  map<string, double> mP3R1_w;
  map<string, double> u_mP3R1_w; // unlabeled
  map<string, double> u_rate_mP3R1; // unlabeled
  
  // Protein 3 - RNA 3
  char* cP3R3;
  char* cP3R3_r;
  map<string, int> mcorrect_P3R3;
  map<string, int> mP3R3;
  map<string, double> u_mP3R3; // unlabeled
  map<string, double> mP3R3_w;
  map<string, double> u_mP3R3_w; // unlabeled
  map<string, double> u_rate_mP3R3; // unlabeled
  
  // Protein 5 - RNA 1
  char* cP5R1;
  char* cP5R1_r;
  map<string, int> mcorrect_P5R1;
  map<string, int> mP5R1;
  map<string, double> u_mP5R1; // unlabeled
  map<string, double> mP5R1_w;
  map<string, double> u_mP5R1_w; // unlabeled
  map<string, double> u_rate_mP5R1; // unlabeled
  
  // Protein 1 - RNA 5
  char* cP1R5;
  char* cP1R5_r;
  map<string, int> mcorrect_P1R5;
  map<string, int> mP1R5;
  map<string, double> u_mP1R5; // unlabeled
  map<string, double> mP1R5_w;
  map<string, double> u_mP1R5_w; // unlabeled
  map<string, double> u_rate_mP1R5; // unlabeled

  // Protein 5 - RNA 5
  char* cP5R5;
  char* cP5R5_r;
  map<string, int> mcorrect_P5R5;
  map<string, int> mP5R5;
  map<string, double> u_mP5R5; // unlabeled
  map<string, double> mP5R5_w;
  map<string, double> u_mP5R5_w; // unlabeled
  map<string, double> u_rate_mP5R5; // unlabeled

  ssvm(vector<string>, vector<string>, vector<string>, vector<string>, vector<string>, vector<string>, vector<string>, vector<string>, int, double, double, double, double, int); // training
  
  ~ssvm(){
    for(int j = 0; j < MAX_len*2; j++){
      delete[] cost[j];
    }delete[] cost;
    
    for(int j = 0; j < MAX_len*2; j++){
      delete[] correct_match[j];
    }delete[] correct_match;
    
    delete[] protein; delete[] RNA;
    delete[] rna_ss; delete[] pro_ss;

    delete[] cRss5; delete[] cRss5_r;
    delete[] cPss5; delete[] cPss5_r;
    delete[] cPsRs;
    delete[] cPs3Rs3; delete[] cPs3Rs3_r;
    delete[] mcorrect_AAp; delete[] mAAp;
    delete[] mcorrect_AAab; delete[] mAAab;
    delete[] mcorrect_AAh; delete[] mAAh;
    delete[] mcorrect_Rpp; delete[] mRpp; 
    delete[] cP1R1;
    delete[] cP3R1; delete[] cP3R1_r;
    delete[] cP3R3; delete[] cP3R3_r;
    delete[] cP5R1; delete[] cP5R1_r;
    delete[] cP1R5; delete[] cP1R5_r;
    delete[] cP5R5; delete[] cP5R5_r;
  };

  void label_learn(string, string, string);
  void unlabel_learn(string, string);
  void semi_learn(string, string, string);
  void predict(string, string, string);

  void new_array();
  void memset_seq();
  void set_parameter(double, double, double, double);
  void set_accuracy();

  void init_weight();
  void init_graph_parameter();
  void init_match_parameter();

  void update_feature_parameter(int);//0:label, 1:label+unlabel
  void unlabel_update();
  void fobos();
  void unlabel_fobos();
  double clip(double);

  void update_cost();
  void unlabeled_count();
  void correct_matching(string);
  void Hamming_distance();
  void add_edge(int, int, double);
  int BellmanFord(int, int);
  int max_match(int ,int);
  void bm(string, string);

  int get_pro_seq(string);
  int get_rna_seq(string);
  int get_pro_2nd_pdb(string);
  int get_pro_2nd_ss2(string);
  int AA_7group(char);
  int AA_3ab(char);
  int AA_3hydro(char);
  int RNA_2pp(char);

  void output_result();
  void output_train_weight();
  void output_cost();
};

#endif
