#include "TEST_ssvm_semi_r.h"

ssvm::
ssvm(vector<string> pro_tr_seq, vector<string> rna_tr_seq, vector<string> match_tr_data, vector<string> pro_te_seq, vector<string> rna_te_seq, vector<string> match_te_data, vector<string> un_pro_seq, vector<string> un_rna_seq, int unlabel_data, double H_bonus_, double f_clip_, double eta_, double rate_, int output_)
{
  output = output_;
  new_array();
  set_parameter(H_bonus_, f_clip_, eta_, rate_);
  set_accuracy();
  un_predict_edge = 0; // unlabeled

  /// labeled learning 
  unlabel = 0;
  int train_size = pro_tr_seq.size();
  int test_size = pro_te_seq.size();
  
  for(int i = 0; i < 10; i++) {
  int k = 0;
    vector<int> random;
    for (int j = 0; j < train_size; j++) random.push_back(j);
    
    random_shuffle ( random.begin(), random.end() );
    
    for (vector<int>::iterator it=random.begin(); it!=random.end(); ++it){
      
      if(output>0){ cout << "\n\n---- learning "<< i+1 << " - " << k+1 << endl;}
      label_learn(pro_tr_seq[*it], rna_tr_seq[*it], match_tr_data[*it]);
      learn_num++;
      k++;
    }
 }
  ////// unlabeled learning ////////////////////////////
  if(unlabel_data > 0) {
    set_accuracy();
    unlabel = 1; // unlabeled
    
    for(int i = 0; i < unlabel_data; i++) {
      if(output>0){ cout << "\n\n---- unlabel_learning " << i+1 << endl;}
      unlabel_learn(un_pro_seq[i], un_rna_seq[i]);
    }
    // if(output>0){ output_result();}
    
    if(predict_edge > 0){
      unlabeled_count();//unlabeled
    }
    ////// labeled learning
    set_accuracy();
    unlabel = 0; // semi

    for(int i = 0; i < train_size; i++) {
      if(output>0){ cout << "\n\n---- semi_learning " << i+1 << endl;}
      semi_learn(pro_tr_seq[i], rna_tr_seq[i], match_tr_data[i]);
      learn_num++;
    }
    // if(output>0){ output_result();}
  }
  /////////// labeled prediction /////////////////////////
  set_accuracy();
  unlabel = 2; // cut cost  
  for(int i = 0; i < test_size; i++) {
    if(output>0){ cout << "\n\n---- learned_prediction " << i+1 << endl;}
    predict(pro_te_seq[i], rna_te_seq[i], match_te_data[i]);
    if(out_cost>0){ output_cost();}
  }
  output_result();

  if(out_train_weight>0){ output_train_weight();}
}

void
ssvm::
label_learn(string p_seq, string r_seq, string match)
{
  memset_seq();
  Pro_len = get_pro_seq(p_seq);
  RNA_len = get_rna_seq(r_seq);
  init_weight();
  init_graph_parameter();
  init_match_parameter();
  update_cost();
  correct_matching(match);
  Hamming_distance();
  bm(p_seq, r_seq);
  update_feature_parameter(0);
  fobos();
}

void
ssvm::
unlabel_learn(string un_p_seq, string un_r_seq)
{
  memset_seq();
  Pro_len = get_pro_seq(un_p_seq);
  RNA_len = get_rna_seq(un_r_seq);
  init_weight();
  init_graph_parameter();
  init_match_parameter();
  update_cost();
  bm(un_p_seq, un_r_seq);
}

void
ssvm::
semi_learn(string p_seq, string r_seq, string match)
{
  memset_seq();
  Pro_len = get_pro_seq(p_seq);
  RNA_len = get_rna_seq(r_seq);
  init_weight();
  init_graph_parameter();
  init_match_parameter();
  update_cost();
  correct_matching(match);
  Hamming_distance();
  bm(p_seq, r_seq);
  update_feature_parameter(1);
  fobos();
  unlabel_fobos();
}

void
ssvm::
predict(string p_seq, string r_seq, string match)
{
  Pro_len = get_pro_seq(p_seq);
  RNA_len = get_rna_seq(r_seq);  
  
  init_graph_parameter();
  init_match_parameter();
  update_cost();
  correct_matching(match);
  bm(p_seq, r_seq);
  all_edge += Pro_len * RNA_len;
  all_protein += Pro_len;
  all_RNA += RNA_len;
}

void
ssvm::
new_array()
{
  protein = new char[MAX_len]; RNA = new char[MAX_len];
  rna_ss = new char[MAX_len]; pro_ss = new char[MAX_len];

  cost = new double*[MAX_len*2];
  for(int j = 0; j < MAX_len*2; j++){
    cost[j] = new double[MAX_len*2];}
  correct_match = new int*[MAX_len*2];
  for(int j = 0; j < MAX_len*2; j++){
    correct_match[j] = new int[MAX_len*2];}
  
  correct_protein = new int[MAX_len];  correct_rna = new int[MAX_len];

  cRss5 = new char[5]; cRss5_r = new char[5];
  cPss5 = new char[5]; cPss5_r = new char[5];
  cPsRs = new char[2];
  cPs3Rs3 = new char[6]; cPs3Rs3_r = new char[6];

  mcorrect_AAp = new int[7]; mAAp = new int[7];
  mcorrect_AAab = new int[3]; mAAab = new int[3];
  mcorrect_AAh = new int[3]; mAAh = new int[3];
  mcorrect_Rpp = new int[2]; mRpp = new int[2];

  cP1R1 = new char[2];
  cP3R1 = new char[4]; cP3R1_r = new char[4];
  cP3R3 = new char[6]; cP3R3_r = new char[6];
  cP5R1 = new char[6]; cP5R1_r = new char[6];
  cP1R5 = new char[6]; cP1R5_r = new char[6];
  cP5R5 = new char[10]; cP5R5_r = new char[10];
}

void
ssvm::
memset_seq()
{
  memset(protein, '\0', MAX_len);
  memset(RNA, '\0', MAX_len);
  memset(rna_ss, '\0', MAX_len);
  memset(pro_ss, '\0', MAX_len);
}

void
ssvm::
set_parameter(double H_bonus_, double f_clip_, double eta_, double rate_)
{
  H_bonus = H_bonus_; // Hamming_bonus
  f_clip = f_clip_;  // fobos_clip
  eta = eta_;    // weight bonus : default 0.001
  rate = rate_;  // weight threshold rate

  W = 1;
  all_edge = 0;
  all_protein = 0;
  all_RNA = 0;
  learn_num = 0;

  f_Rss5 = 1;/////
  f_Pss5 = 1;
  f_PsRs = 1;
  f_Ps3Rs3 = 1;
  f_AAp = 1;
  f_AAab = 1;
  f_AAh = 1;
  f_Rpp = 1;
  f_P1 = 1;
  f_P1R1 = 1;
  f_P3R1 = 1;
  f_P3R3 = 1;
  f_P5R1 = 1;
  f_P1R5 = 1;
  f_P5R5 = 1;

  out_train_weight = 0;
  out_cost = 0;
}

void
ssvm::
set_accuracy()
{
  predict_edge = 0;
  all_correct_edge = 0;
  correct_predict = 0;
  correct_predict_protein = 0;
  correct_predict_rna = 0;
}
 
void
ssvm::
init_weight()
{
  double WW, UW;
  WW = W - learn_num * f_clip;
  UW = W - learn_num * f_clip;

  // RNA ss 5
  if(f_Rss5>0){
    for(int i = 0; i < RNA_len; i++) {  
      sprintf(cRss5, "%c%c%c%c%c", rna_ss[i], rna_ss[i+1], rna_ss[i+2], rna_ss[i+3], rna_ss[i+4]);
      mRss5_w.insert(pair<string, double>(cRss5, WW));
      u_mRss5_w.insert(pair<string, double>(cRss5, UW));//unlabel
      sprintf(cRss5_r, "%c%c%c%c%c", rna_ss[i+4], rna_ss[i+3], rna_ss[i+2], rna_ss[i+1], rna_ss[i]);
      mRss5_w.insert(pair<string, double>(cRss5_r, WW));
      u_mRss5_w.insert(pair<string, double>(cRss5_r, UW));//unlabel
    }
  }
  // Protein ss
  if(f_Pss5>0){
    for(int i = 0; i < Pro_len; i++) {  
      sprintf(cPss5, "%c%c%c%c%c", pro_ss[i], pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], pro_ss[i+4]);
      mPss5_w.insert(pair<string, double>(cPss5, WW));
      u_mPss5_w.insert(pair<string, double>(cPss5, UW));//unlabel
      sprintf(cPss5_r, "%c%c%c%c%c", pro_ss[i+4], pro_ss[i+3], pro_ss[i+2], pro_ss[i+1], pro_ss[i]);
      mPss5_w.insert(pair<string, double>(cPss5_r, WW));
      u_mPss5_w.insert(pair<string, double>(cPss5_r, UW));//unlabel
    }
  }
  // Pro-RNA ss
  if(f_PsRs>0){
    for(int i = 0; i < Pro_len; i++) {  
      for(int j = 0; j < RNA_len; j++) {
	sprintf(cPsRs, "%c%c", pro_ss[i+2], rna_ss[j+2]);
	mPsRs_w.insert(pair<string, double>(cPsRs, WW));
	u_mPsRs_w.insert(pair<string, double>(cPsRs, UW));//unlabel
      }}
  }
  // Pro-RNA ss 3-3
  if(f_Ps3Rs3>0){  
    for(int i = 0; i < Pro_len; i++) {  
      for(int j = 0; j < RNA_len; j++) {
	sprintf(cPs3Rs3, "%c%c%c%c%c%c", pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], rna_ss[j+1], rna_ss[j+2], rna_ss[j+3]);
	mPs3Rs3_w.insert(pair<string, double>(cPs3Rs3, WW));
	u_mPs3Rs3_w.insert(pair<string, double>(cPs3Rs3, UW));//unlabel
	sprintf(cPs3Rs3_r, "%c%c%c%c%c%c", pro_ss[i+3], pro_ss[i+2], pro_ss[i+1], rna_ss[j+3], rna_ss[j+2], rna_ss[j+1]);
	mPs3Rs3_w.insert(pair<string, double>(cPs3Rs3_r, WW));
	u_mPs3Rs3_w.insert(pair<string, double>(cPs3Rs3_r, UW));//unlabel
      }}
  }
  // AA property
  if(f_AAp>0){ 
    for(int i = 0; i < 7; i++) {
      mAAp_w.insert(pair<int, double>(i, WW));
      u_mAAp_w.insert(pair<int, double>(i, UW));//unlabel
    }
  }
  // AA ab
  if(f_AAab>0){
    for(int i = 0; i < 3; i++) {
      mAAab_w.insert(pair<int, double>(i, WW));
      u_mAAab_w.insert(pair<int, double>(i, UW));//unlabel
    }
  }
  // AA h
  if(f_AAh>0){
    for(int i = 0; i < 3; i++) {  
      mAAh_w.insert(pair<int, double>(i, WW));
      u_mAAh_w.insert(pair<int, double>(i, UW));//unlabel
    }
  }
  // R pp
  if(f_Rpp>0){
    for(int i = 0; i < 2; i++) {  
      mRpp_w.insert(pair<int, double>(i, WW));
      u_mRpp_w.insert(pair<int, double>(i, UW));//unlabel
    }
  }
  // Protein 1 -
  if(f_P1>0){
    for(int i = 0; i < 20; i++) {
      mP1_w.insert(pair<char, double>(aa[i], WW));
      u_mP1_w.insert(pair<char, double>(aa[i], UW));//unlabel
    }
  }
  // Protein 1 - RNA 1
  if(f_P1R1>0){ 
    for(int i = 0; i < Pro_len; i++) {  
      for(int j = 0; j < 4; j++) {
	sprintf(cP1R1, "%c%c", protein[i+2], base[j]);
	mP1R1_w.insert(pair<string, double>(cP1R1, WW));
	u_mP1R1_w.insert(pair<string, double>(cP1R1, UW));//unlabel
      }}
  }
  // Protein 3 - RNA 1
  if(f_P3R1>0){ 
    for(int i = 0; i < Pro_len; i++) {  
      for(int j = 0; j < 4; j++) {
	sprintf(cP3R1, "%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], base[j]);
	mP3R1_w.insert(pair<string, double>(cP3R1, WW));
	u_mP3R1_w.insert(pair<string, double>(cP3R1, UW));//unlabel
	sprintf(cP3R1_r, "%c%c%c%c", protein[i+3], protein[i+2], protein[i+1], base[j]);
	mP3R1_w.insert(pair<string, double>(cP3R1_r, WW));
	u_mP3R1_w.insert(pair<string, double>(cP3R1_r, UW));//unlabel
      }}
  }
  // Protein 3 - RNA 3
  if(f_P3R3>0){
    for(int i = 0; i < Pro_len; i++) {  
      for(int j = 0; j < RNA_len; j++) {
	sprintf(cP3R3, "%c%c%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], RNA[j+1], RNA[j+2], RNA[j+3]);
	mP3R3_w.insert(pair<string, double>(cP3R3, WW));
	u_mP3R3_w.insert(pair<string, double>(cP3R3, UW));//unlabel
	sprintf(cP3R3_r, "%c%c%c%c%c%c", protein[i+3], protein[i+2], protein[i+1], RNA[j+3], RNA[j+2], RNA[j+1]);
	mP3R3_w.insert(pair<string, double>(cP3R3_r, WW));
	u_mP3R3_w.insert(pair<string, double>(cP3R3_r, UW));//unlabel
      }}
  }
  // Protein 5 - RNA 1
  if(f_P5R1>0){
    for(int i = 0; i < Pro_len; i++) {
      for(int j = 0; j < 4; j++) {
	sprintf(cP5R1, "%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], base[j]);
	mP5R1_w.insert(pair<string, double>(cP5R1, WW));
	u_mP5R1_w.insert(pair<string, double>(cP5R1, UW));//unlabel
	sprintf(cP5R1_r, "%c%c%c%c%c%c", protein[i+4], protein[i+3], protein[i+2], protein[i+1], protein[i], base[j]);
	mP5R1_w.insert(pair<string, double>(cP5R1_r, WW));
	u_mP5R1_w.insert(pair<string, double>(cP5R1_r, UW));//unlabel
      }}
  }
  // Protein 1 - RNA 5
  if(f_P1R5>0){
    for(int i = 0; i < Pro_len; i++) {
      for(int j = 0; j < RNA_len; j++) {
	sprintf(cP1R5, "%c%c%c%c%c%c", protein[i+2], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	mP1R5_w.insert(pair<string, double>(cP1R5, WW));
	u_mP1R5_w.insert(pair<string, double>(cP1R5, UW));//unlabel
	sprintf(cP1R5_r, "%c%c%c%c%c%c", protein[i+2], RNA[j+4], RNA[j+3], RNA[j+2], RNA[j+1], RNA[j]);
	mP1R5_w.insert(pair<string, double>(cP1R5_r, WW));
	u_mP1R5_w.insert(pair<string, double>(cP1R5_r, UW));//unlabel
      }} 
  }
  // Protein 5 - RNA 5
  if(f_P5R5>0){
    for(int i = 0; i < Pro_len; i++) {
      for(int j = 0; j < RNA_len; j++) {
	sprintf(cP5R5, "%c%c%c%c%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	mP5R5_w.insert(pair<string, double>(cP5R5, WW));
	u_mP5R5_w.insert(pair<string, double>(cP5R5, UW));//unlabel
	sprintf(cP5R5_r, "%c%c%c%c%c%c%c%c%c%c", protein[i+4], protein[i+3], protein[i+2], protein[i+1], protein[i], RNA[j+4], RNA[j+3], RNA[j+2], RNA[j+1], RNA[j]);
	mP5R5_w.insert(pair<string, double>(cP5R5_r, WW));
	u_mP5R5_w.insert(pair<string, double>(cP5R5_r, UW));//unlabel
      }} 
  }
}

void
ssvm::
init_graph_parameter()
{
  for(int i = 0; i < MAX_len*2; i++) {
    for(int j = 0; j < MAX_len*2; j++) {
      cost[i][j] = -INF;
      correct_match[i][j] = 0;   
    }}
}

void
ssvm::
init_match_parameter()
{
  for(int i = 0; i < MAX_len; i++) {
    correct_protein[i] = 0;
    correct_rna[i] = 0;
  }

  mRss5.clear(); mcorrect_Rss5.clear();
  mPss5.clear(); mcorrect_Pss5.clear();
  mPsRs.clear(); mcorrect_PsRs.clear();
  mPs3Rs3.clear(); mcorrect_Ps3Rs3.clear();

  for(int i = 0; i < 7; i++) {
    mAAp[i] = 0; mcorrect_AAp[i] = 0;
  }
  for(int i = 0; i < 3; i++) {
    mAAab[i] = 0; mcorrect_AAab[i] = 0;
  }
  for(int i = 0; i < 3; i++) {
    mAAh[i] = 0; mcorrect_AAh[i] = 0;
  }
  for(int i = 0; i < 2; i++) {
    mRpp[i] = 0; mcorrect_Rpp[i] = 0;
  }
  mP1.clear(); mcorrect_P1.clear();
  mP1R1.clear(); mcorrect_P1R1.clear();
  mP3R1.clear(); mcorrect_P3R1.clear();
  mP3R3.clear(); mcorrect_P3R3.clear();
  mP5R1.clear(); mcorrect_P5R1.clear();
  mP1R5.clear(); mcorrect_P1R5.clear();
  mP5R5.clear(); mcorrect_P5R5.clear();
}

void
ssvm::
update_feature_parameter(int un)
{
  if(result > 0) { // incorrect match > 0

    // RNA secondary structure : centroidfold
    if(f_Rss5>0){ 
   map<string, int>::iterator it_Rss5 = mRss5.begin();
    while( it_Rss5 != mRss5.end() ) {
      mRss5_w[(*it_Rss5).first] -= eta;
      ++it_Rss5;}
    map<string, int>::iterator it_coRss5 = mcorrect_Rss5.begin();
    while( it_coRss5 != mcorrect_Rss5.end() ) {
      mRss5_w[(*it_coRss5).first] += eta;
      ++it_coRss5;}
    }
    // Protein secondary structure : PDB or psipred
    if(f_Pss5>0){   
      map<string, int>::iterator it_Pss5 = mPss5.begin();
      while( it_Pss5 != mPss5.end() ) {
	mPss5_w[(*it_Pss5).first] -= eta;
	++it_Pss5;}
      map<string, int>::iterator it_coPss5 = mcorrect_Pss5.begin();
      while( it_coPss5 != mcorrect_Pss5.end() )	{
	mPss5_w[(*it_coPss5).first] += eta;
	++it_coPss5;}
    }
    // Pro-RNA ss
    if(f_PsRs>0){   
      map<string, int>::iterator it_PsRs = mPsRs.begin();
      while( it_PsRs != mPsRs.end() ) {
	mPsRs_w[(*it_PsRs).first] -= eta;
	++it_PsRs;}
      map<string, int>::iterator it_coPsRs = mcorrect_PsRs.begin();
      while( it_coPsRs != mcorrect_PsRs.end() )	{
	mPsRs_w[(*it_coPsRs).first] += eta;
	++it_coPsRs;}
    }
    // Pro-RNA ss 3-3
    if(f_Ps3Rs3>0){   
      map<string, int>::iterator it_Ps3Rs3 = mPs3Rs3.begin();
      while( it_Ps3Rs3 != mPs3Rs3.end() ) {
	mPs3Rs3_w[(*it_Ps3Rs3).first] -= eta;
	++it_Ps3Rs3;}
      map<string, int>::iterator it_coPs3Rs3 = mcorrect_Ps3Rs3.begin();
      while( it_coPs3Rs3 != mcorrect_Ps3Rs3.end() )	{
	mPs3Rs3_w[(*it_coPs3Rs3).first] += eta;
	++it_coPs3Rs3;}
    }
    // AA property
    if(f_AAp>0){  
      for(int i = 0; i < 7; i++) {
	if(mAAp[i] > 0) {
	  mAAp_w[i] -= eta;}
	if(mcorrect_AAp[i] > 0) {
	  mAAp_w[i] += eta;}
      }
    }
    // AA ab
    if(f_AAab>0){  
      for(int i = 0; i < 3; i++) {
	if(mAAab[i] > 0) {
	  mAAab_w[i] -= eta;}
	if(mcorrect_AAab[i] > 0) {
	  mAAab_w[i] += eta;}
      }
    }
    // AAh
    if(f_AAh>0){
      for(int i = 0; i < 3; i++) {
	if(mAAh[i] > 0) {
	  mAAh_w[i] -= eta;}
	if(mcorrect_AAh[i] > 0) {
	  mAAh_w[i] += eta;}
      }
    }
    // Rpp
    if(f_Rpp>0){ 
      for(int i = 0; i < 2; i++) {
	if(mRpp[i] > 0) {
	  mRpp_w[i] -= eta;}
	if(mcorrect_Rpp[i] > 0) {
	  mRpp_w[i] += eta;}
      }
    }
    // Protein 1 -
    if(f_P1>0){  
      for(int i = 0; i < 20; i++) {
	if(mP1[aa[i]] > 0) {
	  mP1_w[aa[i]] -= eta;}
	if(mcorrect_P1[aa[i]] > 0) {
	  mP1_w[aa[i]] += eta;}
      }
    }
    // Protein 1 - RNA 1
    if(f_P1R1>0){ 
      map<string, int>::iterator it_P1R1 = mP1R1.begin();
      while( it_P1R1 != mP1R1.end() ) {
	mP1R1_w[(*it_P1R1).first] -= eta;
	++it_P1R1;}
      map<string, int>::iterator it_coP1R1 = mcorrect_P1R1.begin();
      while( it_coP1R1 != mcorrect_P1R1.end() ) {
      mP1R1_w[(*it_coP1R1).first] += eta;
      ++it_coP1R1;}
    }
    // Protein 3 - RNA 1
    if(f_P3R1>0){   
      map<string, int>::iterator it_P3R1 = mP3R1.begin();
      while( it_P3R1 != mP3R1.end() ) {
	mP3R1_w[(*it_P3R1).first] -= eta;
	++it_P3R1;}
      map<string, int>::iterator it_coP3R1 = mcorrect_P3R1.begin();
      while( it_coP3R1 != mcorrect_P3R1.end() )	{
	mP3R1_w[(*it_coP3R1).first] += eta;
	++it_coP3R1;}
    }
    // Protein 3 - RNA 3
    if(f_P3R3>0){ 
      map<string, int>::iterator it_P3R3 = mP3R3.begin();
      while( it_P3R3 != mP3R3.end() ) {
	mP3R3_w[(*it_P3R3).first] -= eta;
	++it_P3R3;}
      map<string, int>::iterator it_coP3R3 = mcorrect_P3R3.begin();
      while( it_coP3R3 != mcorrect_P3R3.end() )	{
	mP3R3_w[(*it_coP3R3).first] += eta;
	++it_coP3R3;}
    }
    // Protein 5 - RNA 1
    if(f_P5R1>0){    
      map<string, int>::iterator it_P5R1 = mP5R1.begin();
      while( it_P5R1 != mP5R1.end() ) {
	mP5R1_w[(*it_P5R1).first] -= eta;
	++it_P5R1;}
      map<string, int>::iterator it_coP5R1 = mcorrect_P5R1.begin();
      while( it_coP5R1 != mcorrect_P5R1.end() )	{
	mP5R1_w[(*it_coP5R1).first] += eta;
	++it_coP5R1;} 
    }
    // Protein 1 - RNA 5
  if(f_P1R5>0){    
    map<string, int>::iterator it_P1R5 = mP1R5.begin();
    while( it_P1R5 != mP1R5.end() ) {
      mP1R5_w[(*it_P1R5).first] -= eta;
      ++it_P1R5;}
    map<string, int>::iterator it_coP1R5 = mcorrect_P1R5.begin();
    while( it_coP1R5 != mcorrect_P1R5.end() )	{
      mP1R5_w[(*it_coP1R5).first] += eta;
      ++it_coP1R5;}
  }
    // Protein 5 - RNA 5
  if(f_P5R5>0){    
    map<string, int>::iterator it_P5R5 = mP5R5.begin();
    while( it_P5R5 != mP5R5.end() ) {
      mP5R5_w[(*it_P5R5).first] -= eta;
      ++it_P5R5;}
    map<string, int>::iterator it_coP5R5 = mcorrect_P5R5.begin();
    while( it_coP5R5 != mcorrect_P5R5.end() )	{
      mP5R5_w[(*it_coP5R5).first] += eta;
      ++it_coP5R5;}
  }
  if(un == 1){unlabel_update();}
  }
  result = 0;
} 

//unlabel
void
ssvm::
unlabel_update()
{
  // RNA secondary structure : centroidfold
  if(f_Rss5>0){ 
  map<string, int>::iterator it_Rss5 = mRss5.begin();
  while( it_Rss5 != mRss5.end() ) {
    u_mRss5_w[(*it_Rss5).first] -= eta;
    ++it_Rss5;}
  map<string, int>::iterator it_coRss5 = mcorrect_Rss5.begin();
  while( it_coRss5 != mcorrect_Rss5.end() ) {
    u_mRss5_w[(*it_coRss5).first] += eta;
    ++it_coRss5;}
  }
  // Protein secondary structure : PDB or psipred
  if(f_Pss5>0){ 
  map<string, int>::iterator it_Pss5 = mPss5.begin();
  while( it_Pss5 != mPss5.end() ) {
    u_mPss5_w[(*it_Pss5).first] -= eta;
    ++it_Pss5;}
  map<string, int>::iterator it_coPss5 = mcorrect_Pss5.begin();
  while( it_coPss5 != mcorrect_Pss5.end() ) {
    u_mPss5_w[(*it_coPss5).first] += eta;
    ++it_coPss5;}
  }
  // Pro-RNA ss
  if(f_PsRs>0){ 
    map<string, int>::iterator it_PsRs = mPsRs.begin();
    while( it_PsRs != mPsRs.end() ) {
      u_mPsRs_w[(*it_PsRs).first] -= eta;
      ++it_PsRs;}
    map<string, int>::iterator it_coPsRs = mcorrect_PsRs.begin();
    while( it_coPsRs != mcorrect_PsRs.end() ) {
      u_mPsRs_w[(*it_coPsRs).first] += eta;
      ++it_coPsRs;}
  }
  // Pro-RNA ss 3-3
  if(f_Ps3Rs3>0){ 
    map<string, int>::iterator it_Ps3Rs3 = mPs3Rs3.begin();
    while( it_Ps3Rs3 != mPs3Rs3.end() ) {
      u_mPs3Rs3_w[(*it_Ps3Rs3).first] -= eta;
      ++it_Ps3Rs3;}
    map<string, int>::iterator it_coPs3Rs3 = mcorrect_Ps3Rs3.begin();
    while( it_coPs3Rs3 != mcorrect_Ps3Rs3.end() ) {
      u_mPs3Rs3_w[(*it_coPs3Rs3).first] += eta;
      ++it_coPs3Rs3;}
  }
  // AA property
  if(f_AAp>0){  
    for(int i = 0; i < 7; i++) {
      if(mAAp[i] > 0) {
	u_mAAp_w[i] -= eta;}
      if(mcorrect_AAp[i] > 0) {
	u_mAAp_w[i] += eta;}
    }
  }
  // AA ab
  if(f_AAab>0){  
    for(int i = 0; i < 3; i++) {
      if(mAAab[i] > 0) {
	u_mAAab_w[i] -= eta;}
      if(mcorrect_AAab[i] > 0) {
	u_mAAab_w[i] += eta;}
    }
  }
  // AAh
  if(f_AAh>0){ 
    for(int i = 0; i < 3; i++) {
      if(mAAh[i] > 0) {
	u_mAAh_w[i] -= eta;}
      if(mcorrect_AAh[i] > 0) {
	u_mAAh_w[i] += eta;}
    }
  }
  // Rpp
  if(f_Rpp>0){  
    for(int i = 0; i < 2; i++) {
      if(mRpp[i] > 0) {
	u_mRpp_w[i] -= eta;}
      if(mcorrect_Rpp[i] > 0) {
	u_mRpp_w[i] += eta;}
    }
  }
  // Protein 1 -
  if(f_P1>0){  
    for(int i = 0; i < 20; i++) {
      if(mP1[aa[i]] > 0) {
	u_mP1_w[aa[i]] -= eta;}
      if(mcorrect_P1[aa[i]] > 0) {
	u_mP1_w[aa[i]] += eta;}
    }
  }
  // Protein 1 - RNA 1
  if(f_P1R1>0){ 
    map<string, int>::iterator it_P1R1 = mP1R1.begin();
    while( it_P1R1 != mP1R1.end() ) {
      u_mP1R1_w[(*it_P1R1).first] -= eta;
      ++it_P1R1;}
    map<string, int>::iterator it_coP1R1 = mcorrect_P1R1.begin();
    while( it_coP1R1 != mcorrect_P1R1.end() ) {
      u_mP1R1_w[(*it_coP1R1).first] += eta;
      ++it_coP1R1;}
  }
  // Protein 3 - RNA 1
  if(f_P3R1>0){  
    map<string, int>::iterator it_P3R1 = mP3R1.begin();
    while( it_P3R1 != mP3R1.end() ) {
      u_mP3R1_w[(*it_P3R1).first] -= eta;
      ++it_P3R1;}
    map<string, int>::iterator it_coP3R1 = mcorrect_P3R1.begin();
    while( it_coP3R1 != mcorrect_P3R1.end() )	{
      u_mP3R1_w[(*it_coP3R1).first] += eta;
      ++it_coP3R1;}
  }
  // Protein 3 - RNA 3
  if(f_P3R3>0){ 
    map<string, int>::iterator it_P3R3 = mP3R3.begin();
    while( it_P3R3 != mP3R3.end() ) {
      u_mP3R3_w[(*it_P3R3).first] -= eta;
      ++it_P3R3;}
    map<string, int>::iterator it_coP3R3 = mcorrect_P3R3.begin();
    while( it_coP3R3 != mcorrect_P3R3.end() )	{
      u_mP3R3_w[(*it_coP3R3).first] += eta;
    ++it_coP3R3;}
  }
  // Protein 5 - RNA 1
  if(f_P5R1>0){ 
    map<string, int>::iterator it_P5R1 = mP5R1.begin();
    while( it_P5R1 != mP5R1.end() ) {
      u_mP5R1_w[(*it_P5R1).first] -= eta;
      ++it_P5R1;}
    map<string, int>::iterator it_coP5R1 = mcorrect_P5R1.begin();
    while( it_coP5R1 != mcorrect_P5R1.end() )	{
      u_mP5R1_w[(*it_coP5R1).first] += eta;
      ++it_coP5R1;}
  }
  // Protein 1 - RNA 5
  if(f_P1R5>0){ 
    map<string, int>::iterator it_P1R5 = mP1R5.begin();
    while( it_P1R5 != mP1R5.end() ) {
      u_mP1R5_w[(*it_P1R5).first] -= eta;
      ++it_P1R5;}
    map<string, int>::iterator it_coP1R5 = mcorrect_P1R5.begin();
    while( it_coP1R5 != mcorrect_P1R5.end() )	{
      u_mP1R5_w[(*it_coP1R5).first] += eta;
      ++it_coP1R5;} 
  }
  // Protein 5 - RNA 5
  if(f_P5R5>0){ 
    map<string, int>::iterator it_P5R5 = mP5R5.begin();
    while( it_P5R5 != mP5R5.end() ) {
      u_mP5R5_w[(*it_P5R5).first] -= eta;
      ++it_P5R5;}
    map<string, int>::iterator it_coP5R5 = mcorrect_P5R5.begin();
    while( it_coP5R5 != mcorrect_P5R5.end() )	{
      u_mP5R5_w[(*it_coP5R5).first] += eta;
      ++it_coP5R5;}
  } 
}

void
ssvm::
fobos()
{
  // RNA ss
  if(f_Rss5>0){ 
    map<string, double>::iterator it_Rss5w = mRss5_w.begin();
    while( it_Rss5w != mRss5_w.end() ) {
      mRss5_w[(*it_Rss5w).first] = clip(mRss5_w[(*it_Rss5w).first]);
      ++it_Rss5w;}
  }
  // Protein ss
  if(f_Pss5>0){ 
    map<string, double>::iterator it_Pss5w = mPss5_w.begin();
    while( it_Pss5w != mPss5_w.end() ) {
      mPss5_w[(*it_Pss5w).first] = clip(mPss5_w[(*it_Pss5w).first]);
      ++it_Pss5w;}
  }   
  // Pro-RNA ss
  if(f_PsRs>0){ 
    map<string, double>::iterator it_PsRsw = mPsRs_w.begin();
    while( it_PsRsw != mPsRs_w.end() ) {
      mPsRs_w[(*it_PsRsw).first] = clip(mPsRs_w[(*it_PsRsw).first]);
      ++it_PsRsw;}
  }   
  // Pro-RNA ss 3-3
  if(f_Ps3Rs3>0){ 
    map<string, double>::iterator it_Ps3Rs3w = mPs3Rs3_w.begin();
    while( it_Ps3Rs3w != mPs3Rs3_w.end() ) {
      mPs3Rs3_w[(*it_Ps3Rs3w).first] = clip(mPs3Rs3_w[(*it_Ps3Rs3w).first]);
      ++it_Ps3Rs3w;}   
  }
  // AA property
  if(f_AAp>0){ 
    for(int i = 0; i < 7; i++) {
      mAAp_w[i] = clip(mAAp_w[i]);}
  }
  // AA ab
  if(f_AAab>0){ 
    for(int i = 0; i < 3; i++) {
      mAAab_w[i] = clip(mAAab_w[i]);}
  }
  // AA h
  if(f_AAh>0){ 
    for(int i = 0; i < 3; i++) {
      mAAh_w[i] = clip(mAAh_w[i]);}
  }
  // Rpp
  if(f_Rpp>0){ 
    for(int i = 0; i < 2; i++) {
      mRpp_w[i] = clip(mRpp_w[i]);}
  }
  // Protein 1 -
  if(f_P1>0){ 
    for(int i = 0; i < 20; i++) {
      mP1_w[aa[i]] = clip(mP1_w[aa[i]]);}
  }
  // Protein 1 - RNA 1
  if(f_P1R1>0){ 
    map<string, double>::iterator it_P1R1w = mP1R1_w.begin();
    while( it_P1R1w != mP1R1_w.end() ) {
      mP1R1_w[(*it_P1R1w).first] = clip(mP1R1_w[(*it_P1R1w).first]);
      ++it_P1R1w;}   
  }
  // Protein 3 - RNA 1
  if(f_P3R1>0){ 
    map<string, double>::iterator it_P3R1w = mP3R1_w.begin();
    while( it_P3R1w != mP3R1_w.end() ) {
      mP3R1_w[(*it_P3R1w).first] = clip(mP3R1_w[(*it_P3R1w).first]);
      ++it_P3R1w;}   
  }
  // Protein 3 - RNA 3
  if(f_P3R3>0){
    map<string, double>::iterator it_P3R3w = mP3R3_w.begin();
    while( it_P3R3w != mP3R3_w.end() ) {
      mP3R3_w[(*it_P3R3w).first] = clip(mP3R3_w[(*it_P3R3w).first]);
      ++it_P3R3w;}   
  }
  // Protein 5 - RNA 1
  if(f_P5R1>0){ 
    map<string, double>::iterator it_P5R1w = mP5R1_w.begin();
    while( it_P5R1w != mP5R1_w.end() ) {
      mP5R1_w[(*it_P5R1w).first] = clip(mP5R1_w[(*it_P5R1w).first]);
      ++it_P5R1w;}   
  }
  // Protein 1 - RNA 5
  if(f_P1R5>0){ 
    map<string, double>::iterator it_P1R5w = mP1R5_w.begin();
    while( it_P1R5w != mP1R5_w.end() ) {
      mP1R5_w[(*it_P1R5w).first] = clip(mP1R5_w[(*it_P1R5w).first]);
      ++it_P1R5w;}   
  }
  // Protein 5 - RNA 5
  if(f_P5R5>0){
    map<string, double>::iterator it_P5R5w = mP5R5_w.begin();
    while( it_P5R5w != mP5R5_w.end() ) {
      mP5R5_w[(*it_P5R5w).first] = clip(mP5R5_w[(*it_P5R5w).first]);
      ++it_P5R5w; }   
  }
}
//unlabel
void
ssvm::
unlabel_fobos()
{
  // RNA ss
  if(f_Rss5>0){ 
    map<string, double>::iterator it_u_Rss5w = u_mRss5_w.begin();
    while( it_u_Rss5w != u_mRss5_w.end() ) {
      u_mRss5_w[(*it_u_Rss5w).first] = clip(u_mRss5_w[(*it_u_Rss5w).first]);
      ++it_u_Rss5w;}
  }
  // Protein ss
  if(f_Pss5>0){ 
    map<string, double>::iterator it_u_Pss5w = u_mPss5_w.begin();
    while( it_u_Pss5w != u_mPss5_w.end() ) {
      u_mPss5_w[(*it_u_Pss5w).first] = clip(u_mPss5_w[(*it_u_Pss5w).first]);
      ++it_u_Pss5w;}   
  }
  // Pro-RNA ss
  if(f_PsRs>0){ 
    map<string, double>::iterator it_u_PsRsw = u_mPsRs_w.begin();
    while( it_u_PsRsw != u_mPsRs_w.end() ) {
      u_mPsRs_w[(*it_u_PsRsw).first] = clip(u_mPsRs_w[(*it_u_PsRsw).first]);
      ++it_u_PsRsw;}   
  }
  // Pro-RNA ss 3-3
  if(f_Ps3Rs3>0){ 
    map<string, double>::iterator it_u_Ps3Rs3w = u_mPs3Rs3_w.begin();
    while( it_u_Ps3Rs3w != u_mPs3Rs3_w.end() ) {
      u_mPs3Rs3_w[(*it_u_Ps3Rs3w).first] = clip(u_mPs3Rs3_w[(*it_u_Ps3Rs3w).first]);
      ++it_u_Ps3Rs3w;}   
  }
  // AA property
  if(f_AAp>0){ 
    for(int i = 0; i < 7; i++) {
      u_mAAp_w[i] = clip(u_mAAp_w[i]);}
  }
  // AA ab
  if(f_AAab>0){ 
    for(int i = 0; i < 3; i++) {
      u_mAAab_w[i] = clip(u_mAAab_w[i]);}
  }
  // AA h
  if(f_AAh>0){ 
    for(int i = 0; i < 3; i++) {
      u_mAAh_w[i] = clip(u_mAAh_w[i]);}
  }
  // Rpp
  if(f_Rpp>0){ 
    for(int i = 0; i < 2; i++) {
      u_mRpp_w[i] = clip(u_mRpp_w[i]);}
  }
  // Protein 1 -
  if(f_P1>0){ 
    for(int i = 0; i < 20; i++) {
      u_mP1_w[aa[i]] = clip(u_mP1_w[aa[i]]);}
  }
  // Protein 1 - RNA 1
  if(f_P1R1>0){ 
    map<string, double>::iterator it_u_P1R1w = u_mP1R1_w.begin();
    while( it_u_P1R1w != u_mP1R1_w.end() ) {
      u_mP1R1_w[(*it_u_P1R1w).first] = clip(u_mP1R1_w[(*it_u_P1R1w).first]);
      ++it_u_P1R1w;}   
  }
  // Protein 3 - RNA 1
  if(f_P3R1>0){ 
    map<string, double>::iterator it_u_P3R1w = u_mP3R1_w.begin();
    while( it_u_P3R1w != u_mP3R1_w.end() ) {
      u_mP3R1_w[(*it_u_P3R1w).first] = clip(u_mP3R1_w[(*it_u_P3R1w).first]);
      ++it_u_P3R1w;}   
  }
  // Protein 3 - RNA 3
  if(f_P3R3>0){
    map<string, double>::iterator it_u_P3R3w = u_mP3R3_w.begin();
    while( it_u_P3R3w != u_mP3R3_w.end() ) {
      u_mP3R3_w[(*it_u_P3R3w).first] = clip(u_mP3R3_w[(*it_u_P3R3w).first]);
      ++it_u_P3R3w;}   
  }
  // Protein 5 - RNA 1
  if(f_P5R1>0){
    map<string, double>::iterator it_u_P5R1w = u_mP5R1_w.begin();
    while( it_u_P5R1w != u_mP5R1_w.end() ) {
      u_mP5R1_w[(*it_u_P5R1w).first] = clip(u_mP5R1_w[(*it_u_P5R1w).first]);
      ++it_u_P5R1w;}   
  }
  // Protein 1 - RNA 5
  if(f_P1R5>0){
    map<string, double>::iterator it_u_P1R5w = u_mP1R5_w.begin();
    while( it_u_P1R5w != u_mP1R5_w.end() ) {
      u_mP1R5_w[(*it_u_P1R5w).first] = clip(u_mP1R5_w[(*it_u_P1R5w).first]);
      ++it_u_P1R5w;}   
  }
  // Protein 5 - RNA 5
  if(f_P5R5>0){ 
    map<string, double>::iterator it_u_P5R5w = u_mP5R5_w.begin();
    while( it_u_P5R5w != u_mP5R5_w.end() ) {
      u_mP5R5_w[(*it_u_P5R5w).first] = clip(u_mP5R5_w[(*it_u_P5R5w).first]);
      ++it_u_P5R5w; }   
  } 
}

double
ssvm::
clip(double w)
{
  if(w >= 0){
    if(w > f_clip){
      return w - f_clip;}
    else {return 0;}}
  else {return -1 * clip(-1 * w);}
}

void
ssvm::
update_cost()
{
  for(int i = 0; i < Pro_len; i++) {
    for(int j = 0; j < RNA_len; j++) {
      double w_all = 0;
      // RNAss 
      double wRss5 = 0;     
      if(f_Rss5>0){
	w_all += W;
	sprintf(cRss5, "%c%c%c%c%c", rna_ss[j], rna_ss[j+1], rna_ss[j+2], rna_ss[j+3], rna_ss[j+4]);
	if(mRss5_w[cRss5] > 0 ) { wRss5 = mRss5_w[cRss5];}
	if(u_rate_mRss5[cRss5] > 0) { wRss5 += u_rate_mRss5[cRss5] * u_mRss5_w[cRss5];} // unlabeled
      }
      // Proteinss 
      double wPss5 = 0;     
      if(f_Pss5>0){     
	w_all += W;
	sprintf(cPss5, "%c%c%c%c%c", pro_ss[i], pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], pro_ss[i+4]);
	if(mPss5_w[cPss5] > 0) { wPss5 = mPss5_w[cPss5];}
	if(u_rate_mPss5[cPss5] > 0) { wPss5 += u_rate_mPss5[cPss5] * u_mPss5_w[cPss5];} // unlabeled
      }
      // Pro-RNAss 
      double wPsRs = 0;     
      if(f_PsRs>0){     
	w_all += W;
	sprintf(cPsRs, "%c%c", pro_ss[i+2], rna_ss[j+2]);
	if(mPsRs_w[cPsRs] > 0) { wPsRs = mPsRs_w[cPsRs];}
	if(u_rate_mPsRs[cPsRs] > 0) { wPsRs += u_rate_mPsRs[cPsRs] * u_mPsRs_w[cPsRs];} // unlabeled
      }
      // Pro-RNA ss 3-3
      double wPs3Rs3 = 0;     
      if(f_Ps3Rs3>0){    
	w_all += W;
	sprintf(cPs3Rs3, "%c%c%c%c%c%c", pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], rna_ss[j+1], rna_ss[j+2], rna_ss[j+3]);
	if(mPs3Rs3_w[cPs3Rs3] > 0) { wPs3Rs3 = mPs3Rs3_w[cPs3Rs3];}
	if(u_rate_mPs3Rs3[cPs3Rs3] > 0) { wPs3Rs3 += u_rate_mPs3Rs3[cPs3Rs3] * u_mPs3Rs3_w[cPs3Rs3];} // unlabeled
      }
      // AAp
      double wAAp = 0;
      if(f_AAp>0){     
	w_all += W;
	int Apn = AA_7group(protein[i + 2]);
	if(mAAp_w[Apn] > 0) { wAAp = mAAp_w[Apn];}
	if(u_rate_mAAp[Apn] > 0) { wAAp += u_rate_mAAp[Apn] * u_mAAp_w[Apn];} // unlabeled
      }
      // AAab
      double wAAab = 0;
      if(f_AAab>0){     
	w_all += W;
	int Aabn = AA_3ab(protein[i + 2]);
	if(mAAab_w[Aabn] > 0) { wAAab = mAAab_w[Aabn];}
	if(u_rate_mAAab[Aabn] > 0) { wAAab += u_rate_mAAab[Aabn] * u_mAAab_w[Aabn];} // unlabeled
      }
      // AAh
      double wAAh = 0;
      if(f_AAh>0){
	w_all += W;
	int Ahn = AA_3hydro(protein[i + 2]);
	if(mAAh_w[Ahn] > 0) { wAAh = mAAh_w[Ahn];}
	if(u_rate_mAAh[Ahn] > 0) { wAAh += u_rate_mAAh[Ahn] * u_mAAh_w[Ahn];} // unlabeled
      }
      // Rpp
      double wRpp = 0;
      if(f_Rpp>0){     
	w_all += W;
	int Rpn = RNA_2pp(RNA[j + 2]);
	if(mRpp_w[Rpn] > 0) { wRpp = mRpp_w[Rpn];}
	if(u_rate_mRpp[Rpn] > 0) { wRpp += u_rate_mRpp[Rpn] * u_mRpp_w[Rpn];} // unlabeled
      }
      // Protein 1 -
      double wP1 = 0;
      if(f_P1>0){     
	w_all += W;
	if(mP1_w[protein[i+2]] > 0) { wP1 = mP1_w[protein[i+2]];}
	if(u_rate_mP1[protein[i+2]] > 0) { wP1 += u_rate_mP1[protein[i+2]] * u_mP1_w[protein[i+2]];} // unlabeled
      }
      // Protein 1 - RNA 1
      double wP1R1 = 0;
      if(f_P1R1>0){ 
	w_all += W;
	sprintf(cP1R1, "%c%c", protein[i+2], RNA[j+2]);
	if(mP1R1_w[cP1R1] > 0) { wP1R1 = mP1R1_w[cP1R1];}
	if(u_rate_mP1R1[cP1R1] > 0) { wP1R1 += u_rate_mP1R1[cP1R1] * u_mP1R1_w[cP1R1];} // unlabeled
      }
      // Protein 3 - RNA 1
      double wP3R1 = 0;    
      if(f_P3R1>0){    
	w_all += W;
	sprintf(cP3R1, "%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], RNA[j+2]);
	if(mP3R1_w[cP3R1] > 0) { wP3R1 = mP3R1_w[cP3R1];}
	if(u_rate_mP3R1[cP3R1] > 0) { wP3R1 += u_rate_mP3R1[cP3R1] * u_mP3R1_w[cP3R1];} // unlabeled
      }
      // Protein 3 - RNA 3
      double wP3R3 = 0;    
      if(f_P3R3>0){    
	w_all += W;
	sprintf(cP3R3, "%c%c%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], RNA[j+1], RNA[j+2], RNA[j+3]);
	if(mP3R3_w[cP3R3] > 0) { wP3R3 = mP3R3_w[cP3R3];}
	if(u_rate_mP3R3[cP3R3] > 0) { wP3R3 += u_rate_mP3R3[cP3R3] * u_mP3R3_w[cP3R3];} // unlabeled
      }
      // Protein 5 - RNA 1
      double wP5R1 = 0;     
      if(f_P5R1>0){     
	w_all += W;
	sprintf(cP5R1, "%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], RNA[j+2]);
	if(mP5R1_w[cP5R1] > 0) { wP5R1 = mP5R1_w[cP5R1];}
	if(u_rate_mP5R1[cP5R1] > 0) { wP5R1 += u_rate_mP5R1[cP5R1] * u_mP5R1_w[cP5R1];} // unlabeled
      }
      // Protein 1 - RNA 5
      double wP1R5 = 0;    
      if(f_P1R5>0){     
	w_all += W;
	sprintf(cP1R5, "%c%c%c%c%c%c", protein[i+2], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	if(mP1R5_w[cP1R5] > 0) { wP1R5 = mP1R5_w[cP1R5];}
	if(u_rate_mP1R5[cP1R5] > 0) { wP1R5 += u_rate_mP1R5[cP1R5] * u_mP1R5_w[cP1R5];} // unlabeled
      }
      // Protein 5 - RNA 5
      double wP5R5 = 0;   
      if(f_P5R5>0){     
	w_all += W;
	sprintf(cP5R5, "%c%c%c%c%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	if(mP5R5_w[cP5R5] > 0) { wP5R5 = mP5R5_w[cP5R5];}
	if(u_rate_mP5R5[cP5R5] > 0) { wP5R5 += u_rate_mP5R5[cP5R5] * u_mP5R5_w[cP5R5];} //unlabeled
      }
      cost[i + 1][j + Pro_len + 1] = 
	w_all - (wRss5 + wPss5 + wPsRs + wPs3Rs3 + wAAp + wAAab + wAAh + wRpp + wP1 + wP1R1 + wP3R1 + wP3R3 + wP5R1 + wP1R5 + wP5R5);
      if(unlabel == 2){
	if(cost[i + 1][j + Pro_len + 1] >= (w_all - (w_all/W) * learn_num * f_clip) * (1 - rate)) {
	  cost[i + 1][j + Pro_len + 1] = -INF;
	}
      }
      if(-INF < cost[i + 1][j + Pro_len + 1] && cost[i + 1][j + Pro_len + 1] <= 0) {
	cost[i + 1][j + Pro_len + 1] = 1;
      }
    }}
}

// unlabeled
void
ssvm::
unlabeled_count()
{
  // RNA secondary structure : centroidfold
  if(f_Rss5>0){ 
    map<string, int>::iterator it_u_Rss5 = u_mRss5.begin();
    while( it_u_Rss5 != u_mRss5.end() ) {
      if((*it_u_Rss5).second > 0) {
	u_rate_mRss5[(*it_u_Rss5).first] =  u_mRss5[(*it_u_Rss5).first] / un_predict_edge;}
      ++it_u_Rss5;
    }
  }
  // Protein secondary structure : PDB or psipred
  if(f_Pss5>0){ 
    map<string, int>::iterator it_u_Pss5 = u_mPss5.begin();
    while( it_u_Pss5 != u_mPss5.end() ) {
      if((*it_u_Pss5).second > 0) {
	u_rate_mPss5[(*it_u_Pss5).first] = u_mPss5[(*it_u_Pss5).first] / un_predict_edge;}
      ++it_u_Pss5;
    }
  }
  // Pro-RNA ss
  if(f_PsRs>0){ 
    map<string, int>::iterator it_u_PsRs = u_mPsRs.begin();
    while( it_u_PsRs != u_mPsRs.end() ) {
      if((*it_u_PsRs).second > 0) {
	u_rate_mPsRs[(*it_u_PsRs).first] = u_mPsRs[(*it_u_PsRs).first] / un_predict_edge;}
      ++it_u_PsRs;
    }
  }
  // Pro-RNA ss 3-3
  if(f_Ps3Rs3>0){ 
    map<string, int>::iterator it_u_Ps3Rs3 = u_mPs3Rs3.begin();
    while( it_u_Ps3Rs3 != u_mPs3Rs3.end() ) {
      if((*it_u_Ps3Rs3).second > 0) {
	u_rate_mPs3Rs3[(*it_u_Ps3Rs3).first] = u_mPs3Rs3[(*it_u_Ps3Rs3).first] / un_predict_edge;}
      ++it_u_Ps3Rs3;
    }
  }
  // AA property
  if(f_AAp>0){ 
    for(int i = 0; i < 7; i++) {
      if(u_mAAp[i] > 0) {
	u_rate_mAAp[i] = u_mAAp[i] / un_predict_edge;}
    }
  }
  // AA ab
  if(f_AAab>0){ 
    for(int i = 0; i < 3; i++) {
      if(u_mAAab[i] > 0) {
	u_rate_mAAab[i] = u_mAAab[i] / un_predict_edge;}
    }
  }
  // AAh
  if(f_AAh>0){ 
    for(int i = 0; i < 3; i++) {
      if(u_mAAh[i] > 0) {
	u_rate_mAAh[i] = u_mAAh[i] / un_predict_edge;}
    }
  }
  // Rpp
  if(f_Rpp>0){ 
    for(int i = 0; i < 2; i++) {
      if(u_mRpp[i] > 0) {
	u_rate_mRpp[i] = u_mRpp[i] / un_predict_edge;}
    }
  }
  // Protein 1 -
  if(f_P1>0){ 
    for(int i = 0; i < 20; i++) {
      if(u_mP1[aa[i]] > 0) {
	u_rate_mP1[aa[i]] = u_mP1[aa[i]] / un_predict_edge;}
    }
  }
  // Protein 1 - RNA 1 // unlabeled
  if(f_P1R1>0){ 
    map<string, double>::iterator it_u_P1R1 = u_mP1R1.begin();
    while( it_u_P1R1 != u_mP1R1.end() ) {
      if((*it_u_P1R1).second > 0) {
	u_rate_mP1R1[(*it_u_P1R1).first] = u_mP1R1[(*it_u_P1R1).first] / un_predict_edge;}
      ++it_u_P1R1;
    }
  }
  // Protein 3 - RNA 1 // unlabeled
  if(f_P3R1>0){ 
    map<string, double>::iterator it_u_P3R1 = u_mP3R1.begin();
    while( it_u_P3R1 != u_mP3R1.end() ) {
      if((*it_u_P3R1).second > 0) {
	u_rate_mP3R1[(*it_u_P3R1).first] = u_mP3R1[(*it_u_P3R1).first] / un_predict_edge;}
      ++it_u_P3R1;
    }
  }
  // Protein 3 - RNA 3 // unlabeled
  if(f_P3R3>0){ 
    map<string, double>::iterator it_u_P3R3 = u_mP3R3.begin();
    while( it_u_P3R3 != u_mP3R3.end() ) {
      if((*it_u_P3R3).second > 0) {
	u_rate_mP3R3[(*it_u_P3R3).first] = u_mP3R3[(*it_u_P3R3).first] / un_predict_edge;}
      ++it_u_P3R3;
    }
  }
  // Protein 5 - RNA 1 // unlabeled
  if(f_P5R1>0){ 
    map<string, double>::iterator it_u_P5R1 = u_mP5R1.begin();
    while( it_u_P5R1 != u_mP5R1.end() ) {
      if((*it_u_P5R1).second > 0) {
	u_rate_mP5R1[(*it_u_P5R1).first] = u_mP5R1[(*it_u_P5R1).first] / un_predict_edge;}
      ++it_u_P5R1;
    }
  }
  // Protein 1 - RNA 5 // unlabeled
  if(f_P1R5>0){ 
    map<string, double>::iterator it_u_P1R5 = u_mP1R5.begin();
    while( it_u_P1R5 != u_mP1R5.end() ) {
      if((*it_u_P1R5).second > 0) {
	u_rate_mP1R5[(*it_u_P1R5).first] = u_mP1R5[(*it_u_P1R5).first] / un_predict_edge;}
      ++it_u_P1R5;
    }
  }
  // Protein 5 - RNA 5 // unlabeled
  if(f_P5R5>0){ 
    map<string, double>::iterator it_u_P5R5 = u_mP5R5.begin();
    while( it_u_P5R5 != u_mP5R5.end() ) {
      if((*it_u_P5R5).second > 0) {
	u_rate_mP5R5[(*it_u_P5R5).first] = u_mP5R5[(*it_u_P5R5).first] / un_predict_edge;}
      ++it_u_P5R5;
    }
  }
}

void
ssvm::
correct_matching(string match)
{
  int F, T;
  FILE *fp = fopen(match.c_str(), "r");
  while(fscanf(fp, "%d %d", &F, &T) != EOF) {
    correct_match[F + 1][T + Pro_len + 1]++;
    all_correct_edge++;
    correct_protein[F]++;
    correct_rna[T]++;
    if(output>0){ printf("\n");}

    // RNAss
    if(f_Rss5>0){ 
      sprintf(cRss5, "%c%c%c%c%c", rna_ss[T], rna_ss[T+1], rna_ss[T+2], rna_ss[T+3], rna_ss[T+4]);
      if(output>0){ printf("%c%c%c%c%c (%d) | RNA secondary structure\n", rna_ss[T], rna_ss[T+1], rna_ss[T+2], rna_ss[T+3], rna_ss[T+4], T);}
      mcorrect_Rss5.insert(pair<string, double>(cRss5, 0));
      mcorrect_Rss5[cRss5]++;
      sprintf(cRss5_r, "%c%c%c%c%c", rna_ss[T+4], rna_ss[T+3], rna_ss[T+2], rna_ss[T+1], rna_ss[T]);
      if(memcmp(cRss5, cRss5_r, 5) != 0){
	mcorrect_Rss5.insert(pair<string, double>(cRss5_r, 0));
	mcorrect_Rss5[cRss5_r]++;}  
    }
    // Proteinss
    if(f_Pss5>0){ 
      sprintf(cPss5, "%c%c%c%c%c", pro_ss[F], pro_ss[F+1], pro_ss[F+2], pro_ss[F+3], pro_ss[F+4]);
      if(output>0){ printf("%c%c%c%c%c (%d) | Protein secondary structure\n", pro_ss[F], pro_ss[F+1], pro_ss[F+2], pro_ss[F+3], pro_ss[F+4], F);}
      mcorrect_Pss5.insert(pair<string, double>(cPss5, 0));
      mcorrect_Pss5[cPss5]++;
      sprintf(cPss5_r, "%c%c%c%c%c", pro_ss[F+4], pro_ss[F+3], pro_ss[F+2], pro_ss[F+1], pro_ss[F]);
      if(memcmp(cPss5, cPss5_r, 5) != 0){
	mcorrect_Pss5.insert(pair<string, double>(cPss5_r, 0));
	mcorrect_Pss5[cPss5_r]++;}
    }
    // Pro-RNA ss
    if(f_PsRs>0){ 
      sprintf(cPsRs, "%c%c", pro_ss[F+2], rna_ss[T+2]);
      if(output>0){ printf("%c%c  (%d - %d) | Pro-RNA secondary structure\n", pro_ss[F+2], rna_ss[T+2], F, T);}
      mcorrect_PsRs.insert(pair<string, double>(cPsRs, 0));
      mcorrect_PsRs[cPsRs]++;
    }
    // Pro-RNA ss 3-3
    if(f_Ps3Rs3>0){ 
      sprintf(cPs3Rs3, "%c%c%c%c%c%c", pro_ss[F+1], pro_ss[F+2], pro_ss[F+3], rna_ss[T+1], rna_ss[T+2], rna_ss[T+3]);
      if(output>0){  printf("%c%c%c %c%c%c (%d - %d)     | P-R ss 3-3\n", pro_ss[F+1], pro_ss[F+2], pro_ss[F+3], rna_ss[T+1], rna_ss[T+2], rna_ss[T+3], F, T);}
      mcorrect_Ps3Rs3.insert(pair<string, double>(cPs3Rs3, 0));     
      mcorrect_Ps3Rs3[cPs3Rs3]++;
      sprintf(cPs3Rs3_r, "%c%c%c%c%c%c", pro_ss[F+3], pro_ss[F+2], pro_ss[F+1], rna_ss[T+3], rna_ss[T+2], rna_ss[T+1]);
      if(memcmp(cPs3Rs3, cPs3Rs3_r, 6) != 0){
	mcorrect_Ps3Rs3.insert(pair<string, double>(cPs3Rs3_r, 0));     
	mcorrect_Ps3Rs3[cPs3Rs3_r]++;}
    }
    // AAp
    if(f_AAp>0){ 
      int Apn = AA_7group(protein[F+2]);
      mcorrect_AAp[Apn]++;
      if(output>0){ printf("%c:%d (%d) | amino 7 groups : AGV, ILFP, YMTS, HNQW, RK, DE, C\n", protein[F+2], Apn+1, F);}
    }
    // AAab
    if(f_AAab>0){    
      int Aabn = AA_3ab(protein[F+2]);
      mcorrect_AAab[Aabn]++;
      if(output>0){ printf("%c:%d (%d)    | amino 0:acid, 1:base, 2:neutral\n", protein[F+2], Aabn, F);}
    }
    // AAh
    if(f_AAh>0){    
      int Ahn = AA_3hydro(protein[F+2]);
      if(output>0){ printf("%c:%d (%d)       | amino 0:hydrophobic, 1:hydrophilic, 2:neutral\n", protein[F+2], Ahn, F);}
      mcorrect_AAh[Ahn]++;
    }
    // Rpp
    if(f_Rpp>0){ 
      int Rpn = RNA_2pp(RNA[T+2]);
      if(output>0){ printf("%c:%d (%d)       | nucleotide 0:purine, 1:pyrimidine\n", RNA[T+2], Rpn, T);}
      mcorrect_Rpp[Rpn]++;
    }
    // Protein 1 -
    if(f_P1>0){ 
      mcorrect_P1.insert(pair<char, double>(protein[F+2], 0));    
      mcorrect_P1[protein[F+2]]++;
      if(output>0){ printf("%c (%d)           | protein 1-mer\n", protein[F+2], F);}
    }
    // Protein 1 - RNA 1
    if(f_P1R1>0){ 
      sprintf(cP1R1, "%c%c", protein[F+2], RNA[T+2]);
      if(output>0){ printf("%c %c (%d - %d)         | P-R : 1-1\n", protein[F+2], RNA[T+2], F, T);}
      mcorrect_P1R1.insert(pair<string, double>(cP1R1, 0));
      mcorrect_P1R1[cP1R1]++;
    }
    // Protein 3 - RNA 1
    if(f_P3R1>0){ 
      sprintf(cP3R1, "%c%c%c%c", protein[F+1], protein[F+2], protein[F+3], RNA[T+2]);
      mcorrect_P3R1.insert(pair<string, double>(cP3R1, 0));     
      mcorrect_P3R1[cP3R1]++;
      if(output>0){ printf("%c%c%c %c (%d - %d)       | P-R : 3-1\n", protein[F+1], protein[F+2], protein[F+3], RNA[T+2], F, T);}
      sprintf(cP3R1_r, "%c%c%c%c", protein[F+3], protein[F+2], protein[F+1], RNA[T+2]);
      if(memcmp(cP3R1, cP3R1_r, 4) != 0){
	mcorrect_P3R1.insert(pair<string, double>(cP3R1_r, 0));     
	mcorrect_P3R1[cP3R1_r]++;} 
    }
    // Protein 3 - RNA 3
    if(f_P3R3>0){ 
      sprintf(cP3R3, "%c%c%c%c%c%c", protein[F+1], protein[F+2], protein[F+3], RNA[T+1], RNA[T+2], RNA[T+3]);
      if(output>0){ printf("%c%c%c %c%c%c (%d - %d)     | P-R : 3-3\n", protein[F+1], protein[F+2], protein[F+3], RNA[T+1], RNA[T+2], RNA[T+3], F, T);}
      mcorrect_P3R3.insert(pair<string, double>(cP3R3, 0));     
      mcorrect_P3R3[cP3R3]++;
      sprintf(cP3R3_r, "%c%c%c%c%c%c", protein[F+3], protein[F+2], protein[F+1], RNA[T+3], RNA[T+2], RNA[T+1]);
      if(memcmp(cP3R3, cP3R3_r, 6) != 0){
	mcorrect_P3R3.insert(pair<string, double>(cP3R3_r, 0));     
	mcorrect_P3R3[cP3R3_r]++;}
    }
    // Protein 5 - RNA 1
    if(f_P5R1>0){  
      sprintf(cP5R1, "%c%c%c%c%c%c", protein[F], protein[F+1], protein[F+2], protein[F+3], protein[F+4], RNA[T+2]);
      if(output>0){ printf("%c%c%c%c%c %c (%d - %d)     | P-R : 5-1\n", protein[F], protein[F+1], protein[F+2], protein[F+3], protein[F+4], RNA[T+2], F, T);}
      mcorrect_P5R1.insert(pair<string, double>(cP5R1, 0));     
      mcorrect_P5R1[cP5R1]++;
      sprintf(cP5R1_r, "%c%c%c%c%c%c", protein[F+4], protein[F+3], protein[F+2], protein[F+1], protein[F], RNA[T+2]);
      if(memcmp(cP5R1, cP5R1_r, 6) != 0){
	mcorrect_P5R1.insert(pair<string, double>(cP5R1_r, 0));     
	mcorrect_P5R1[cP5R1_r]++;}
    }
    // Protein 1 - RNA 5
    if(f_P1R5>0){   
      sprintf(cP1R5, "%c%c%c%c%c%c", protein[F+2], RNA[T], RNA[T+1], RNA[T+2], RNA[T+3], RNA[T+4]);
      if(output>0){ printf("%c %c%c%c%c%c (%d - %d)     | P-R : 1-5\n", protein[F+2], RNA[T], RNA[T+1], RNA[T+2], RNA[T+3], RNA[T+4], F, T);}
      mcorrect_P1R5.insert(pair<string, double>(cP1R5, 0));     
      mcorrect_P1R5[cP1R5]++;
      sprintf(cP1R5_r, "%c%c%c%c%c%c", protein[F+2], RNA[T+4], RNA[T+3], RNA[T+2], RNA[T+1], RNA[T]);
      if(memcmp(cP1R5, cP1R5_r, 6) != 0){
	mcorrect_P1R5.insert(pair<string, double>(cP1R5_r, 0));     
	mcorrect_P1R5[cP1R5_r]++;}
    }
    // Protein 5 - RNA 5
    if(f_P5R5>0){   
      sprintf(cP5R5, "%c%c%c%c%c%c%c%c%c%c", protein[F], protein[F+1], protein[F+2], protein[F+3], protein[F+4], RNA[T], RNA[T+1], RNA[T+2], RNA[T+3], RNA[T+4]);
      if(output>0){ printf("%c%c%c%c%c %c%c%c%c%c (%d - %d) | P-R : 5-5\n", protein[F], protein[F+1], protein[F+2], protein[F+3], protein[F+4], RNA[T], RNA[T+1], RNA[T+2], RNA[T+3], RNA[T+4], F, T);}
      mcorrect_P5R5.insert(pair<string, double>(cP5R5, 0));        
      mcorrect_P5R5[cP5R5]++;
      sprintf(cP5R5_r, "%c%c%c%c%c%c%c%c%c%c", protein[F+4], protein[F+3], protein[F+2], protein[F+1], protein[F], RNA[T+4], RNA[T+3], RNA[T+2], RNA[T+1], RNA[T]);
      if(memcmp(cP5R5, cP5R5_r, 10) != 0){
	mcorrect_P5R5.insert(pair<string, double>(cP5R5, 0));     
	mcorrect_P5R5[cP5R5_r]++;}
    }
  }
  fclose(fp);
  if(output>0){ cout << endl;}
}
 
void
ssvm::
Hamming_distance()
{
  for(int i = 1; i < Pro_len; i++) {
    for(int j = 1 + Pro_len; j < Pro_len + RNA_len + 1; j++) {
      if(correct_match[i][j] > 0 && cost[i][j] > (-1) * H_bonus) {
	cost[i][j] += H_bonus;
      }
      if(correct_match[i][j] == 0 && cost[i][j] > H_bonus)
	cost[i][j] -= H_bonus;
    }}
}

void
ssvm::
add_edge(int from, int to, double w)
{
  cost[from][to] = w;
}

int
ssvm::
BellmanFord(int s, int t)
{
  int V = t + 1;
  double dist[MAX_len*2];       // min_distance from s
  int prev[MAX_len*2];      // previous node
  int path[MAX_len*2];
  fill(dist, dist + MAX_len*2, INF);
  fill(prev, prev + MAX_len*2, -1);
  fill(path, path + MAX_len*2, -1);
  dist[s] = 0;

  int limit = 0;
  bool update = true;
  while(update) {     //get min_dis from s
    update = false;
    for(int v = 0; v < V - 1; v++) {
      if(dist[v] == INF) {continue;}
      for(int u = 0; u < V; u++) {
       	if(dist[v] + cost[v][u] > 0 && dist[u] > dist[v] + cost[v][u]) {
	  dist[u] = dist[v] + cost[v][u];
	  prev[u] = v;
	  update = true;
	}}}
    limit++;
    if(limit > MAX_len){
      if(output>0){ cout << "loop --- bad parameter" << endl; return -1;}
    }
  }
  if(limit == 0){return -1;}
  limit = 0;
  if(prev[t] == -1 || dist[t] == INF) {return -1;} // no more path
  int p = 0;
  while(t >= 0) {
    path[p] = t;
    p++;
    t = prev[t];
    if(p > MAX_len) {
      if(output>0){ cout << "loop --- bad parameter" << endl; return -1;}
    }
  }

  for(int i = 0; i < p - 1; i++) {   
    cost[path[i]][path[i + 1]] = (-1) * cost[path[i + 1]][path[i]];
    cost[path[i + 1]][path[i]] = -INF;
  }
  return 1;
}

int
ssvm::
max_match(int s, int t)
{
  int flow = 0;
  for(;;) {  
    int f = BellmanFord(s, t);
    if(f == -1) {return flow;} // no_more_path
    flow += 1;
  }
}

void
ssvm::
bm(string pseq, string rseq)
{
  for(int i = 1; i < Pro_len + 1; i++) // start(0) - protein(1~P)
    add_edge(0, i, 0.00000000001);
  for(int i = 1; i < RNA_len + 1; i++) // RNA(P+1~P+R) - goal(P+R+1)
    add_edge(Pro_len + i, Pro_len + RNA_len + 1, 0.00000000001);
  
  int m = max_match(0, Pro_len + RNA_len + 1);
  if(m == 0) {
    if(output>0){ cout << "no matching" << endl;}
    result++;
  }
  if(m > 0) {
    if(output>0){ cout << "match = " << m << endl;}
  
    for(int i = 0; i < Pro_len; i++) { // cost(1 ~ P+R)
      for(int j = 0; j < RNA_len; j++) {
	if(-INF + H_bonus < cost[j + Pro_len + 1][i + 1]) {
	  if(output>0){ printf("%c %d - %d %c\n", protein[i+2], i, j, RNA[j+2]);}
	  predict_edge++;
	  if(unlabel == 1){un_predict_edge++;}
	  if(correct_match[i + 1][j + Pro_len + 1] == 0) {
	    result++; 
	    if(correct_protein[i] > 0){correct_predict_protein++;}
	    if(correct_rna[j] > 0){correct_predict_rna++;}
	  }
	  if(correct_match[i + 1][j + Pro_len + 1] > 0) {
	    correct_predict_protein++; correct_predict_rna++; correct_predict++;
	  }
	  
	  // RNAss
	  if(f_Rss5>0){ 
	    sprintf(cRss5, "%c%c%c%c%c", rna_ss[j], rna_ss[j+1], rna_ss[j+2], rna_ss[j+3], rna_ss[j+4]);
	    mRss5[cRss5]++;
	    sprintf(cRss5_r, "%c%c%c%c%c", rna_ss[j+4], rna_ss[j+3], rna_ss[j+2], rna_ss[j+1], rna_ss[j]);
	    if(memcmp(cRss5, cRss5_r, 5) != 0){mRss5[cRss5_r]++;}	  
	    if(unlabel == 1){  // unlabeled
	      u_mRss5.insert(pair<string, double>(cRss5, 0));
	      u_rate_mRss5.insert(pair<string, double>(cRss5, 0));
	      u_mRss5[cRss5]++;
	      if(memcmp(cRss5, cRss5_r, 5) != 0){
		u_mRss5.insert(pair<string, double>(cRss5_r, 0));
		u_rate_mRss5.insert(pair<string, double>(cRss5_r, 0));
		u_mRss5[cRss5_r]++;}}
	  }
	  // Proteinss
	  if(f_Pss5>0){ 
	    sprintf(cPss5, "%c%c%c%c%c", pro_ss[i], pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], pro_ss[i+4]);	   
	    mPss5[cPss5]++;
	    sprintf(cPss5_r, "%c%c%c%c%c", pro_ss[i+4], pro_ss[i+3], pro_ss[i+2], pro_ss[i+1], pro_ss[i]);   
	    if(memcmp(cPss5, cPss5_r, 5) != 0){mPss5[cPss5_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mPss5.insert(pair<string, double>(cPss5, 0));
	      u_rate_mPss5.insert(pair<string, double>(cPss5, 0));
	      u_mPss5[cPss5]++;
	      if(memcmp(cPss5, cPss5_r, 5) != 0){
		u_mPss5.insert(pair<string, double>(cPss5_r, 0));
		u_rate_mPss5.insert(pair<string, double>(cPss5_r, 0));
		u_mPss5[cPss5_r]++;}}
	  }
	  // Pro-RNA ss
	  if(f_PsRs>0){ 
	    sprintf(cPsRs, "%c%c", pro_ss[i+2], rna_ss[j+2]);	   
	    mPsRs[cPsRs]++;
	    if(unlabel == 1){  // unlabeled
	      u_mPsRs.insert(pair<string, double>(cPsRs, 0));
	      u_rate_mPsRs.insert(pair<string, double>(cPsRs, 0));
	      u_mPsRs[cPsRs]++;}
	  }
	  // Pro-RNA ss 3-3
	  if(f_Ps3Rs3>0){ 
	    sprintf(cPs3Rs3, "%c%c%c%c%c%c", pro_ss[i+1], pro_ss[i+2], pro_ss[i+3], rna_ss[j+1], rna_ss[j+2], rna_ss[j+3]);
	    mPs3Rs3[cPs3Rs3]++;
	    sprintf(cPs3Rs3_r, "%c%c%c%c%c%c", pro_ss[i+3], pro_ss[i+2], pro_ss[i+1], rna_ss[j+3], rna_ss[j+2], rna_ss[j+1]);
	    if(memcmp(cPs3Rs3, cPs3Rs3_r, 6) != 0){mPs3Rs3[cPs3Rs3_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mPs3Rs3.insert(pair<string, double>(cPs3Rs3, 0));
	      u_rate_mPs3Rs3.insert(pair<string, double>(cPs3Rs3, 0));
	      u_mPs3Rs3[cPs3Rs3]++;
	      if(memcmp(cPs3Rs3, cPs3Rs3_r, 6) != 0){
		u_mPs3Rs3.insert(pair<string, double>(cPs3Rs3_r, 0));
		u_rate_mPs3Rs3.insert(pair<string, double>(cPs3Rs3_r, 0));
		u_mPs3Rs3[cPs3Rs3_r]++;}}
	  }
	  // AAp
	  if(f_AAp>0){ 
	    int Apn = AA_7group(protein[i+2]);
	    mAAp[Apn]++;
	    if(unlabel == 1){  // unlabeled
	      u_mAAp.insert(pair<int, double>(Apn, 0));
	      u_rate_mAAp.insert(pair<int, double>(Apn, 0));
	      u_mAAp[Apn]++;}
	  }
	  // AAab
	  if(f_AAab>0){ 
	    int Aabn = AA_3ab(protein[i+2]);
	    mAAab[Aabn]++;
	    if(unlabel == 1){  // unlabeled
	      u_mAAab.insert(pair<int, double>(Aabn, 0));
	      u_rate_mAAab.insert(pair<int, double>(Aabn, 0));
	      u_mAAab[Aabn]++;}
	  }
	  // AAh
	  if(f_AAh>0){ 
	    int Ahn = AA_3hydro(protein[i+2]);
	    mAAh[Ahn]++;
	    if(unlabel == 1){  // unlabeled
	      u_mAAh.insert(pair<int, double>(Ahn, 0));
	      u_rate_mAAh.insert(pair<int, double>(Ahn, 0));
	      u_mAAh[Ahn]++;}
	  }
	  // Rpp
	  if(f_Rpp>0){ 
	    int Rpn = RNA_2pp(RNA[j+2]);
	    mRpp[Rpn]++;
	    if(unlabel == 1){  // unlabeled
	      u_mRpp.insert(pair<int, double>(Rpn, 0));
	      u_rate_mRpp.insert(pair<int, double>(Rpn, 0));
	      u_mRpp[Rpn]++;}
	  }
	  // Protein 1 -
	  if(f_P1>0){ 
	    mP1[protein[i+2]]++;
	    if(unlabel == 1){  // unlabeled
	      u_mP1.insert(pair<char, double>(protein[i+2], 0));
	      u_rate_mP1.insert(pair<char, double>(protein[i+2], 0));
	      u_mP1[protein[i+2]]++;}
	  }
	  // Protein 1 - RNA 1
	  if(f_P1R1>0){ 
	    sprintf(cP1R1, "%c%c", protein[i+2], RNA[j+2]);
	    mP1R1[cP1R1]++;
	    if(unlabel == 1){  // unlabeled
	      u_mP1R1.insert(pair<string, double>(cP1R1, 0));
	      u_rate_mP1R1.insert(pair<string, double>(cP1R1, 0));
	      u_mP1R1[cP1R1]++;}
	  }
	  // Protein 3 - RNA 1
	  if(f_P3R1>0){ 
	    sprintf(cP3R1, "%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], RNA[j+2]);
	    mP3R1[cP3R1]++;
	    sprintf(cP3R1_r, "%c%c%c%c", protein[i+3], protein[i+2], protein[i+1], RNA[j+2]);
	    if(memcmp(cP3R1, cP3R1_r, 4) != 0){mP3R1[cP3R1_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mP3R1.insert(pair<string, double>(cP3R1, 0));
	      u_rate_mP3R1.insert(pair<string, double>(cP3R1, 0));	   
	      u_mP3R1[cP3R1]++;
	      if(memcmp(cP3R1, cP3R1_r, 4) != 0){
		u_mP3R1.insert(pair<string, double>(cP3R1_r, 0));
		u_rate_mP3R1.insert(pair<string, double>(cP3R1_r, 0));
		u_mP3R1[cP3R1_r]++;}}
	  }
	  // Protein 3 - RNA 3
	  if(f_P3R3>0){
	    sprintf(cP3R3, "%c%c%c%c%c%c", protein[i+1], protein[i+2], protein[i+3], RNA[j+1], RNA[j+2], RNA[j+3]);
	    mP3R3[cP3R3]++;
	    sprintf(cP3R3_r, "%c%c%c%c%c%c", protein[i+3], protein[i+2], protein[i+1], RNA[j+3], RNA[j+2], RNA[j+1]);
	    if(memcmp(cP3R3, cP3R3_r, 6) != 0){mP3R3[cP3R3_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mP3R3.insert(pair<string, double>(cP3R3, 0));
	      u_rate_mP3R3.insert(pair<string, double>(cP3R3, 0));
	      u_mP3R3[cP3R3]++;
	      if(memcmp(cP3R3, cP3R3_r, 6) != 0){
		u_mP3R3.insert(pair<string, double>(cP3R3_r, 0));
		u_rate_mP3R3.insert(pair<string, double>(cP3R3_r, 0));
		u_mP3R3[cP3R3_r]++;}}	 
	  }
	  // Protein 5 - RNA 1
	  if(f_P5R1>0){	 
	    sprintf(cP5R1, "%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], RNA[j+2]);
	    mP5R1[cP5R1]++;
	    sprintf(cP5R1_r, "%c%c%c%c%c%c", protein[i+4], protein[i+3], protein[i+2], protein[i+1], protein[i], RNA[j+2]);
	    if(memcmp(cP5R1, cP5R1_r, 6) != 0){mP5R1[cP5R1_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mP5R1.insert(pair<string, double>(cP5R1, 0));
	      u_rate_mP5R1.insert(pair<string, double>(cP5R1, 0));
	      u_mP5R1[cP5R1]++;
	      if(memcmp(cP5R1, cP5R1_r, 6) != 0){
		u_mP5R1.insert(pair<string, double>(cP5R1_r, 0));
		u_rate_mP5R1.insert(pair<string, double>(cP5R1_r, 0));
		u_mP5R1[cP5R1_r]++;}}	
	  }
	  // Protein 1 - RNA 5
	  if(f_P1R5>0){	 
	    sprintf(cP1R5, "%c%c%c%c%c%c", protein[i+2], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	    mP1R5[cP1R5]++;
	    sprintf(cP1R5_r, "%c%c%c%c%c%c", protein[i+2], RNA[j+4], RNA[j+3], RNA[j+2], RNA[j+1], RNA[j]);
	    if(memcmp(cP1R5, cP1R5_r, 6) != 0){mP1R5[cP1R5_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mP1R5.insert(pair<string, double>(cP1R5, 0));
	      u_rate_mP1R5.insert(pair<string, double>(cP1R5, 0));
	      u_mP1R5[cP1R5]++;
	      if(memcmp(cP1R5, cP1R5_r, 6) != 0){
		u_mP1R5.insert(pair<string, double>(cP1R5_r, 0));
		u_rate_mP1R5.insert(pair<string, double>(cP1R5_r, 0));	     
		u_mP1R5[cP1R5_r]++;}}
	  }
	  // Protein 5 - RNA 5
	  if(f_P5R5>0){	 
	    sprintf(cP5R5, "%c%c%c%c%c%c%c%c%c%c", protein[i], protein[i+1], protein[i+2], protein[i+3], protein[i+4], RNA[j], RNA[j+1], RNA[j+2], RNA[j+3], RNA[j+4]);
	    mP5R5[cP5R5]++;
	    sprintf(cP5R5_r, "%c%c%c%c%c%c%c%c%c%c", protein[i+4], protein[i+3], protein[i+2], protein[i+1], protein[i], RNA[j+4], RNA[j+3], RNA[j+2], RNA[j+1], RNA[j]);
	    if(memcmp(cP5R5, cP5R5_r, 10) != 0){mP5R5[cP5R5_r]++;}
	    if(unlabel == 1){  // unlabeled
	      u_mP5R5.insert(pair<string, double>(cP5R5, 0));
	      u_rate_mP5R5.insert(pair<string, double>(cP5R5, 0));	    
	      u_mP5R5[cP5R5]++;
	      if(memcmp(cP5R5, cP5R5_r, 10) != 0){
		u_mP5R5.insert(pair<string, double>(cP5R5_r, 0));
		u_rate_mP5R5.insert(pair<string, double>(cP5R5_r, 0));
		u_mP5R5[cP5R5_r]++;}}	  
	  }
	}}}
  }
  if(output>0){ cout << endl;}
}

int
ssvm::
get_pro_seq(string seq) {
  
  string seqfile;
  seqfile = seq;
  seqfile = seqfile.append(".seq");
  protein[0] = '-'; protein[1] = '-';
  char p_char;
  int p_num = 1;
  FILE *fp = fopen(seqfile.c_str(), "r");
  while(!feof(fp)) {
    p_char = fgetc(fp);
    if(p_char == EOF) {break;} // end of file
    if(p_num == 1) { // ">NAME" line
      if(p_char == '\n') {p_num++;continue;}}
    if(p_num > 1) { // seq lines
      if(p_char == '\n') {continue;} // line end
      if(p_char != '\n'){
	protein[p_num] = p_char;// protein[2]~[P+1]
	p_num++;}}
  }
  fclose(fp);
  protein[p_num] = '-'; protein[p_num+1] = '-';

  /*
    string pro2nd_pdb = seq;
    pro2nd_pdb = pro2nd_pdb.append(".2nd");
    cout << pro2nd_pdb << endl;
    get_pro_2nd_pdb(pro2nd_pdb);
  */
  string pro2nd_ss2 = seq;
  pro2nd_ss2 = pro2nd_ss2.append(".ss2");
  if(output>0){ cout << pro2nd_ss2 << endl;}
  get_pro_2nd_ss2(pro2nd_ss2);
   
  return p_num-2;
}

int
ssvm::
get_rna_seq(string seqfile) {

  RNA[0] = '-';  RNA[1] = '-';
  char r_char;
  int r_num = 1;
  FILE *fp = fopen(seqfile.c_str(), "r");
  while(!feof(fp)) {
    r_char = fgetc(fp);
    if(r_char == EOF) {break;} // end of file
    if(r_num == 1) { // ">NAME" line
      if(r_char == '\n') {r_num++; continue;}}
    if(r_num > 1) { // seq lines
      if(r_char == '\n') {continue;} // line end
      if(r_char != '\n'){
	RNA[r_num] = r_char;// RNA[2]~[R+1]
	r_num++;}}
  }
  fclose(fp);
  RNA[r_num] = '-'; RNA[r_num+1] = '-';
  
  rna_ss[0] = '-'; rna_ss[1] = '-';
  FILE *fpR;
  char rss_char;
  int rss_num = 1;
  const char *cmdline = "/home/slab/kashiwagi/centroid_fold-0.0.9-x86_64-linux/centroid_fold"; 
  char cmd[1000];
  sprintf(cmd, "%s %s" ,cmdline, seqfile.c_str());
  fpR = popen(cmd,"r");
  while(!feof(fpR)) {
    rss_char = fgetc(fpR);
    if(rss_char == EOF) {break;}
    if(rss_num == 1) {
      if(rss_char == '\n') {rss_num++; continue;}
    }
    if(rss_num >1) {
      if(rss_char == '(' || rss_char == ')' || rss_char == '.') {
	rna_ss[rss_num] = rss_char;
	rss_num++;}
      if(rss_char == ' '){break;}
    }
  }
  pclose(fpR);
  rna_ss[rss_num] = '-'; rna_ss[rss_num+1] = '-';
  if(output>0){ cout << "\n ------------ predicted RNA structure by CentroidFold" << endl;}
  for(int i = 0; i < rss_num-2; i++){
    if(output>0){ cout << RNA[i+2];}
  }
  if(output>0){ cout << endl;}

  for(int i = 0; i < rss_num-2; i++){
    if(output>0){ cout << rna_ss[i+2];}
    if(rna_ss[i+2] == '(' || rna_ss[i+2] == ')') {
      rna_ss[i+2] = '|';
    }
  }
  if(output>0){ cout << "\n\nRNA length = " << rss_num-2 << endl << endl;}
  return r_num-2;
}

int
ssvm::
get_pro_2nd_pdb(string seq) {
  char pss_char;
  int pss_num = 1;
  pro_ss[0] = '-'; pro_ss[1] = '-';
  FILE *fp = fopen(seq.c_str(), "r");
  while(!feof(fp)) {
    pss_char = fgetc(fp);
    if(pss_char == EOF) {break;} // end of file
    if(pss_num == 1) { // ">NAME" line
      if(pss_char == '\n') {pss_num++; continue;}
    }
    if(pss_num > 1) { // ss lines
      if(pss_char == '\n') {continue;} // line end
      else{
	pro_ss[pss_num] = pss_char;
	pss_num++;}// pss[2]~[P+1]
    }
  }
  fclose(fp);
  pro_ss[pss_num] = '-'; pro_ss[pss_num+1] = '-';

  for(int i = 0; i < pss_num-2; i++){
    if(output>0){ cout << protein[i+2];}
  }
  if(output>0){ cout << endl;}
  for(int i = 0; i < pss_num-2; i++){
    if(output>0){ cout << pro_ss[i+2];}
  }
  if(output>0){ cout << endl;}
}

int
ssvm::
get_pro_2nd_ss2(string seq) {
  int pss_num = 0;
  string ss2_line;
  pro_ss[0] = '-'; pro_ss[1] = '-';
  ifstream ss2file(seq.c_str());
  while( !ss2file.eof() ) {
    getline(ss2file, ss2_line);
    const char *LINE = ss2_line.c_str(); 
    if(pss_num > 1 && ss2_line.size() > 3) {
      if(LINE[7] == 'C') {
	pro_ss[pss_num] =  '.';}
      if(LINE[7] == 'E') {
	pro_ss[pss_num] =  '>';}
      if(LINE[7] == 'H') {
	pro_ss[pss_num] =  '=';}
    }
    pss_num++;
  }
  pro_ss[pss_num-1] = '-'; pro_ss[pss_num] = '-';
  
  for(int i = 0; i < pss_num-3; i++){
    if(output>0){ cout << protein[i+2];}
  }
  if(output>0){ cout << endl;}
  for(int i = 0; i < pss_num-3; i++){
    if(output>0){ cout << pro_ss[i+2];}
  }
  if(output>0){ cout << endl;}
}

int
ssvm::
AA_7group(char a)
{
  if(a == 'A' || a == 'G' || a == 'V') {return 0;}
  if(a == 'I' || a == 'L' || a == 'F' || a == 'P') {return 1;}
  if(a == 'Y' || a == 'M' || a == 'T' || a == 'S') {return 2;}
  if(a == 'H' || a == 'N' || a == 'Q' || a == 'W') {return 3;}
  if(a == 'R' || a == 'K') {return 4;}
  if(a == 'D' || a == 'E') {return 5;}
  if(a == 'C') {return 6;}
}

int
ssvm::
AA_3ab(char a)
{
  if(a == 'D' || a == 'E') {return 0;}
  if(a == 'R' || a == 'H' || a == 'K') {return 1;}
  else{return 2;}
}

int
ssvm::
AA_3hydro(char a)
{
  if(a == 'A' || a == 'M' || a == 'C' || a == 'F' || a == 'L' || a == 'V' || a == 'I'){return 0;}
  if(a == 'R' || a == 'K' || a == 'Q' || a == 'N' || a == 'E' || a == 'D' || a == 'H'){return 1;}
  else{return 2;}
}

int
ssvm::
RNA_2pp(char b)
{
  if(b == 'A' || b == 'G'){return 0;}
  else{return 1;}
}

void
ssvm::
output_result()
{
  cout << "all_edge: " << all_edge << " Protein: " << all_protein << " RNA: " << all_RNA << endl;
  cout << "all_correct_edge: " << all_correct_edge << endl;
  cout << "predict_edge = " << predict_edge << endl;
  cout << "correct_predict_edge: " << correct_predict << endl;
  cout << "correct_predict_protein: " << correct_predict_protein << endl;
  cout << "correct_predict_rna: " << correct_predict_rna << endl << endl;
  double ALL = (double)all_edge; double ALLp = (double)all_protein; double ALLr = (double)all_RNA;
  double AC = (double)all_correct_edge;
  double pe = (double)predict_edge;
  double TPe = (double)correct_predict; double TPp = (double)correct_predict_protein; double TPr = (double)correct_predict_rna;
  double FPe = pe - TPe; double FPp = pe - TPp; double FPr = pe - TPr;
  double TNe = ALL - AC - FPe; double TNp = ALLp - AC - FPp; double TNr = ALLr - AC - FPr;
  double FNe = AC - TPe; double FNp = AC - TPp; double FNr = AC - TPr;
  cout << "edge_PPV = " << TPe/(TPe + FPe) << endl;
  cout << "edge_Sensitivity = " << TPe/(TPe + FNe) << endl;
  cout << "edge_Specificity = " << TNe/(FPe + TNe) << endl;
  cout << "Protein_PPV = " <<  TPp/(TPp + FPp) << endl;
  cout << "Protein_Sensitivity = " <<  TPp/(TPp + FNp) << endl;
  cout << "Protein_Specificity = " << TNp/(FPp + TNp) << endl;
  cout << "RNA_PPV = " <<  TPr/(TPr + FPr) << endl;
  cout << "RNA_Sensitivity = " <<  TPr/(TPr + FNr) << endl;
  cout << "RNA_Specificity = " << TNr/(FPr + TNr) << endl;
 
  PPV_edge += TPe/(TPe + FPe);
  PPV_Pro += TPp/(TPp + FPp);
  PPV_RNA += TPr/(TPr + FPr);
  SEN_edge += TPe/(TPe + FNe);
  SEN_Pro += TPp/(TPp + FNp);
  SEN_RNA += TPr/(TPr + FNr);
  SPE_edge += TNe/(FPe + TNe);
  SPE_Pro += TNp/(FPp + TNp);
  SPE_RNA += TNr/(FPr + TNr);
}

void
ssvm::
output_train_weight()
{
  ofstream t_out("train_output");
  int f_n = 0; // fearture_num

  // Rss5 // MAX 162
  if(f_Rss5>0){ 
    t_out << "Rss5" << endl;
    map<string, double>::iterator it_Rss5_w = mRss5_w.begin();
    while( it_Rss5_w != mRss5_w.end() ){
      double a = mRss5_w[(*it_Rss5_w).first], b = u_mRss5_w[(*it_Rss5_w).first], c = u_rate_mRss5[(*it_Rss5_w).first];
      t_out << (*it_Rss5_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_Rss5_w;
      f_n++;}
  }
  // Pss5 // MAX 768
  if(f_Pss5>0){ 
    t_out << "Pss5" <<endl;
    map<string, double>::iterator it_Pss5_w = mPss5_w.begin();
    while( it_Pss5_w != mPss5_w.end() ){
      double a = mPss5_w[(*it_Pss5_w).first], b = u_mPss5_w[(*it_Pss5_w).first], c = u_rate_mPss5[(*it_Pss5_w).first];
      t_out << (*it_Pss5_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_Pss5_w;
      f_n++;}
  }
  // Pro-RNA ss // MAX 12
  if(f_PsRs>0){
    t_out << "PsRs" <<endl;
    map<string, double>::iterator it_PsRs_w = mPsRs_w.begin();
    while( it_PsRs_w != mPsRs_w.end() ){
      double a = mPsRs_w[(*it_PsRs_w).first], b = u_mPsRs_w[(*it_PsRs_w).first], c = u_rate_mPsRs[(*it_PsRs_w).first];
      t_out << (*it_PsRs_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_PsRs_w;
      f_n++;}
  }
  // Pro-RNA ss 3-3 // MAX 1728
  if(f_Ps3Rs3>0){
    t_out << "Ps3Rs3" <<endl;
    map<string, double>::iterator it_Ps3Rs3_w = mPs3Rs3_w.begin();
    while( it_Ps3Rs3_w != mPs3Rs3_w.end() ){
      double a = mPs3Rs3_w[(*it_Ps3Rs3_w).first], b = u_mPs3Rs3_w[(*it_Ps3Rs3_w).first], c = u_rate_mPs3Rs3[(*it_Ps3Rs3_w).first];
      t_out << (*it_Ps3Rs3_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_Ps3Rs3_w;
      f_n++;}
  }
  // AAp // 7
  if(f_AAp>0){
    t_out << "AAp" << endl;
    for(int i = 0; i < 7; i++) {
      double a = mAAp_w[i], b = u_mAAp_w[i], c = u_rate_mAAp[i];
      t_out << i << " " << a << " + " <<  b << " * " << c << " = " << a+b*c << endl;
      f_n++;}
  }
  // AAab // 3
  if(f_AAab>0){
    t_out << "AAab" << endl;
    for(int i = 0; i < 3; i++) {
      double a = mAAab_w[i], b = u_mAAab_w[i], c = u_rate_mAAab[i];
      t_out << i << " " << a << " + " <<  b << " * " << c << " = " << a+b*c << endl;
      f_n++;}
  }
  // AAh // 3
  if(f_AAh>0){
    t_out << "AAh" << endl;
    for(int i = 0; i < 3; i++) {
      double a = mAAh_w[i], b = u_mAAh_w[i], c = u_rate_mAAh[i];
      t_out << i << " " << a << " + " <<  b << " * " << c << " = " << a+b*c << endl;
      f_n++;}
  }
  // Rpp // 2
  if(f_Rpp>0){
    t_out << "Rpp" << endl;
    for(int i = 0; i < 2; i++) {
      double a = mRpp_w[i], b = u_mRpp_w[i], c = u_rate_mRpp[i];
      t_out << i << " " << a << " + " <<  b << " * " << c << " = " << a+b*c << endl;
      f_n++;}
  }
  // Protein 1 - // 20
  if(f_P1>0){
    t_out << "P1" << endl;
    for(int i = 0; i < 20; i++) {
      double a = mP1_w[aa[i]], b = u_mP1_w[aa[i]], c = u_rate_mP1[aa[i]];
      t_out << aa[i] << " " << a << " + " <<  b << " * " << c << " = " << a+b*c << endl;
      f_n++;}
  }
  // Protein 1 - RNA 1 // MAX 80
  if(f_P1R1>0){
    t_out << "P1R1" <<endl;
    map<string, double>::iterator it_P1R1_w = mP1R1_w.begin();
    while( it_P1R1_w != mP1R1_w.end() ) {
      double a = mP1R1_w[(*it_P1R1_w).first], b = u_mP1R1_w[(*it_P1R1_w).first], c = u_rate_mP1R1[(*it_P1R1_w).first];
      t_out << (*it_P1R1_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P1R1_w;
      f_n++;}
  }
  // Protein 3 - RNA 1 // MAX 8820
  if(f_P3R1>0){
    t_out << "P3R1" <<endl;
    map<string, double>::iterator it_P3R1_w = mP3R1_w.begin();
    while( it_P3R1_w != mP3R1_w.end() ) {
      double a = mP3R1_w[(*it_P3R1_w).first], b = u_mP3R1_w[(*it_P3R1_w).first], c = u_rate_mP3R1[(*it_P3R1_w).first];
      t_out << (*it_P3R1_w).first << " " << a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P3R1_w;
      f_n++;}
  }
  // Protein 3 - RNA 3 // MAX 882,000
  if(f_P3R3>0){
    t_out << "P3R3" <<endl;
    map<string, double>::iterator it_P3R3_w = mP3R3_w.begin();
    while( it_P3R3_w != mP3R3_w.end() ) {
      double a = mP3R3_w[(*it_P3R3_w).first], b = u_mP3R3_w[(*it_P3R3_w).first], c = u_rate_mP3R3[(*it_P3R3_w).first];
      t_out << (*it_P3R3_w).first << " " <<  a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P3R3_w;
      f_n++;}
  }
  // Protein 5 - RNA 1 // MAX 15,558,480
  if(f_P5R1>0){ 
    t_out << "P5R1" <<endl;
    map<string, double>::iterator it_P5R1_w = mP5R1_w.begin();
    while( it_P5R1_w != mP5R1_w.end() ) {
      double a = mP5R1_w[(*it_P5R1_w).first], b = u_mP5R1_w[(*it_P5R1_w).first], c = u_rate_mP5R1[(*it_P5R1_w).first];
      t_out << (*it_P5R1_w).first << " " <<  a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P5R1_w;
      f_n++;}
  }
  // Protein 1 - RNA 5 // MAX 50000
  if(f_P1R5>0){ 
    t_out << "P1R5" <<endl;
    map<string, double>::iterator it_P1R5_w = mP1R5_w.begin();
    while( it_P1R5_w != mP1R5_w.end() ) {
      double a = mP1R5_w[(*it_P1R5_w).first], b = u_mP1R5_w[(*it_P1R5_w).first], c = u_rate_mP1R5[(*it_P1R5_w).first];
      t_out << (*it_P1R5_w).first << " " <<  a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P1R5_w;
      f_n++;}
  }
  // Protein 5 - RNA 5 // MAX 25,204,737,600
  if(f_P5R5>0){
    t_out << "P5R5" <<endl;
    map<string, double>::iterator it_P5R5_w = mP5R5_w.begin();
    while( it_P5R5_w != mP5R5_w.end() ) {
      double a = mP5R5_w[(*it_P5R5_w).first], b = u_mP5R5_w[(*it_P5R5_w).first], c = u_rate_mP5R5[(*it_P5R5_w).first];
      t_out << (*it_P5R5_w).first << " " <<  a << " + " << b << " * " << c << " = " << a+b*c << endl;
      ++it_P5R5_w;
      f_n++;}
  }

  cout << "feature nun = " << f_n << endl;
}

void
ssvm::
output_cost()
{
  for(int i=0; i<Pro_len; i++) {
    for(int j=0; j<RNA_len; j++) {
      cout << cost[i + 1][j + Pro_len + 1] << endl;
    }}
}

void
printUsage()
{
  cout << "Usage: [-option] labeled_data [unlabeled_data]" << endl;
  cout << "\t-b : Hamming bonus \tdefault = 0.01" << endl;
  cout << "\t-e : learning eta \tdefault = 0.01" << endl;
  cout << "\t-f : FOBOS clip \tdefault= 0.001" << endl;
  cout << "\t-r : threshold rate \tdefault = 0.9" << endl;
  cout << "\t-c : cross-validation \tdefault = 4" << endl;
  cout << "-o : output only result" << endl;
}

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
