#ifndef _INC_AA_H__
#define _INC_AA_H__

#include <string>

struct AA {
  std::string name;
  std::string seq;
  std::string ss;
  std::string g10;
  std::string g8;
  std::string g4;
  std::string g2;
    
  AA() {}
  AA(const std::string& fa_name, const std::string& ss_name) : AA() { read(fa_name, ss_name); }
  int read(const std::string& fa_name, const std::string& ss_name);
  int read_fa(const std::string& fa_name);
  int read_ss(const std::string& ss_name);

  static char group10(char a);
  static char group8(char a);
  static char group4(char a);
  static char group2(char a);

  static uint max_intraction(char x);
};

#endif

// Local Variables:
// mode: C++
// End:
