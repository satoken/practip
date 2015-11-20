#ifndef _INC_RNA_H_
#define _INC_RNA_H_

#include <string>
#include "typedef.h"

struct RNA {
  std::string name;
  std::string seq;
  std::string ss;
  std::string g2;
    
  RNA() {}
  RNA(const std::string& fa_name, const std::string& ss_name) : RNA() { read(fa_name, ss_name); }
  int read(const std::string& fa_name, const std::string& ss_name);
  int read_fa(const std::string& fa_name);
  int read_ss(const std::string& ss_name);

  static char group2(char);
  static void structural_profile(const std::string& ss, std::string& profile);

  static uint max_intraction(char x);
};

#endif

// Local Variables:
// mode: C++
// End:
