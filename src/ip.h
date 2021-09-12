#ifndef __INC_IP_H__
#define __INC_IP_H__

class IPimpl;

class IP
{
public:
  typedef enum {MIN, MAX} DirType;
  typedef enum {FR, LO, UP, DB, FX} BoundType;

public:
  IP(DirType dir, int n_th);
  ~IP();
  int make_variable(double coef);
  int make_variable(double coef, int lo, int hi);
  int make_constraint(BoundType bnd, double l, double u);
  void add_constraint(int row, int col, double val);
  void update();
  double solve();
  double get_value(int col) const;

private:
  IPimpl* impl_;
};

#endif  // __INC_IP_H__

// Local Variables:
// mode: C++
// End:
