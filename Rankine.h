// Rankine Class Definition

#ifndef RANKINE_H
#define RANKINE_H

#include <string>
using namespace std;

// We are following CoolProp's units
// http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
// Go through the documentation of propsSI function if you cannot access the
// link above

// Pressure P: Pa
// Temperature T: K
// Enthalpy H: J/kg
// Entropy S: J/(kg-K)
//

// Rankine class definition
class Rankine {
public:
  explicit Rankine(string w_f, double T_e, double T_c, double T_3, double x_3);
  void setStates();
  string printStates() const;
  double calc_eff() const;
  double calc_power() const;
  double get_T(int istate) const;
  double get_x(int istate) const;
  double get_h(int istate) const;
  double get_s(int istate) const;

private:
  string wf;           // working fluid
  double Tevap, Pevap; // evaporator temp., press.
  double Tcond, Pcond; // condensor temp., press.
  double T3, x3;
  double T[5];
  double x[5];
  double h[5];
  double s[5];
};

#endif
