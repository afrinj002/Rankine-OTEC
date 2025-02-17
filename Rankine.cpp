#include "Rankine.h"
#include "CoolProp.h"
#include <iomanip>
#include <sstream>
#include <string>
#define EPS 1.e-3 // useful for relational operators, e.g., >, <, =
#define EPS2 1.e-5 // to check inequality of entropy
using namespace std;

// For coolprop units please refer to their documentation:
// http://www.coolprop.org/coolprop/HighLevelAPI.html#parameter-table
// Pressure P: Pa
// Temperature T: K
// Enthalpy H: J/kg
// Entropy S: J/(kg-K)
// Values with these units should be entered as arguments to CoolProp library
// functions

Rankine::Rankine(string w_f, double T_e, double T_c, double T_3, double x_3)
    : wf{w_f}, Tevap{T_e}, Tcond{T_c}, T3{T_3}, x3{x_3} {
  Pevap = CoolProp::PropsSI("P", "T", Tevap, "Q", 0.5, wf);
  Pcond = CoolProp::PropsSI("P", "T", Tcond, "Q", 0.5, wf);
  setStates();
}

void Rankine::setStates() {
  // State 1: Exit of Condensor (saturated liquid)
  T[1] = Tcond;
  x[1] = 0.0;
  h[1] = CoolProp::PropsSI("H", "T", T[1], "Q", x[1], wf);
  s[1] = CoolProp::PropsSI("S", "T", T[1], "Q", x[1], wf);

  // State 2: Exit of pump
  s[2] = s[1];
  x[2] = -0.01; // negative value represents subcooled liquid (mimicing EES)
  T[2] = CoolProp::PropsSI("T", "P", Pevap, "S", s[2], wf);
  h[2] = CoolProp::PropsSI("H", "P", Pevap, "S", s[2], wf);

  // State 3: Exit of evaporator
  T[3] = T3;
  x[3] = x3;
  if (T3 > Tevap + EPS) { // superheated state 3
    x[3] = 1.01;          // >1 value denotes superheated vapor (mimicing EES)
    h[3] = CoolProp::PropsSI("H", "T", T[3], "P", Pevap, wf);
    s[3] = CoolProp::PropsSI("S", "T", T[3], "P", Pevap, wf);
  } else { // mixture state 3
    h[3] = CoolProp::PropsSI("H", "P", Pevap, "Q", x3, wf);
    s[3] = CoolProp::PropsSI("S", "P", Pevap, "Q", x3, wf);
  }

  // State 4: Exit of Turbine
  s[4] = s[3];
  double s4sat = CoolProp::PropsSI("S", "P", Pcond, "Q", 1.0, wf);
  T[4] = CoolProp::PropsSI("T", "P", Pcond, "S", s[4], wf);
  h[4] = CoolProp::PropsSI("H", "P", Pcond, "S", s[4], wf);

  x[4] = 1.0;
  if (s[4] > s4sat + EPS2) { // superheated vapor state 4
    x[4] = 1.01;
  } else if (s[4] < s4sat - EPS2) { // mixture state 4
    x[4] = CoolProp::PropsSI("Q", "P", Pcond, "S", s[4], wf);
  }
}

string Rankine::printStates() const {
  ostringstream os;
  for (int i = 1; i < 5; ++i) {
    os << setfill(' ') << fixed << setprecision(2) << "State " << i << ": T "
       << setw(7) << T[i] << " K, x " << setw(5) << x[i] << ", h " << setw(7)
       << h[i] / 1000.0 << " kJ/kg, s " << setw(5) << s[i] / 1000.0
       << " kJ/(kg-K)" << endl;
  }
  return os.str();
}

double Rankine::calc_eff() const {
  double eff = ((h[3] - h[4]) - (h[2] - h[1])) / (h[3] - h[2]);
  return eff;
}

double Rankine::calc_power() const {
  // calculates net power output per unit mass flow rate of refrigerator (J/kg)
  double power = h[3] - h[4] - (h[2] - h[1]);
  return power;
}

double Rankine::get_T(int istate) const { return T[istate]; }

double Rankine::get_x(int istate) const { return x[istate]; }

double Rankine::get_h(int istate) const { return h[istate]; }

double Rankine::get_s(int istate) const { return s[istate]; }
