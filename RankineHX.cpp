// Implementation file of HX around heat source & heat sink

#include "RankineHX.h"
#include "CoolProp.h"
#include "Rankine.h"
#include <cmath>
#include <string>
const double patm = 1.e5;

RankineHX::RankineHX(string wwf_i, string cwf_i, string wf_i, double Twwi_i,
                     double Twwo_i, double Tcwi_i, double Tcwo_i, double Te_i,
                     double Tc_i, double T3_i, double x3_i, double m_wwi)

    : wwf{wwf_i}, cwf{cwf_i}, wf{wf_i}, Twwi{Twwi_i}, Twwo{Twwo_i},
      Tcwi{Tcwi_i}, Tcwo{Tcwo_i}, Te{Te_i}, Tc{Tc_i}, T3{T3_i}, x3{x3_i},
      m_ww{m_wwi}, ranCycle{wf, Te, Tc, T3, x3} {

  calc_mwf(); // calculates Q_ww as side-effect
  calc_mcw(); // calculates Q_cw as side-effect
  calc_UA_ww();
  calc_UA_cw();
}

void RankineHX::calc_mwf() {
  // find the enthalpies at inlet and outlet warm water
  Q_ww = calc_Qww();
  
  // enthalpies of working fluid at entry (2) and exit (3) of evaporator
  double h3 = ranCycle.get_h(3);
  double h2 = ranCycle.get_h(2);

  m_wf = Q_ww / (h3 - h2);
}

void RankineHX::calc_mcw() {
  // find the enthalpies at inlet and outlet cold water
  double h_cwi = CoolProp::PropsSI("H", "P", patm, "T", Tcwi, cwf);
  double h_cwo = CoolProp::PropsSI("H", "P", patm, "T", Tcwo, cwf);

  Q_cw = calc_Qcw();

  m_cw = Q_cw / (h_cwo - h_cwi);
}

void RankineHX::calc_UA_ww() {
  // first calculate LMTD at warm water HX
  double dT1 = Twwi - ranCycle.get_T(3);
  double dT2 = Twwo - ranCycle.get_T(2);

  double LMTD = (dT1 - dT2) / std::log(dT1 / dT2);

  UA_ww = Q_ww / LMTD;
}

void RankineHX::calc_UA_cw() {
  // calculate LMTD at cold water HX
  double dT1 = ranCycle.get_T(4) - Tcwo;
  double dT2 = ranCycle.get_T(1) - Tcwi;

  double LMTD = (dT1 - dT2) / std::log(dT1 / dT2);

  UA_cw = Q_cw / LMTD;
}

double RankineHX::calc_Qww() const {
  double h_wwi = CoolProp::PropsSI("H", "P", patm, "T", Twwi, wwf);
  double h_wwo = CoolProp::PropsSI("H", "P", patm, "T", Twwo, wwf);

  double Qww = m_ww * (h_wwi - h_wwo);
  return Qww;
}

double RankineHX::calc_Qcw() const {
  // enthalpies of working fluid at entry (4) and exit (1) of condenser
  double h4 = ranCycle.get_h(4);
  double h1 = ranCycle.get_h(1);

  double Qcw = m_wf * (h4 - h1);
  return Qcw;
}

double RankineHX::get_Pnet() const {
  double Pnet = m_wf * ranCycle.calc_power();
  return Pnet;
}

double RankineHX::get_eff() const { return ranCycle.calc_eff(); }

double RankineHX::calc_Pturb() const {
  double h4 = ranCycle.get_h(4);
  double h3 = ranCycle.get_h(3);

  double Pturb = m_wf * (h3 - h4);
  return Pturb;
}

double RankineHX::calc_Ppump() const {
  double h2 = ranCycle.get_h(2);
  double h1 = ranCycle.get_h(1);

  double Ppump = m_wf * (h2 - h1);
  return Ppump;
}

// get functions
double RankineHX::get_UA_ww() const { return UA_ww; }

double RankineHX::get_UA_cw() const { return UA_cw; }

double RankineHX::get_mwf() const { return m_wf; }

double RankineHX::get_mcw() const { return m_cw; }

double RankineHX::get_Qww() const { return Q_ww; }

double RankineHX::get_Qcw() const { return Q_cw; }
