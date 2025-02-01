// Heat source and Sink Heat Exchangers around the Ranking Cycle

#ifndef RANKINEHX_H
#define RANKINEHX_H

#include "Rankine.h"
#include <string>
using namespace std;

// Class definitions
class RankineHX {
public:
  explicit RankineHX(string wwf_i, string cwf_i, string wf_i, double Twwi_i,
                     double Twwo_i, double Tcwi_i, double Tcwo_i, double Te_i,
                     double Tc_i, double T3_i, double x3_i, double m_wwi = 1.0);

  double get_Pnet() const;   // fetches net power output of the Rankine cycle (object)
  double get_eff() const;    // fetches the efficiency of the Rankine cycle (object)
  
  double calc_Pturb() const; // Turbine power produced (W)
  double calc_Ppump() const; // Pump power consumed (W)
  double get_UA_ww() const;
  double get_UA_cw() const;
  double get_mwf() const;
  double get_mcw() const;
  double get_Qww() const;
  double get_Qcw() const;
  
private:
  // utility functions
  void calc_mwf();     // calculates mass flow rate of working fluid
  void calc_mcw();     // calculates mass flow rate of cold water
  double calc_Qww() const;   // Heat transferred at warm water HX
  double calc_Qcw() const;   // Heat transferred at cold water HX
  void calc_UA_ww();   // calculates UA_ww (see below) of warm water HX
  void calc_UA_cw();   // calculates UA_cw (see below) of cold water HX

  // data members
  string wwf;        // fluid for carrying warm water heat (seawater) - string
  string cwf;        // fluid for carrying cold water heat (seawater) - string
  string wf;         // working fluid of the Rankine cycle
  double Twwi, Twwo; // Warm water inlet and outlet temperature (K)
  double Tcwi, Tcwo; // Cold water inlet and outlet temperature (K)
  double Te;         // Evaporation (saturation) temperature (K)
  double Tc;         // Condenser (saturation) temperature (K)
  double T3;         // Temp. of evap. exit (K)
  double x3;         // Quality of evap. exit
  double m_ww;       // warm water flow rate (kg/s)
  double m_cw;       // cold water flow rate (kg/s)
  double m_wf;       // work fluid mass flow rate (kg/s)
  double Q_ww;       // Heat transferred at the warm water HX (W)
  double Q_cw;       // Heat transferred at the cold water HX (W)
  double UA_ww;      // (W/K) overall heat transfer coeff multiplied by HX Area
  double UA_cw;      // (W/K)

  const Rankine ranCycle;
};

#endif
