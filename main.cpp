#include "Rankine.h"
#include "RankineHX.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "hdf5.h"

#define EPS 1.e-3
using namespace std;

#define HDF5_WRITE

int main() {
  // Checking the Rankine Cycle .. borrowed from Cengel and Boles
  // Thermodynamics book's solved problem on the simple ideal Rankine
  // Cycle
  // Rankine rn{"Water", 507, 364.91, 623.15, 1.01};
  // cout << rn.printStates();
  // cout << "Efficiency: " << rn.calc_eff() << endl;
  // cout << "Net Power: " << rn.calc_power() / 1000.0 << " kJ/kg" << endl;

  Rankine rn{"Ammonia", 285.15, 284.15, 288.15, 1.01};
  cout << rn.printStates();
  cout << "Efficiency: " << rn.calc_eff() << endl;
  cout << "Net Power: " << rn.calc_power() / 1000.0 << " kJ/kg" << endl;


  // For the following algorithm to work correctly, ensure that the
  // dTpinch should be a multiple of the smallest unit on the discrete
  // temperature scale. The same goes for the difference (Twwi -Tcwi).
  // Parametric study - Problem input
  
  const double Tref = 273.15;
  const double Twwi = 30.0 + Tref;  // K
  const double Tcwi = 5.0 + Tref;   // K

  const double dTpinch = 4.0;       // K or (deg. C)
  const string wwf{"Water"};
  const string cwf{"Water"};
  const string wf{"Ammonia"};

  const double m_ww = 1.0; // kg/s (mass flow rate of warm water)

  // increase n to 2, 4, 8 for refining the discrete temperature scale
  const int n = 1;

  int N = round(Twwi - Tcwi) * n; // used round for taking care of
                                  // round-off errors in calculating the
                                   // difference, The temp. difference (Twwi-Tcwi) should
                                   // be an integer multiple of dT.
  int Npinch = round(dTpinch) * n; // round-off errors in storing
                                   // dTpinch, dTpinch should be an integer multiple of dT.
  

  double dT = (Twwi - Tcwi) / N;

#ifndef HDF5_WRITE
  
  ofstream outClientFile{"ORC_parametric.txt", ios::out};
  if (!outClientFile) {
    cerr << "File could not be opened" << endl;
    exit(EXIT_FAILURE);
  }

  // set all sticky settings of printing
  outClientFile << right; // right justified for numbers (sticky)
  outClientFile << fixed; // fixed point notation as opposed to scientific (sticky)

  const int wd = 13;
  const int precl = 4;
  const int precs = 2;
  
  outClientFile << "Input parameters are:\n";
  outClientFile << setw(wd) << "Twwi = " << (Twwi - Tref)
                << " C, Tcwi = " << (Tcwi - Tref) << " C, m_ww = " << m_ww
                << " kg/s \n";
  outClientFile << "Working fluids are:\n";
  outClientFile << "\twarm HX fluid " << wwf << "\n";
  outClientFile << "\tcold HX fluid " << cwf << "\n";
  outClientFile << "\tworking fluid " << wf << "\n";

  outClientFile << "\n\n";

  outClientFile << setw(wd) << "Te (K)" << setw(wd) << "Twwo (K)" << setw(wd)
                << "Tc (K)" << setw(wd) << "Tcwo (K)" << setw(wd) << "T3 (K)"
                << setw(wd) << "x3" << setw(wd) << "UA_ww (W/K)" << setw(wd)
                << "UA_cw (W/K)" << setw(wd) << "Qww (W)" << setw(wd)
                << "Qcw (W)" << setw(wd) << "Pnet (W)" << setw(wd)
                << "Pturb (W)" << setw(wd) << "Ppump (W)" << setw(wd)
                << "m_wf (kg/s)" << setw(wd) << "m_cw (kg/s)" << setw(wd)
                << "eta (%)" << endl;
  // nested loop where the indices are mapped to variables as follows:
  // i: Te
  // j: Tww,o
  // k: Tc
  // l: Tcw,o
  // m: T3
  // on discrete temp. scale - index 0 representes Tcwi, & index N + 1
  // represents Twwi
  int count = 0;
  for (int i = 2 + Npinch; i <= N - Npinch; ++i) {
    double Te = Tcwi + i * dT;
    for (int j = i + Npinch; j <= N; ++j) {
      double Twwo = Tcwi + j * dT;
      for (int k = 1 + Npinch; k <= i - 1; ++k) {
        double Tc = Tcwi + k * dT;
        for (int l = 1; l <= k - Npinch; ++l) {
          double Tcwo = Tcwi + l * dT;
          for (int m = i; m <= N + 1 - Npinch; ++m) {
            double T3 = Tcwi + m * dT;
            double x3 = 1.0;
            if (T3 > Te + EPS) {
              x3 = 1.01;
            }
            RankineHX rnhx{wwf,  cwf, wf, Twwi, Twwo, Tcwi,
                           Tcwo, Te,  Tc, T3,   x3,   m_ww};
            outClientFile << setprecision(precs) << setw(wd) << Te << setw(wd)
                          << Twwo << setw(wd) << Tc << setw(wd) << Tcwo
                          << setw(wd) << T3 << setw(wd)
			  << setprecision(precl)
                          << x3
			  << setprecision(precs) << setw(wd)
                          << rnhx.get_UA_ww() << setw(wd) << rnhx.get_UA_cw()
                          << setw(wd) << rnhx.get_Qww() << setw(wd)
                          << rnhx.get_Qcw() << setw(wd)
			  << setprecision(precl)
                          << rnhx.get_Pnet() << setw(wd) << rnhx.calc_Pturb()
                          << setw(wd) << rnhx.calc_Ppump() << setw(wd)
                          << rnhx.get_mwf() << setw(wd) << rnhx.get_mcw()
                          << setw(wd) << rnhx.get_eff() * 100 << "\n";
            ++count; if (count == 1000) {outClientFile << flush; exit(EXIT_SUCCESS);}; // this line should be commented after testing
          } // endif m
        } // endif l
      } // endif k
    } // endif j
  } // endif i

#endif

  // outClientFile automatically gets closed

#ifdef HDF5_WRITE
  vector<double> arTe, arTwwo, arTc, arTcwo, arT3, arx3;
  vector<double> arUAww, arUAcw, arQww, arQcw;
  vector<double> arPnet, arPturb, arPpump, armwf, armcw, areff;

  // nested loop where the indices are mapped to variables as follows:
  // i: Te
  // j: Tww,o
  // k: Tc
  // l: Tcw,o
  // m: T3
  // on discrete temp. scale - index 0 representes Tcwi, & index N + 1
  // represents Twwi
  int count = 0;
  for (int i = 2 + Npinch; i <= N - Npinch; ++i) {
    double Te = Tcwi + i * dT;
    for (int j = i + Npinch; j <= N; ++j) {
      double Twwo = Tcwi + j * dT;
      for (int k = 1 + Npinch; k <= i - 1; ++k) {
        double Tc = Tcwi + k * dT;
        for (int l = 1; l <= k - Npinch; ++l) {
          double Tcwo = Tcwi + l * dT;
          for (int m = i; m <= N + 1 - Npinch; ++m) {
            double T3 = Tcwi + m * dT;
            double x3 = 1.0;
            if (T3 > Te + EPS) {
              x3 = 1.01;
            }
            RankineHX rnhx{wwf,  cwf, wf, Twwi, Twwo, Tcwi,
                           Tcwo, Te,  Tc, T3,   x3,   m_ww};
	    
            arTe.push_back(Te);
	    arTwwo.push_back(Twwo);
	    arTc.push_back(Tc);
	    arTcwo.push_back(Tcwo);
	    arT3.push_back(T3);
	    arx3.push_back(x3);

	    arUAww.push_back(rnhx.get_UA_ww());
	    arUAcw.push_back(rnhx.get_UA_cw());
	    arQww.push_back(rnhx.get_Qww());
	    arQcw.push_back(rnhx.get_Qcw());

	    arPnet.push_back(rnhx.get_Pnet());
	    arPturb.push_back(rnhx.calc_Pturb());
	    arPpump.push_back(rnhx.calc_Ppump());
	    armwf.push_back(rnhx.get_mwf());
	    armcw.push_back(rnhx.get_mcw());
	    areff.push_back(rnhx.get_eff() * 100);
	    
            // ++count; if (count % 1000 == 0) {cout << "count is: " << count << " " << flush;}
	    ++count; if (count == 1000) {goto out_loop;}
          } // endif m
        } // endif l
      } // endif k
    } // endif j
  } // endif i

 out_loop:
  
  cout << "Total number of rows at the end are: " << count << endl;

    // variable declarations
    hid_t file_id, dataspace_id, dataset_id, attr_space_id, attr_id, attr_type_id;
    herr_t status;
    hsize_t dims[1] = {static_cast<hsize_t>(count)};
    hsize_t attr_dims[1] = {1};

    // file, dataspace & attribute space ids
    file_id = H5Fcreate("ORC_Ammonia.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    dataspace_id = H5Screate_simple(1, dims, NULL);
    attr_space_id = H5Screate_simple(1, attr_dims, NULL);

    // string attribute data type and id (variable size)
    attr_type_id = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(attr_type_id, H5T_VARIABLE);

    // file attribute labmda
    auto setAttribute = [&attr_id, &status, attr_space_id] (auto dest_id, hid_t atype_id, const char* attr_name, auto aval) {
      attr_id = H5Acreate(dest_id, attr_name, atype_id, attr_space_id, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Awrite(attr_id, atype_id, &aval);
      status = H5Aclose(attr_id);
    };

    setAttribute(file_id, H5T_NATIVE_DOUBLE, "Twwi (K)", Twwi);
    setAttribute(file_id, H5T_NATIVE_DOUBLE, "Tcwi (K)", Tcwi);
    setAttribute(file_id, H5T_NATIVE_DOUBLE, "dTpinch (K)", dTpinch);
    setAttribute(file_id, H5T_NATIVE_DOUBLE, "m_ww (kg/s)", m_ww);
    const char *wf_ptr = wf.c_str();
    setAttribute(file_id, attr_type_id, "Working fluid", wf_ptr);

    auto setDatasetWithAttribute = [&dataset_id, &status, dataspace_id, file_id, setAttribute] (const char* dname, hid_t dtype_id, auto ddata, hid_t atype_id, auto unitval) {
      dataset_id = H5Dcreate(file_id, dname, dtype_id, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      // dataset Te attribute
      setAttribute(dataset_id, atype_id, "Unit", unitval);

      status = H5Dwrite(dataset_id, dtype_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, ddata);
      status = H5Dclose(dataset_id);      
    };

    setDatasetWithAttribute("/Te", H5T_NATIVE_DOUBLE, arTe.data(), attr_type_id, "K");
    setDatasetWithAttribute("/Twwo", H5T_NATIVE_DOUBLE, arTwwo.data(), attr_type_id, "K");
    setDatasetWithAttribute("/Tc", H5T_NATIVE_DOUBLE, arTc.data(), attr_type_id, "K");
    setDatasetWithAttribute("/Tcwo", H5T_NATIVE_DOUBLE, arTcwo.data(), attr_type_id, "K");
    setDatasetWithAttribute("/T3", H5T_NATIVE_DOUBLE, arT3.data(), attr_type_id, "K");
    setDatasetWithAttribute("/x3", H5T_NATIVE_DOUBLE, arx3.data(), attr_type_id, " ");

    setDatasetWithAttribute("/UAww", H5T_NATIVE_DOUBLE, arUAww.data(), attr_type_id, "W/K");
    setDatasetWithAttribute("/UAcw", H5T_NATIVE_DOUBLE, arUAcw.data(), attr_type_id, "W/K");
    setDatasetWithAttribute("/Qww", H5T_NATIVE_DOUBLE, arQww.data(), attr_type_id, "W");
    setDatasetWithAttribute("/Qcw", H5T_NATIVE_DOUBLE, arQcw.data(), attr_type_id, "W");

    setDatasetWithAttribute("/Pnet", H5T_NATIVE_DOUBLE, arPnet.data(), attr_type_id, "W");
    setDatasetWithAttribute("/Pturb", H5T_NATIVE_DOUBLE, arPturb.data(), attr_type_id, "W");
    setDatasetWithAttribute("/Ppump", H5T_NATIVE_DOUBLE, arPpump.data(), attr_type_id, "W");
    setDatasetWithAttribute("/mwf", H5T_NATIVE_DOUBLE, armwf.data(), attr_type_id, "kg/s");
    setDatasetWithAttribute("/mcw", H5T_NATIVE_DOUBLE, armcw.data(), attr_type_id, "kg/s");
    setDatasetWithAttribute("/Eff", H5T_NATIVE_DOUBLE, areff.data(), attr_type_id, "%");
    
    status = H5Tclose(attr_type_id);
    status = H5Sclose(attr_space_id);
    status = H5Sclose(dataspace_id);
    status = H5Fclose(file_id);

#endif
}
