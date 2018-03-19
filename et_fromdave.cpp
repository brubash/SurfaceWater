/*  Copyright 2017 Lambert Rubash

    This file is part of TopNetCpp, a translation and enhancement of
    Fortran TopNet.

    TopNetCpp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    TopNetCpp is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with TopNetCpp.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "topnet.hh"
#include <dirent.h>

//   subroutines to compute Evapotranspiration following methods in
//   Handbook of Hydrology, 1993, Maidment editor, Chapter 4 by
//   Shuttleworth on ET.

using namespace std;
using namespace Eigen;

int et(const int iyear,  const int month,     const int iday,      const double hour1,
	const double dt,     const double albedo, const double elevsb, const double elevtg,
    const double rlapse, const double temp,   const double dewp,   const double trange,
    const double xlat,   const double xlong,  const double stdlon, const ArrayXd &dtbar,
    double &ret,         const double tmin,   const double tmax,   const double wind2m, const int method);


// **************************** etall() ****************************
int etall(const double xlat, const double xlong, const double stdlon, const double elevtg,
	const ArrayXd &dtbar, double &evap, const double temp, const double dewp, const double trange,
	const double elevsb, const double albedo, const double rlapse, const int sdate, const int shour,
	const int dtsec, const int m, const int istep, int &iyear, int &month, int &iday,
	int &ihr, int &imm, int &isec, double &hour1, const double tmin, const double tmax,
	const double wind2m, const int method)
{
	double dt;
	//  tmin, tmax: Daily minimum temperature and max temperature in deg C
	//  wind2m:  Wind speed at 2 m
	//  method: A flag to indicate the reference ET method
	//     method=0 refers to the Priestly Taylor Method
	//     method=1 refers to Penman Monteith ASCE Standard Reference for short grass h=0.12 m known as ETo
	//     method=2 refers to Penman Monteith ASCE Standard Reference for tall alfalfa h=0.5 m known as ETr
#if TRACE
	static int ncalls = 0;
    string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> etall(" << ncalls << ")" << std::endl;
    }
    caller = "etall";
#endif
	if (istep == 1) {
		iyear = sdate/10000;
		month = sdate/100 - iyear*100;
		iday  = sdate - iyear*10000 - month*100;
		ihr   = (shour)/10000;
		imm   = shour/100 - ihr*100;
		isec  = shour - ihr*10000 - imm*100;
		hour1 = (double)(ihr) + (double)(imm)/60.0 + (double)(isec)/3600.0;
	}
	dt = (double)(dtsec)/3600.0;   // Time step is in hours
	et(iyear, month, iday, hour1, dt, albedo, elevsb, elevtg, rlapse,
		temp, dewp, trange, xlat, xlong, stdlon, dtbar, evap,
		tmin, tmax, wind2m, method);
	// want to update time outside here, after the snow routine has run
	// call updatetime(iyear,month,iday,hour1,dt)
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving etall(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

// **************************** et() ****************************
//    This subroutine does one time step

int et(const int iyear,  const int month,     const int iday,      const double hour1,
	const double dt,     const double albedo, const double elevsb, const double elevtg,
    const double rlapse, const double temp,   const double dewp,   const double trange,
    const double xlat,   const double xlong,  const double stdlon, const ArrayXd &dtbar,
    double &ret,         const double tmin,   const double tmax,   const double wind2m, const int method)
{
	const double solcon = 4914.0;  // solar constant (4914 kJ/m^2/hr)
	const double ae = 0.34, be = -0.14;  	// Parameters in Brunt like net emissivity function
								//        See equation 4.2.8 in Handbook of hydrology
	const double ac = 1.0, bc = 0.0;  	// Parameters in cloudiness factor equation
								//        for Humid areas.  See equation 4.2.10 in Handbook of hydrology
	const double sigma = 2.0747e-7; // Stefan Boltzman constant  kJ/m^2/hr/K
	const double alpha = 1.26;  // Priestly - Taylor coefficient for humid climates
	const double bca = 0.8, bcc = 2.4;  // Bristow Campbell Parameters
	const double cn1 = 37, cd1 = 0.34;  // ASCE ETo constants
	const double cn2 = 66, cd2 = 0.38;  // ASCE ETr constants

	double shift, sn, t, dewpt;
	double hri, coszen, tf;
    double ed, emnet, cf, xln, es, del, p;
    double xlv, Gamma, num1, num2, den;
	double tmaxl, tminl, vpd;

	// double tmin, tmax, Daily minimum temperature and max temperature in deg C
	// double wind2m, Wind speed at 2 m
	// int method, A flag to indicate the reference ET method
	//     method=0 refers to the Priestly Taylor Method
	//     method=1 refers to Penman Monteith ASCE Standard Reference for short grass h=0.12 m known as ETo
	//     method=2 refers to Penman Monteith ASCE Standard Reference for tall alfalfa h=0.5 m known as ETr

	//   hour1 = start hour of time step
	//   dt = time step in hours
	//   albedo = surface reflectivity coefficient
	//   elevsb = elevation of the basin (m)
	//   elevtg = elevation of temperature guage (m)
	//   rlapse = lapse rate in degrees C per meter (positive for decreases with height)
	//   xlat = basin lattitude in degrees, negative for southern hemisphere
	//   xlong = basin longitude in degrees, negative for western hemisphere
	//   stdlon = reference time longitude, negative for western hemisphere
	//   dtbar = 12 dimensional array of mean monthly temperature ranges
	//   bca = A parameter in Bristow Campbell
	//   bcc = C parameter in Bristow Campbell

	//    shift time for std longitude versus basin longitude
#if TRACE
	static int ncalls = 0;
    string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> et(" << ncalls << ")" << std::endl;
    }
    caller = "et";
#endif
	shift = (xlong - stdlon)/15.0;
	hyri(iyear, month, iday, hour1 + shift, dt, 0.0, 0.0, xlat, hri, coszen);
	atf(tf, trange, month, dtbar, bca, bcc);
	sn = (1.0 - albedo)*tf*hri*solcon;
	t = temp - (elevsb - elevtg)*rlapse;
	tmaxl = tmax - (elevsb - elevtg)*rlapse;
	tminl = tmin - (elevsb - elevtg)*rlapse;
	dewpt = min(dewp, t); // Dewpoint to use is the minimum of lapsed temperature
	//   and recorded dewpoint.  This assumes dewpoint (and therefore vapour pressure)
	//   is roughly spatially constant, and therefore relative humidity increases as
	//   temperature drops
	ed = svp(dewpt)/1000.0;  // Saturation vapour pressure at dew point  (kPa)
	emnet = ae + be*sqrt(ed);  // Net emmissivity for longwave radiation
	cf = ac * (tf/bca) + bc;  // Cloudiness factor based on equation 4.2.10
	es = svp(t)/1000.0;   //  Saturation vapour pressure at actual temperature (kPa)
	del = 4096.0*es/pow(237.3 + t, 2);   // Equation 4.2.3
	p = 101.3*pow((293.0 - 0.0065*elevsb)/(293.0), 5.256);   // Equation 4.4.12
	//           Atmospheric pressure kPa
	xlv = 2501.0 - 2.361*t;   // Eqn 4.2.1 for latent heat of vaporization kJ/kg
	Gamma = 1.6286*p/xlv;    //  Psychometric constant eqn 4.2.28 kPa/C
	if (method == 1 || method == 2) {
		//        Penman Monteith ASCE Standard Reference
	    xln = -cf*emnet*sigma*(pow(tmaxl+273.15, 4) + pow(tminl + 273.15, 4))*0.5; // REF-ET Equation 26
		es = (svp(tmin)/1000.0 + svp(tmax)/1000.0)*0.5; //convert Pa to kPa
		vpd = max(es - ed,0.0); //kPa
		num1 = del*(sn + xln)/xlv;
		if (method == 1) { // short grass ETo
			num2 = cn1*Gamma*wind2m*vpd/(273.15 + t);
			den = del + Gamma*(1.0 + cd1*wind2m);
		} else {  // taller alfalfa ETr
			num2 = cn2*Gamma*wind2m*vpd/(273.15 + t);
			den = del + Gamma*(1.0 + cd2*wind2m);
		}
		ret = (num1 + num2)/den;
		ret = max(dt*ret, 0.0);
		//	  The units are formally kg/m^2 but this is the same as mm
	} else { // Priestly Taylor is the default if method is anything else
		xln = -cf*emnet*sigma*pow(t + 273.15, 4);  // Net Longwave (eqn 4.2.7)
		ret = alpha*del/(del + Gamma)*(sn + xln);  // Radiation based reference ET
		//       Above is energy units  kJ/m^2/hr
		ret = max(dt*ret /(xlv), 0.0);  //  Priestly Taylor
		cout << ret << " else " << endl;
		//	  The units are now formally kg/m^2 but this is the same as mm
	}
	//      Set ET = 0 when radiation is negative

#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving et(" << ncalls << ")" << "\n";
    }
    ncalls++;
#endif

	return 0;
}

//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;

// **************************** cliparam () ****************************
//     To read in the parameters needed for ET calculations
//
int cliParam(double *xlat, double *xlong, double &stdlon, double *elevtg, double **dtBar, int &ns_temper, int *temper_id)
{
	// 8-dec-2004 special nooksack version
    int i, j;
	string testStr;
	bool exist;
	struct dirent *dirp;
	DIR *dp;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> cliParam(" << ncalls << ")" << std::endl;
    }
	caller = "cliParam";
#endif

	//------------ clipar.dat -------------------
	exist = false;
	if ((dp = opendir(".")) == NULL) {
		cerr << "Can't open directory\n";
		exit(EXIT_FAILURE);
	}
	while ((dirp = readdir(dp)) != NULL) {
		testStr = dirp->d_name;
		if (testStr == "clipar.dat") {
			exist = true;
			closedir(dp);
			break;
		}
	}
	if (exist == false) {
		cerr << " No 'clipar.dat' file found\n";
		exit(EXIT_FAILURE);
	}

	ifstream cliparFile("clipar.dat");	// Fortran unit 88
	if (!cliparFile.is_open()) {
		cerr << "Failed to open clipar.dat\n";
		exit(EXIT_FAILURE);
	}

	//-------------------------------------------

	getline(cliparFile, testStr, '\n');		// skip a line
	getline(cliparFile, testStr, '\n');
	cliparFile >> ns_temper;
	getline(cliparFile, testStr, '\n');		// read the rest of the line
	cliparFile >> stdlon;					//  One std longitude for all locations
	getline(cliparFile, testStr, '\n');		// read the rest of the line
	getline(cliparFile, testStr, '\n');		// read the headings line
	for (i = 0; i < ns_temper; i++) {
		for (j = 0; j < 12; j++) {
			dtBar[j][i] = 10.0;   // default diurnal temperature range
		}
 		cliparFile >> temper_id[i] >> xlat[i] >> xlong[i] >> elevtg[i];
 		for (j = 0; j < 12; j++) {
 			cliparFile >> dtBar[j][i];
 		}
	}
	cliparFile.close();
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving cliParam(" << ncalls << ")" << "\n" << endl;;
    }
    ncalls++;
#endif

	return 0;
}
