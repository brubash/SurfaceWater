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
#include "snow.hh"

using namespace Eigen;
using namespace std;

//   Wrapper subroutine for calling by TOPNET
//   David Tarboton  1/10/05
//   This subroutine performs the following functions
//   -  subdivision of the external timestep into a number of internal time steps to represent diurnal forcing
//   -  Bypasses calling of snow subroutine when there is no snow on the ground or no snow falling

int snowueb(const ArrayXd &snowsitev, ArrayXd &snowstatev1, const ArrayXd &snowparam,
    const int ndepletionpoints, double **dfc, const Array<int,Dynamic,1> &snowcontrol, const ArrayXd &dtbar,
	const ArrayXd &snowforcing, ArrayXd &snowsurfacetemp1, ArrayXd &snowaveragetemp1,
    const double timestep, const int nstepday, double &surfacewaterinput, double &snowevaporation,  //outputs (both in m/h)
    double &areafractionsnow, const int modelelement)
{
	// real snowsitev(1)  sitev 1:9  forestcoverfrac, elevsb (m), xlat, xlon, stdlon,
	//      elevtg (m), rlapse (C/m), slope (deg), azimuth (deg)
	// real snowstatev(1)   1:5 point state variables (U, W, dimensionless age, refdepth, totalrefdepth)
	//        plus five depletion curve lumped parameters 6:11 (Wmax, areafrac, meltflag,Wtgt,Wtgtmax,Aftgt)
	// real snowparam(1)    25 snow parameters - definitions and default values as follows
	// (1) tr	  Temperature above which all is rain (3 C)
	// (2) ts	  Temperature below which all is snow (-1 C)
	// (3) es	  emmissivity of snow (nominally 0.99)
	// (4) cg 	  Ground heat capacity (nominally 2.09 KJ/kg/C)
	// (5) z	  Nominal meas. height for air temp. and humidity (2m)
	// (6) zo	  Surface aerodynamic roughness (0.01 m)
	// (7) rho	  Snow Density (Nominally200 kg/m^3)
	// (8) rhog	  Soil Density (nominally 1700 kg/m^3)
	// (9) lc	  Liquid holding capacity of snow (0.05)
	// (10) ks	  Snow Saturated hydraulic conductivity (200 m/hr)
	// (11) de	  Thermally active depth of soil (0.1 m)
	// (12)  abg	  Bare ground albedo  (0.25)
	// (13)  avo	  Visual new snow albedo (0.85)
	// (14)  anir0	  NIR new snow albedo (0.65)
	// (15) lans	  the thermal conductivity of fresh (dry) snow (0.33 kJ/m/k/hr)
	// (16) lang	  the thermal conductivity of soil (6.5 kJ/m/k/hr)
	// (17) wlf	  Low frequency fluctuation in deep snow/soil layer (1/4 w1)
	// (18) rd1	  Amplitude correction coefficient of heat conduction (1)
	// (19) fstab	  Stability correction control parameter 1 - means stability corrections fully used
	// (20) Tref	  Reference temperature of soil layer in ground heat calculation input
	// (21) dNewS	  The threshold depth of for new snow (0.002 m)
	// (22) gsurf	  The fraction of surface melt that runs off (e.g. from a glacier).
	// (23) df       Drift factor
	// (24) G        Ground heat flux (kJ/m2/h)
	// (25) AEF      Albedo extinction parameter

 	// real dfc(ndepletionpoints,2)   Depletion curve  wa/wamax, afrac
 	//      integer snowcontrol(1)   1:6 Control variables and flags
 	//   1 irad (0=estimated, 1=meas inc solar, 2=meas inc solar and longwave, 3=measured net radiation)
 	//   2 snowprintflag  0 do not print, 1 print
 	//   3 snow print unit
 	//   4 number of internal time steps in 1 calling time step
 	//   5 start date yyyymmdd
 	//   6 start time hhmmss

	//      real dtbar(12)   monthly mean diurnal temperature ranges (C)
	//	real snowforcing(1)  1:7  Forcing variables
 	//	1 airtemp (deg C)
 	//	2 precip (m/h)
 	//	3 windspeed (m/s)
 	//	4 dewpt (deg C)
 	//	5 diurnal trange (deg C)
 	//	6 incoming solar (kJ/m2/h)
 	//	7 net radn (kJ/m2/h)

 	//  real snowsurfacetemp(*)    array of model snow surface temperature extending back 1 day for surface temperature average
 	//	real snowaveragetemp(*)    array of model snow average temperature extending back 1 day for average
 	//	real timestep     Number of hours that model should advance
 	//   nstepday - the number of time steps in a day.   The dimension of daily arrays has to be at least this so that
 	//    daily averages can be tracked
 	//  Variables and arrays to pass to snowlsub
 	const int niv = 7;
 	ArrayXXd inpt(niv,1);
 	ArrayXd outv(23);
 	ArrayXd sitev(8);
 	Array<int,Dynamic,1> iflag(5);
 	double elevsb, xlat, xlon, stdlon, elevtg, rlapse, trange, dewptg, pa, pg, edg, eda, hour1, dt;
 	double templocal, shift, t, tmid, es, cump, cume, cummr;
 	int iyear, month, iday, ihr, imm,  isec, nstepsint, jj;
#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> snowueb(" << ncalls << ")" << std::endl;
    }
	caller = "snowueb";
#endif
 	//     naming inputs that are used
	elevsb = snowsitev(1);
	xlat   = snowsitev(2);
	xlon   = snowsitev(3);
	stdlon = snowsitev(4);
	elevtg = snowsitev(5);
	rlapse = snowsitev(6);
	trange = snowforcing(4);	// diurnal temperature range
	dewptg = snowforcing(3); 	// dew point temperature at temperature gage
 	//   Atmospheric pressure
	pa = 101.3*pow((293.0 - 0.0065*elevsb)/(293.0), 5.256);   // Equation 4.4.12 - Handbook of Hydrology
	pa = pa * 1000.0;  //  convert to pascal
 	//	pa = 1012.4 - 11.34*(elevsb/100) + 0.00745*(elevsb/100)**2.4
 	//     pa = pa * 100.  !  convert to pascal
	pg = 101.3*pow((293.0 - 0.0065*elevtg)/(293.0), 5.256);   // Equation 4.4.12 - Handbook of Hydrology
	pg = pg*1000.0;  //  convert to pascal

 	//   Assume constant mixing ratio for adjustment of dewpoint by elevation
	edg = svp(dewptg);  // Vapour pressure at gage (kPa)
	eda = edg *pa/pg;  //  Vapor pressure adjusted to mean elevation of subbasin (kPa)
	iyear   = snowcontrol(4)/10000;
	month   = snowcontrol(4)/100 - iyear*100;
	iday    = snowcontrol(4) - iyear*10000 - month*100;
	ihr     = snowcontrol(5)/10000;
	imm     = snowcontrol(5)/100 - ihr*100;
	isec    = snowcontrol(5) - ihr*10000 - imm*100;
	hour1   = (double)(ihr) + (double)(imm)/60.0 + (double)(isec)/3600.0;

 	//    Site variables
	sitev(0) = snowsitev(0);  // Forest cover fraction
	sitev(1) = snowparam(22); // Drift factor
	sitev(2) = pa;               // Atmospheric pressure
	sitev(3) = snowparam(23);    // Ground heat flux
	sitev(4) = snowparam(24);    // Albedo Extinction parameter
	sitev(5) = snowsitev(7);     // slope
	sitev(6) = snowsitev(8);     // azimuth
	sitev(7) = xlat;   // latitude

	//    set control flags
	iflag(1-1) = snowcontrol(0);   // radiation is shortwave in (5) and longwave in (6)
	iflag(2-1) = snowcontrol(1);
	iflag(3-1) = snowcontrol(2);
	iflag(4-1) = 1;      // how albedo calculations are done - 1 means albedo is calculated
	iflag(5-1) = snowcontrol(6);      //model option for surface temperature approximation
	nstepsint  = snowcontrol(3);   //  Internal time step
	dt = timestep/(double)nstepsint;
	templocal = snowforcing(0) - (elevsb - elevtg)*rlapse;       //temperature adjustment by lapse rate

 	//      initialize water fluxes
	surfacewaterinput = 0.0;
	snowevaporation   = 0.0;
	areafractionsnow  = 0.0;
 	//  Time loop interpolating external forcing to internal time steps assuming diurnal cycle with tmax at 1500 (3pm)
 	//  and tmin at 0300 (3 am)

 	//    shift time for std longitude versus basin longitude
	shift = (xlon - stdlon)/15.0;
	hour1 += shift;
	for (jj = 1; jj <= nstepsint; jj++) {
		if(nstepsint > 1) {  // interpolate
			tmid = hour1 + dt*0.5;  // time at mid point of interval used for diurnal cycle adjustment
			t = templocal + sin((tmid - 9.0)/24.0*6.283185)*trange*0.5;   // diurnal cycle adjusted temperature
		} else {
			t = templocal;
		}
		es = svp(t);
		//    Mapping inputs onto input array
 		inpt(0,0) = t;                      // temperature to subwatershed
	    inpt(1,0) = snowforcing(1);         // precipitation
		inpt(2,0) = snowforcing(2);         // wind speed
		if(inpt(2,0) < 0.0) {
			inpt(2,0) = 2.0;                //FAO Standard for no wind speed measurement
		}
		inpt(3,0) = min(eda/es, 1.0);       //convert to the relative humidity
		inpt(4,0) = snowforcing(5);         // short wave radiation
		inpt(5,0) = snowforcing(6);         // net radiation
		inpt(6,0) = trange;
		// only call snow routine when there is snow - either on ground or falling.
		if(snowstatev1(2-1) > 0.0 || (t <= snowparam(0) && inpt(2-1,0) > 0.0))  {
            outv(22-1) = 0;
			snowLSub(iyear, month, iday, hour1, dt, 1, inpt, sitev, snowstatev1, snowparam,
                iflag, dtbar, nstepday,	cump, cume, cummr, outv, snowsurfacetemp1, snowaveragetemp1, ndepletionpoints,
                dfc, modelelement,jj);  // 4/2/05  DGT added model element
			surfacewaterinput += outv(22-1)/timestep;           // to get answers in m/hr
			snowevaporation   += outv(23-1)/timestep;
			areafractionsnow  += outv(16-1)/(double)nstepsint;  // averaging the snow covered area over the interval
		} else {
			surfacewaterinput += inpt(2-1,0)*dt/timestep;
			updatetime(iyear, month, iday, hour1, dt);
		}
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving snowueb(" << ncalls << ")" << "\n" << endl;
    }
    ncalls++;
#endif

	return 0;
}

