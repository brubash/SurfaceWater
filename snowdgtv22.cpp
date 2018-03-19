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

//  File  snowdgtv22.f
//  Utah Energy Balance Snow Accumulation and Melt Model.

//  Authors
//  David G. Tarboton (dtarb@cc.usu.edu), Utah Water Research Laboratory, Utah State University
//  Charlie Luce cluce (cluce@fs.fed.us), USDA Forest Service, Rocky Mountain Research Station, Boise, ID
//  Graduate Students who have worked on this code
//    Jinsheng You, PhD, 2004 (youjins@hotmail.com)
//    Tanveer Chowdhury, MS, 1993
//    Tom Jackson, PhD, 1994

//  Significant changes
//  5/4/04 Accommodate glacier surface runoff (version 2.1)
//  7/26/01 Jinsheng You dissertation work (version 2.0)
//    - New surface temperature parameterizations comprising
//       1) the Simple Gradient, almost the same as Original UEB,
//       2) Simple Gradient approach with shallow snow correction.
//       3) The classical force-restore approach.
//       4) Modified force-restore approach.
//    - New refreezing scheme


//  This program was written and prepared under agreement with and funding
//  by the U.S. Government , and therefore is in the public domain and not
//  subject to copyright.

int snowueb2(const double dt, const int nt, const ArrayXXd &input, const ArrayXd &sitev, ArrayXd &statev,
	ArrayXd &tsprevday, ArrayXd &taveprevday, const int nstepday, const ArrayXd &param, Array<int,Dynamic,1> &iflag,
	double &cump, double &cume, double &cummr, ArrayXd &outv, const ArrayXd &mtime, const int modelelement, const int jj)
{

	double mr, ub, ub_old, fc, w, refDepth, totalRefDepth, df, aep;
	double cd, rrhoi, rrho, rid, ta, p, rh, qsi, qnetob, coszen, ws, qli, pRain, ps;
	double a, qh, qe, e, qm, q = 0.0, fm, tave, tsurf, qnet, smelt, rkn;
	int i, ii, pflag, iradfl;

	// Definitions
	//  dt  Time step in hours
	//  nt number of time steps
	// input  -- input forcing
	//	  input(1,*) air temperature (C)
	//	  input(2,*) precipitation (m/hr)
	//	  input(3,*) wind speed  (m/s)
	//	  input(4,*) relative humidity (on a 0-1 scale)
	//	  input(5,*) incoming short wave  (kJ/m^2/h)
	//	  input(6,*) net radiation  (kJ/m^2/h)
	//	  input(7,*) Cosine of Zenith angle
	// SITEV -- site variables
	//        site variables (1-5)
	//        sitev(1)  forest cover fraction
	//        sitev(2)  drift factor  (No detailed information give 1)
	//        sitev(3)  air pressure (Pa)
	//        sitev(4) ground heat flux  Kj/m^2/hr (3.6 = 1 W/m^2)
	//        sitev(5) albedo extinction parameter (m)
	// STATEV
	//        statev(1)  Snow Energy Content  (KJ/m^2)
	//        statev(2)  Snow Water Equivalent (m) relative to T = 0 C solid phase
	//        statev(3)  Dimensionless age of snow surface (or albedo - depending on flag 4)
	//        statev(4)  Refreezing depth (m) used as refdepth
	//        statev(5)  Refreezing depth (m) used as totalrefdepth
	//  totalrefdepth is a misnomer.  These are the same quantity - the refreezing depth.  They are repeated because
	//  when refdepth exceeds the diurnal penetration depth it is set to 0 to tell the code to use the regular
	//  surface temperature functions while totalrefdepth is only reset to 0 when there is actual melt generated
	//  or energy content becomes negative indicating freezing of all liquid phase.  This ensures that the regular
	//  surface temperature functions persist when there is liquid water present but the refreezing front has penetrated
	//  to depth greater than diurnal temperature fluctuations.
	//        TsPrevday(1:nstepday)   Surface temperature over the last 24 hours
	// 	   TavePrevday(1:nstepday)   Depth average temperature over the last 24 hours

	// PARAM  --  snowmelt model parameters (see below)
	// iflag  -- flags
	//	   iflag(1) 0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7
	//        iflag(2)        no 0 (/yes 1) printing
	//        iflag(3)  unit to which to print
	//        iflag(4)  how albedo calculations are done (a value 1 means albedo is calculated, otherwise statev(3) is albedo
	// 	   iflag(5)  model option for surface temperature approximation
	//              1) the Simple Gradient, almost the same as Original UEB,
	//              2) Simple Gradient approach with shallow snow correction.
	//              3) The classical force-restore approach.
	//              4) Modified force-restore approach.
	// cump,cume,cummr  -- cumulative precipitation (with df), cumulative evaporation, cumulative melt over time step in m
	// outv  -- output variables
	//       outv(1)=prain   rain  m/hr
	//       outv(2)=ps     snow  m/hr
	//       outv(3)=a     albedo
	//       outv(4)=qh    sensible heat (kJ/m2/hr)
	//       outv(5)=qe    latent heat (kJ/m2/hr)
	//       outv(6)=e     sublimation m/hr
	//       outv(7)=mr    melt outflow m/hr
	//       outv(8)=qm    heat advected with melt
	//       outv(9)=q     Energy rate of change (kJ/m2/hr)
	//       outv(10)=fm   Mass rate of change (m/hr)
	//       outv(11)=tave  Average temperature (C)
	//       outv(12)=tsurf  Surface temperature C
	//       outv(13)=qnet  Net Radiation (kJ/m2/hr)
	//	 outv(14)=smelt   Surface melt  m/hr

	// mtime   4 variable array with time information
	//    mtime(1)  year
	//    mtime(2)  month
	//    mtime(3)  day
	//    mtime(4)  hour (0-24 as a real number)
	//    mtime(5)  model element number

	// Note: in this subroutine, the outv is an array which passes value to this subroutine and back to the snow
	// drive program. The Outv(9) and outv(10) pass the daily average snowpack temperature and daily
	// snow surface temperature to this subroutine but pass the Qin total and combined mass fluxes back.

	int iTsMethod;       // Add model time initialization 09/19/2000

	//  Constant data set
	const double to    = 0.0;			//  Temperature of freezing (0 C)
	const double tk    = 273.15;		//  Temperature to convert C to K (273.15)
	const double hf    = 333.5;			//  Heat of fusion (333.5 KJ/kg)
	const double cs    = 2.09;			//  Ice heat capacity (2.09 KJ/kg/C)
	const double k     = 0.4;			//  Von Karmans constant (0.4)
	const double hff   = 3600.0;		//  Factor to convert /s into /hr (3600)
	const double rhoi  = 917.0;			//  Density of Ice (917 kg/m^3)
	const double rhow  = 1000.0;		//  Density of Water (1000 kg/m^3)

	//  Parameters
	const double tr    = param(1-1);	//  Temperature above which all is rain (3 C)
	const double ts    = param(2-1);	//  Temperature below which all is snow (-1 C)
	const double cg    = param(4-1);	//  Ground heat capacity (nominally 2.09 KJ/kg/C)
	const double z     = param(5-1);	//  Nominal meas. height for air temp. and humidity (2m)
	const double zo    = param(6-1);	//  Surface aerodynamic roughness (m)
	const double rho   = param(7-1);	//  Snow Density (Nominally 450 kg/m^3)
	const double rhog  = param(8-1);	//  Soil Density (nominally 1700 kg/m^3)
	const double lc    = param(9-1);	//  Liquid holding capacity of snow (0.05)
	const double de    = param(11-1);	//  Thermally active depth of soil (0.1 m)
	const double abg   = param(12-1);	//  Bare ground albedo  (0.25)
	const double avo   = param(13-1);	//  Visual new snow albedo (0.95)
	const double anir0 = param(14-1);	//  NIR new snow albedo (0.65)
	const double dNewS = param(21-1);	//  The threshold depth of for new snow (0.001 m)
	//const double gsurf = param(22-1);   //  The fraction of surface melt that runs off (e.g. from a glacier)
	// Bert Rubash: gsurf is passed to predicorr but not used there.

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> snowueb2(" << ncalls << ")" << std::endl;
    }
	caller = "snowueb2";
#endif
	//   debugging
	//      if (mtime(1) .eq. 1948 && mtime(2) .eq. 11 &&
	//     + mtime(3) .eq. 8 && mtime(4)  >  12 &&
	//     + mtime(5) .eq. 169) {
	//	   mtime(1)=mtime(1)
	//	endif
	//      cerr << (input(j,1),j=1,7)
	//   State variables - These serve as initial conditions
	ub = statev(0);    // Snow Energy Content  (KJ/m^2)
	w  = statev(1);    // Snow Water Equivalent (m) relative to T = 0 C solid phase
	if (ub <= 0.0) {
		refDepth = 0.0;
		totalRefDepth = 0.0;
	} else {
		refDepth = statev(4-1);
		totalRefDepth = statev(5-1);
	}

	//	Save old Value 07/23/01
	ub_old = ub;

	//   Site variables
	fc = sitev(1-1);     //  Forest cover fraction (0-1)
	df = sitev(2-1);     //  Drift factor
	aep = sitev(5-1);   //  Albedo extinction parameter to smooth
	//      transition of albedo when snow is shallow. Depends on Veg. height (m)

	//   control flags
	iradfl = iflag(0);
	pflag  = iflag(1);
	//iflag(4) albedo caculation
	iTsMethod = iflag(5-1);	// the method to approximate the surface temperature
							// 1 normal old snow melt model
							// 2 revised direct gradient method (New method for ke) and qe
							// 3 force restore approach
							// 4 modified force restore approach

	//   Calculate constants that need only be calculated once.
	cd = k*k*hff/(log(z/zo)*log(z/zo))*(1.0 - 0.8*fc);   // factor in turbulent fluxes
	//    The FC correction is applied here to hit sensible and latent heat fluxes
	//    and evaporation rate in one place and maintain consistency in the surface
	//     temperature and stability calculations. It is consistent with wind speed
	//     reduction by a factor of 1-0.8 FC which is the most reasonable physical
	//     justification for this approach.
	//     I recognize that this is not a very good way to parameterize vegetation.
	//     This will be an area of future improvements (I hope).
	//    FC is also used below to adjust radiation inputs.
	rrhoi = rhoi/rhow;
	rrho  = rho/rhow;
	rid   = 1.0/rrho - 1.0/rrhoi;

	//   Loop for each time step
	for (i = 1; i <= nt; i++) {   // DGT Time looping disabled as temperature averaging not handled here - must be handled outside
		//   Input variables
        ta = input(1-1,i-1);    // Air temperature input (Degrees C)
        p  = input(2-1,i-1);     // Precipitation rate input (m/hr)
        ws = input(3-1,i-1);    // Wind Speed (m/s)
        rh = input(4-1,i-1);    // Relative humidity (fraction 0-1)
        // DGT 5/27/12 initialize variables to avoid run time check problems
        qsi = 0.0;
        qli = 0.0;
        qnetob = 0.0;
        if (iradfl == 0) {  	// input is incoming short and longwave
			qsi = input(5-1,i-1);  	// Incoming shortwave radiation (KJ/m^2/hr)
			qli = input(6-1,i-1);  	// Incoming longwave radiation (KJ/m^2/hr)
        } else {
          qnetob = input(5-1,i-1); // Net allwave radiation (KJ/m^2/hr)
        }
        coszen = input(7-1,i-1);   // Cos(angle between direct sunlight and surface
		//                             normal).
		//         Representative value over time step used in albedo calculation.
		//         We neglect the difference between direct and diffuse albedo.
		// Daily average temperatures handled internally so that multiple time steps will work
		ts_save::Ts_Ave   = daily_ave(tsprevday, nstepday, -100.0) + tk; 	// (C)  Surface temperature average over last 24 hours
		ts_save::Tave_ave = daily_ave(taveprevday, nstepday, -100.0) + tk; 	// (C)  Depth averaged temperature average over last 24 hours
		ts_save::Ts_old   = tsprevday(nstepday-1) + tk; 					// (C) Surface temperature from previous time step
		ts_save::Tave_old = taveprevday(nstepday-1) + tk; 					// (C) Average temperature from previous time step
		//   If any of these variables are out of range due to any problem set them back to freezing
		if (ts_save::Ts_old < 0.0) {
			cerr << "Invalid previous time step surface temperature " << ts_save::Ts_old << " set to 273 K\n";
			ts_save::Ts_old = tk;
		}
		if (ts_save::Tave_old < 0.0) {
			cerr << "Invalid previous time step average temperature " << ts_save::Tave_old << " set to 273 K\n";
			ts_save::Tave_old = tk;
		}
		if (ts_save::Ts_Ave < 0.0) {
			cerr << "Invalid last 24 hr average surface temperature " << ts_save::Ts_Ave << " set to 273 K\n";
			ts_save::Ts_Ave = tk;
		}
		if (ts_save::Tave_ave < 0.0) {
			cerr << "Invalid last 24 hr average temperature " << ts_save::Tave_ave << " set to 273 K\n";
			ts_save::Tave_ave = tk;
		}

		//  Separate rain and snow
        ps = partsnow(p, ta, tr, ts);
        pRain = p - ps;
		//  Increase precipitation as snow by drift multiplication factor
        ps = ps*df;

		//        if (iflag(4).eq.1) {
		//  Calculate albedo
		//          a=albedo(statev(3),coszen,w/rrho,aep,abg,avo,anir0)
		// Use of this albedo throughout time step neglects the
		//  changes due to new snow within a time step.
		//        else
		//          a=statev(3)
		//        endif

		//  if there is any new snow without rain the snow age is change to 0.

		//       if (ps  >  0.0 && prain .le. 0.0) then
		//	       statev(3) =0.0
		//	  endif

		if (iflag(4-1) == 1) {
			//  Calculate albedo
			a = albedo(statev(3-1), coszen, w/rrho, aep, abg, avo, anir0);
			// Use of this albedo throughout time step neglects the
			//  changes due to new snow within a time step.o
		} else {
			a = statev(3-1);
		}

		//   Calculate neutral mass transfer coefficient
		rkn = cd*ws;
		//   Adjust radiation inputs for the effect of forest canopy.
		qsi = qsi*(1.0 - fc);
		qli = qli*(1.0 - fc);
		qnetob = qnetob*(1.0 - fc);
		//     FC corrections are also in the following subroutines.
		//      qfm where outgoing longwave radiation is calculated.
		//      SNOTMP where outgoing longwave radiation is used in the surface
		//      temperature equilibrium approx.
		//  Debugging
		//      if (mtime(1)== 1950 && mtime(2)==4 &&
		//     + mtime(3)==1 && mtime(4) >12 && mtime(5)==169) {
		//	   mtime(1)=mtime(1)
		//	}
		//   Call predictor corrector subroutine to do all the work
		predicorr(dt, ub, w, a, ta, pRain, ps, ws, rh, qsi, qli, iradfl, rkn,  qnetob, rid, param, sitev,
			iTsMethod, mtime, qh, qe, e, mr, qm, q, fm, tsurf, tave, qnet, refDepth, totalRefDepth, smelt);

		//  DGT 4/2/05   Despite all checks in predicor It can (and does) occur that
		//   we still have ub so large that it results in tave greater than 0, which implies that all the
		//   snow is liquid.  In these cases - just force the snow to disappear and add the energy involved to Qm.
		tave = tavg(ub, w, rhow, cs, to, rhog, de, cg, hf);
		//	if (tave < -1000) {
		//	   tave=tave
		//	}
		if (tave  >  0.0) {   //  all is liquid so snow must disappear
			mr = mr + w/dt;
			qm = qm + w/dt*rhow*hf;
			q = q - w/dt*rhow*hf;
			ub = ub - w/dt*rhow*hf;
			w = 0.0;
		}
		//  To guard against unreasonable UB when there is no snow do not allow bulk temperature to go above 10 C
		if (tave  >  10.0) {
			ub = rhog*de*cg*10.0;
		}
		// surface melt change
		//   Update snow surface age based on snowfall in time step
		//       if (iflag(4).eq.1) call agesn(statev(3),dt,ps,tsurf,tk)
		//   Update snow surface age based on snowfall in time step
		if (iflag(4-1) == 1)
			agesn(statev(3-1), dt, ps, tsurf, tk, dNewS);
		//    accumulate for mass balance
		cump  = cump + (ps + pRain)*dt;
		cume  = cume + e*dt;
		cummr = cummr + mr*dt;
		tave  = tavg(ub, w, rhow, cs, to, rhog, de, cg, hf);   //  this call
		//   necessary to keep track of average internal temperature used in some surface energy algorithms.
		//	if (tave < -1000) {
		//	   tave=tave
		//	}

		// update the total depth of the refreezing depth in the snow pack according the
		// the refreezing depth at time step and the positive energy input. 07/22/01
		//  DGT's revised logic  1/13/05
		if (lc > 0.0) {
			if (refDepth >  0.0) {
				totalRefDepth = refDepth;  // if refreezing depth increases totalrefdepth keeps track
			} else { // here refdepth has gone to 0 somehow
				if (mr > 0.0 || (ub > ub_old && ub  > 0.0))	// If there is melt or an increase in energy refdepth is reset
					totalRefDepth = 0.0;
			}
		} else if (mr > 0.0 || (ub > ub_old && ub  > 0.0)) {
			//   Here lc=0.  If there is melt or an increase in energy refdepth is reset
			//   This is likely redundant because if lc=0 then there is no meltwater to refreeze
			totalRefDepth = 0.0;
		}
		//  Jinsheng's original logic
		//      if (lc > 0.0) then
		//       if ((refDepth-refDepth_old)  >  0.0) then
		//	    totalRefDepth=totalRefDepth +refDepth-refDepth_old
		//	 else
		//	  if (mr > 0.0)  totalRefDepth = 0.0 //totalRefDepth-mr*rhow/rhom
		//	  if ((ub-ub_old)  >  0.0 && ub. gt.0.0) then
		//	    TotalRefDepth = totalRefDepth-(ub-ub_old)/(rhom*hf)
		//	  endif
		//	 endif
		//	elseif (mr > 0.0 .or. (ub > ub_old && ub  > 0.0)) then
		//	 totalRefDepth =0.0
		//	endif
		if (totalRefDepth < 0.0)
			totalRefDepth = 0.0;

		// update tsbackup and tavebackup
		for (ii = 0; ii < nstepday-1; ii++) {
			tsprevday(ii)   = tsprevday(ii+1);
			taveprevday(ii) = taveprevday(ii+1);
		}
		tsprevday(nstepday-1)   = tsurf;
		taveprevday(nstepday-1) = tave;

		if (pflag == 1) {
            snowcontrol3_File << ub    << " " << w     << " " << statev(3-1) << " ";
            snowcontrol3_File << pRain << " " << ps    << " " << a           << " ";
            snowcontrol3_File << qh    << " " << qe    << " " << e           << " ";
            snowcontrol3_File << mr    << " " << qm    << " " << q           << " ";
            snowcontrol3_File << fm    << " " << tave  << " " << tsurf       << " ";
            snowcontrol3_File << cump  << " " << cume  << " " << cummr       << " ";
            snowcontrol3_File << qnet  << " " << smelt << " " << refDepth    << " " << totalRefDepth << '\n';
		// 5/4/04 surface melt smelt
		}
	}

	statev(0) = ub;
	statev(1) = w;
	statev(3) = refDepth;
	statev(4) = totalRefDepth;
	outv(0)   = pRain;
	outv(1)   = ps;
	outv(2)   = a;
	outv(3)   = qh;
	outv(4)   = qe;
	outv(5)   = e;
	outv(6)   = mr;
	outv(7)   = qm;
	outv(8)   = q;
	outv(9)  = fm;
	outv(10)  = tave;
	outv(11)  = tsurf;
	outv(12)  = qnet;
	outv(13)  = smelt;

#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving snowueb2(" << ncalls << ")" << "\n" << endl;
    }
    ncalls++;
#endif

	return 0;
}


// *********************** PREDICORR () **************************
//     Predictor-corrector scheme to update the state variables,
//     U and W for each time step

int predicorr(const double dt, double &ub, double &w, const double a,
	const double ta, const double pRain, const double ps, const double ws, const double rh,
	const double qsi, const double  qli, const double iradfl, const double rkn,
	const double qnetob, const double rid,
	const ArrayXd &param, const ArrayXd &sitev, const double iTsMethod,
	const ArrayXd &mtime,  // pass a modeling time
	//     following variables are output
	double &qh, double &qe, double &e,  double &mr, double &qm, double &q, double &fm,
	double &tsurf, double &tave, double &qnet, double &refDepth, double &totalRefDepth ,
	double &smelt)
{
	// surface melt smelt
	//  Constant data set
	const double hf    = 333.5;			//  Heat of fusion (333.5 KJ/kg)
	const double hneu  = 2834.0;		//  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)f
	const double rhow  = 1000.0;		//  Density of Water (1000 kg/m^3)

	double mr1, qh1, qe1, e1, qm1, tsurf1, qnet1, q1, smelt1, w1;
	double ub1, fm1, w2, ub2, ae, e2, rlf, r;
	int imax, niter;

	const double wtol = 0.025;
	const double utol = 2000.0;

	qfm(ub, w, a, ta, pRain, ps, ws, rh, qsi, qli, rkn, iradfl, qnetob,
		rid, param, sitev, iTsMethod, mtime,
		fm, q, qm, mr, qe, e, tsurf, tave, qh, qnet, dt, refDepth,
		totalRefDepth, smelt);

	// surface melt smelt
	//      PREDICTOR
	w1 = w + dt*fm;
	if (w1 < 0.0) {
		w1 = 0.0;
		prehelp(w1, w, dt, fm, 0.0, 1.0, ps, pRain, e, rhow, hf, q, qm, mr, qe, hneu);
	}
	ub1 = ub + dt*q;
	q1 = q;
	fm1 = fm;
	//   save values so that they can be averaged for output
	qh1 = qh;
	qe1 = qe;
	e1 = e;
	mr1 = mr;
	smelt1 = smelt;  //cdgt 5/4/04 surface melt smelt

	qm1 = qm;
	tsurf1 = tsurf;
	qnet1 = qnet;

	qfm(ub1, w1, a, ta, pRain, ps, ws, rh, qsi, qli, rkn, iradfl, qnetob,
		rid, param, sitev, iTsMethod, mtime,
		fm, q, qm, mr, qe, e, tsurf, tave,  qh, qnet, dt, refDepth,
		totalRefDepth, smelt);
	// surface melt smelt

	//      CORRECTOR
	w2 = w + dt/2.0*(fm1 + fm);
	if (w2 < 0.0) {
		w2 = 0.0;
		prehelp(w2, w, dt, fm, fm1, 2., ps, pRain, e, rhow, hf, q, qm, mr, qe, hneu);
	}
	ub2 = ub + dt/2.0*(q1 + q);
	//   iterate to convergence to enhance stability
	niter = -3;
	imax = 1;
L1: if ((fabs(w2-w1) > wtol || fabs(ub2-ub1) > utol) && niter < imax) {
		w1 = w2;
		ub1 = ub2;
		qfm(ub1, w1, a, ta, pRain, ps, ws, rh, qsi, qli, rkn, iradfl,
			qnetob, rid, param, sitev, iTsMethod, mtime,
			fm, q, qm, mr, qe, e, tsurf, tave, qh, qnet, dt, refDepth,
			totalRefDepth, smelt);
		// surface melt smelt
		//    corrector again
		w2 = w + dt/2.0*(fm1 + fm);
		if (w2 < 0.0) {
			w2 = 0.0;
			prehelp(w2, w, dt, fm, fm1, 2., ps, pRain, e, rhow, hf, q, qm, mr, qe, hneu);
		}
		ub2 = ub + dt/2.0*(q1 + q);
		niter++;
		if (niter >= 1) {  // had * steps to converge now hit it.
			//  What follows is a fix to numerical instability that results from
			//  nonlinearity when the snowpack is shallow and melting rapidly.  If
			//  convergence does not occur when the snowpack is not melting (a very
			//  rare thing) I just accept the predictor corrector solution.

			//  Modified by DGT 7/27/05 to add complete meltout condition
			//  The quantities that this changes are w2, ub2, mr and qm

			//  The fix first checks whether the added energy is sufficient to melt the snow
			//  completely.  If this is the case then the snow disappears.
			//  In other cases we assume that the liquid fraction of the snow remains constant.
			//  This to some extent bypasses the melt outflow estimates.
			//  ae is added energy during the time step.
			ae = (q1 + q + qm1 + qm)*0.5*dt;
			//   This fix is only physically sensible under melting conditions
			//   and when ae is positive and there is snow
			if (ub  >  0.0 && ae  >  0.0 && w  >  0.0) {
				e2 = (e + e1)*0.5;   //  This is the average sublimation
				//   Check liquid fraction with added energy.  If this is greater than 1 then all snow melts
				//   Otherwise implement a solution assuming that the liquid fraction remains constant
				rlf = (ub + ae)/(rhow*w*hf);
				if (rlf >= 1.0) {
					mr = w/dt + (ps + pRain - e2);   // Here snow disappears
					if (mr < 0.0) {
						mr = 0.0;   //  Force this to not go negative
						//      This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
						//      There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
						//      This is ignored because once the snow is gone energy computations are not pertinent.
					}
					qm = mr*rhow*hf;
					w2 = 0.0;
					ub2 = ub + ae - qm*dt;
				} else {
					//   Determine the w/ub ratio at the beginning of the time step.
					//   Physically liquid fraction = ub/(rhow*w*hf) and since rhow and hf are constants
					//   keeping r=w/ub constant is equivalent to keeping liquid fraction constant.
					//   If liquid fraction is constant this is constant.
					r = w/ub;
					//   Solve for ub2 and w2 that satisfy the three equations
					//            r=w2/ub2
					//            ub2=ub+ae-rhow*hf*mr*dt     Energy balance the last term being the energy advected by melt
					//            w2=w+(ps+prain-e2-mr)*dt    Mass balance
					//   The unknowns in the above are ub2, w2 and m and the equations are linear
					//   once the first eqn is multiplied by ub2
					//   The solution is
					ub2 = (rhow*hf*(w + (ps + pRain - e2)*dt) - ae - ub)/(rhow*hf*r - 1.0);
					w2 = r*ub2;
					if (w2 < 0.0) {  // Avoid negative W
						w2 = 0.0;
					}
					mr = (w - w2)/dt - e2 + ps + pRain;
					if (mr < 0.0) {   // Avoid negative mr
						mr=0.0;
						w2 = w + (ps + pRain - e2)/dt;
						if (w2 < 0) {
							w2 = 0.0;
							//      This can only occur if e2 is large compared to other terms.  Setting w2=0 implicitly reduces e2.
							//      There is a resulting energy discrepancy due to a limitation on sublimation and latent heat flux
							//      This is ignored because once the snow is gone energy computations are not pertinent.
						}
					}
					qm = mr*rhow*hf;
					ub2 = ub + ae - qm*dt;   // redundant most of the time but recalc due to exceptions
				}
				//    Check that nothing went wrong
				if (mr < 0.0) {
					cerr << "Error - negative melt rate in snow\n";
				}
				if (w2 < 0.0) {
					cerr << "Error - negative w2 in snow\n";
				}
				q = ae/dt - qm;
				//   Now set first pass values equal to final values to fake the averages below
				qm1 = qm;
				mr1 = mr;
				q1 = q;
			}
		}
		goto L1;
	}
	w = w2;
	ub = ub2;
	//  average values from two time steps for output.  This is done for mr
	//  and e to ensure mass balance and the others for better physical
	//  comparisons
	qh    = (qh    + qh1   )*0.5;
	qe    = (qe    + qe1   )*0.5;
	e     = (e     + e1    )*0.5;
	mr    = (mr    + mr1   )*0.5;
	qm    = (qm    + qm1   )*0.5;
	tsurf = (tsurf + tsurf1)*0.5;
	qnet  = (qnet  + qnet1 )*0.5;
	q     = (q     + q1    )*0.5;
	smelt = (smelt + smelt1)*0.5/(hf*rhow);  // surface melt smelt
	// convert from energy KJ/m^2/hr to water depth m/hr of melt.
	return 0;
}


// ************************ QFM () ********************************
//     Calculates Energy and Mass Flux at any instant

int qfm(const double ub, const double w, const double a, const double ta, const double pRain,
	const double ps, const double ws, const double rh, const double qsi, const double qli,
	const double rkn, const double iradfl, const double qnetob,
	const double rid, const ArrayXd &param, const ArrayXd &sitev, const int iTsMethod, const ArrayXd &mtime,
	double &fm, double &q, double &qm, double &mr, double &qe,
	double &e, double &tsurf, double &tave, double &qh, double &qnet,
	const double dt, double &refDepth, const double totalRefDepth, double &smelt)
{
	double LanS, fc, pr, qg, qsn, tak, tavek, densa, dens, rhom, fKappaS, rs, ds, ea, qp;
	double qc1, qc2, var_a, var_b, x1, qle, qlnet;

	// Constant data set
	const double to    = 0.0;			//  Temperature of freezing (0 C)
	const double tk    = 273.15;		//  Temperature to convert C to K (273.15)
	const double sbc   = 2.0747e-7;		//  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
	const double hf    = 333.5;			//  Heat of fusion (333.5 KJ/kg)
	const double hneu  = 2834.0;		//  Heat of Vaporization (Ice to Vapor, 2834 KJ/kg)f
	const double cw    = 4.18;			//  Water Heat Capacity (4.18 KJ/kg/C)
	const double cs    = 2.09;			//  Ice heat capacity (2.09 KJ/kg/C)
	const double cp    = 1.005;			//  Air Heat Capacity (1.005 KJ/kg/K)
	const double ra    = 287.0;			//  Ideal Gas constant for dry air (287 J/kg/K)
	const double rhow  = 1000.0;		//  Density of Water (1000 kg/m^3)
	const double g     = 9.81;			//  Gravitational acceleration (9.81 m/s^2)
	const double w1day = 0.261799;		//  Daily frequency (2pi/24 hr 0.261799 radians/hr)

	//  Parameters
	const double es    = param(2);	//  emmissivity of snow (nominally 0.99)
	const double cg    = param(3);	//  Ground heat capacity (nominally 2.09 KJ/kg/C)
	const double z     = param(4);	//  Nominal meas. height for air temp. and humidity (2m)
	const double rho   = param(6);	//  Snow Density (Nominally 450 kg/m^3)
	const double rhog  = param(7);	//  Soil Density (nominally 1700 kg/m^3)
	const double lc    = param(8);	//  Liquid holding capacity of snow (0.05)
	const double ks    = param(9);	//  Snow Saturated hydraulic conductivity (20 !160 m/hr)
	const double de    = param(10);	//  Thermally active depth of soil (0.1 m)
    LanS               = param(14);	//  the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	const double rd1   = param(17);	//  Amplitude correction coefficient of heat conduction (1)
	const double fstab = param(18);	//  Stability correction control parameter 0 = no corrections, 1 = full corrections
	const double gsurf = param(21);	//  The fraction of surface melt that runs off (e.g. from a glacier)

	if (isnan(q)) {
        cerr << "q in qfm() is not a number" << endl;
		exit(0);
	}
	//   Site variables
	fc = sitev(1-1);     //  Forest cover fraction (0-1)
	if (fc > 0) {
	   fc = fc;
	}
	//      df=sitev(2)    !  Drift factor
	pr = sitev(3-1);     //  Atmospheric Pressure (Pa)
	qg = sitev(4-1);     //  Ground heat flux (KJ/m^2/hr)  This is more logically an
	//      input variable, but is put here because it is never known at each
	//      time step.  Usually it will be assigned a value 0.

	qsn = qsi*(1.0 - a);
	//     To ensure all temperatures in kelvin
	tak = ta + tk;
	tavek  = tave + tk;
	densa  = pr/(ra*tak);     // Density of Air in kg/m3 if PR in Pa
	dens = rho;
	rhom = lc*rho;

	// calculate the ks=k/ze
	fKappaS = LanS/(rho*cs);
	ds = sqrt(2.0*fKappaS/w1day);
	rs = LanS/(ds*rd1*rho*cs);


	// save the old values
	ts_save::Tave_old = tave;

	ea = svpw(ta)*rh;   // The saturation vapour pressure over water is used
	//       because most instruments report relative humidity relative to water.
	qp = qpf(pRain, ta, to, ps, rhow, hf, cw, cs);
	tave = tavg(ub, w, rhow, cs, to, rhog, de, cg, hf);

	// as described as below, the original UEB does not consider the refreezing. in this
	// change, the model will consider the refreezing effect starting from the top to the bottom.
	// Here is the predictor.

	// the method is: refreezing can be modeled if the ub is greater than 0.0 without any precipiation,
	// meanwhile the total refreezing depth in the snowpack is less than the Refreezing depth times the daily damping of wet snow.

	if (ub  > 0 && ps <= 0.0 && pRain <= 0.0 && totalRefDepth <= rd1*ds && w > 0.0) {
		qc1 = QcEst(tk, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt);
		qc2 = QcEst(tk-0.01, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt);


		Grad(qc1, qc2, 0.0, -0.01,  var_a, var_b);
		x1 = refDep(LanS, var_a, var_b, hf, rhom, dt, refDepth); //refreezing depth
		refDepth = x1;


	} else {
        refDepth = 0.0;
	}

	//       cerr << "Tave in qfm",tave
    tsurf = srftmp(qsi, a, qli, qp, ea, ta, tave, tk, pr, ra, cp, rho,
		rkn, hneu, es, sbc, cs, rs, w, qnetob, iradfl, ws, z, g, fc,
		fstab, mtime, param, iTsMethod, dt, ub, refDepth, smelt);
	// surface melt smelt
	//       cerr << tsurf
	qle = (1.0 - fc)*es*sbc*pow(tsurf + tk, 4);
	qlnet = qli - qle;

	turbflux(pr, ra, ta, tk, tsurf, z, g, cp, rkn, ws, ea, rhow, hneu, qh, qe, e, fstab); // add a fraction to reduce the evaporation after snow gone
	mr = fmelt(ub, rhow, w, hf, lc, rid, ks, pRain);

	//        MR in m/hr
	qm = mr*rhow*(hf + (tave - to)*cw);        //  Include advection of
	//         meltwater/rain that is at tave so that the model does not crash when there is no snow and it rains.
	//        QM in kj/m2/hr
	//  Add surface melt that runs off to melt runoff and melt energy so that energy and mass balances are still correct

	mr += smelt/(hf*rhow)*gsurf;
	qm += smelt*gsurf;
	if (iradfl == 0) {
		qnet = qsi*(1.0 - a) + qlnet;
	} else {
		qnet = qnetob;
	}

	q = qnet + qp + qg + qh + qe - qm;
	fm = pRain + ps - mr - e;

	if (isnan(q)) {
        cerr << "q in qfm() is not a number" << endl;
		exit(0);
	}

	return 0;
}

// **************************** TAVG () ***************************
//     Calculates the average temperature of snow and interacting soil layer

double tavg(const double ub, const double w, const double rhow, const double cs,
	const double to, const double rhog, const double de, const double cg, const double hf)
{
	double ts, snhc, shc, chc, al;

	snhc = rhow*w*cs;
	shc  = rhog*de*cg;
	chc  = snhc + shc;

	//      snhc = Snow heat capacity
	//      shc  = Soil heat capacity
	//      chc  = Combined heat capacity

	if (ub <= 0.0) {
		ts = ub/chc;
	} else {
		al = ub/(rhow*hf);
		if (w > al) {
			ts = to;
		} else {
			ts = (ub - w*rhow*hf)/chc;
		}
	}

	return ts;
}


// **************************** SRFTMP () ***************************
//     Computes the surface temperature of snow
double srftmp(const double qsi, const double a, const double qli, const double qpin,
		const double ea, const double ta, const double tave, const double tk, const double pr,
		const double ra, const double cp,const double rho, const double rkn, const double hneu,
		const double es, const double sbc, const double cs, const double rs, const double w,
		const double qnetob, const double iradfl, const double ws, const double z, const double g,
		const double fc, const double fstab, const ArrayXd &mtime, const ArrayXd &param,
        const int iTsMethod, const double dt, const double ub, const double refDepth, double &smelt)
{
	//   This version written on 4/23/97 by Charlie Luce solves directly the
	//   energy balance equation using Newtons method - avoiding the linearizations
	//   used previously.  The derivative is evaluated numerically over the range
	//   ts to fff*ts  where fff = 0.999
	const double tol = 0.05;
	const int nitermax = 10;
	int niter, ibtowrite, iter;
	double fff, qsn, tak, tavek, densa, dens, qp, ts, er, tslast, f1, f2, tlb, tub, flb, fub;
	double result;

	fff = (273.0 - tol)/273.0;  // Scale the difference used for derivatives by the tolerance
	qsn = qsi*(1.0 - a);
	//     To ensure all temperatures in kelvin
	tak   = ta + tk;
	tavek = tave + tk;
	densa = pr/(ra*tak);     // Density of Air in kg/m3 if PR in Pa
	dens  = rho;

	qp = qpin;            // store input variable locally without changing global value
	//    debugging
	//      if (mtime(1) .eq. 1948 && mtime(2) .eq. 11
	//     + && mtime(3) .eq. 18 && mtime(4)  >  12 &&
	//     + mtime(4) < 13 && mtime(5) .eq. 39) {
	//		mtime(1)=mtime(1)
	//	endif


	if (w <= 0.&& qp > 0.0)
		qp = 0.0;
	//
	//      ignore the effect of precip advected
	//      energy on the calculation of surface temperature when there is no snow.
	//      Without this ridiculously high temperatures can result as the model
	//      tries to balance outgoing radiation with precip advected energy.

	ts = tak;                             // first approximation
	er = tol*2.0;    // so that it does not end on first loop
	niter = 0;
L1: if (er > tol && niter < nitermax) {
		tslast = ts;
        f1 = surfeb(ts, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
		f2 = surfeb(fff*ts, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
		ts = ts - ((1.0 - fff)*ts*f1)/(f1 - f2);
		if (ts < tak - 50)
			goto L11; //  If it looks like it is getting unstable go straight to bisection
		er = fabs(ts - tslast);
		niter++;
		goto L1;   // loop back and iterate
	}

	if (er <= tol)
		goto L10;  // The solution has converged

	//   If still not converged use bisection method
L11: tlb = tak - 20.0;        // First guess at a reasonable range
	tub = tak + 10.0;
	flb = surfeb(tlb, rkn, ws, tak, g, qp, densa, cp, hneu,
		pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
		mtime, param, iTsMethod, w, dt, ub, refDepth);
	fub = surfeb(tub, rkn, ws, tak, g, qp, densa, cp, hneu,
		pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
		mtime, param, iTsMethod, w, dt, ub, refDepth);
	ibtowrite = 0;
	if (flb*fub  >  0.) {   // these are of the same sign so the range needs to be enlarged
		tlb = tak - 150.0;    // an almost ridiculously large range - solution should be in this if it exists
		tub = tak + 100.0;
		flb = surfeb(tlb, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
		fub = surfeb(tub, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
		ibtowrite = 1;
		if (flb*fub  >  0.0) {   // these are of the same sign so no bisection solution
			cerr << "Bisection surface temperature solution failed with large range\n";
			cerr << "Date: " << " " << mtime(1-1) << " " << mtime(2-1) << " " << mtime(3-1) << '\n';
			cerr << "Time: " << " " << mtime(4-1) << '\n';
			cerr << "Model element: " << " " << mtime(5-1) << '\n';
			cerr << "A surface temperature of 273 K assumed\n";
	        ts = tk;
			goto L10;
	    } else {
			cerr << "Bisection surface temperature solution with large range\n";
			cerr << "Date: " << " " << mtime(1-1) << " " << mtime(2-1) << " " << mtime(3-1) << '\n';
			cerr << "Time: " << " " << mtime(4-1) << '\n';
			cerr << "Model element: " << mtime(5-1) << '\n';
			cerr << "This is not a critical problem unless it happens frequently\n";
			cerr << "and solution below appears incorrect\n";
	    }
	} else {
		//		cerr <<
		//     +    "Bisection surface temperature solution"
		//		cerr << "Date: ",mtime(1),mtime(2),mtime(3)
		//		cerr << "Time: ",mtime(4)
		//		cerr << "Model element: ",mtime(5)
		//	    cerr << "This is not a critical problem"
	}
	//     Here do the bisection
    niter = log((tub - tlb)/tol)/log(2.0);   // Number of iterations needed for temperature tolerance
	for (iter = 1; iter <= niter; iter++) {
		ts = 0.5*(tub + tlb);
		f1 = surfeb(ts, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
		if (f1 > 0.0) {  // This relies on the assumption (fact) that this is a monotonically decreasing function
			tlb = ts;
		} else {
			tub = ts;
		}
	}
	if (ibtowrite == 1) {
		cerr << "Surface temperature: " << " " << ts << " " << " K\n";
		cerr << "Energy closure: " << " " << f1 << '\n';
		cerr << "Iterations: " << " " << niter << '\n';
	}

L10: ts = ts - tk;
    if (w > 0.0 && ts > 0.0) {
		// surface melt smelt
		result = 0.0;
	    smelt = surfeb(result+tk, rkn, ws, tak, g, qp, densa, cp, hneu,
			pr, ea, tk, dens, cs, rs, tavek, qsn, qli, fc, sbc, qnetob, iradfl,
			mtime, param, iTsMethod, w, dt, ub, refDepth);
			// surface melt smelt is the energy not accommodated by conduction into the snow
			//  so it results in surface melt which then infiltrates and adds energy to the snowpack
			//  through refreezing
			//  This change added to capture this quantity when the model is used for a glacier where
			//  surface melt may run off rather than infiltrate.
			//  For modeling of glaciers in Antarctica by Mike Gooseff
	} else {
          result = ts;
	    smelt = 0.0;
		//  No surface melt this case
	}

    return result;
}


// *********************************************************************
double surfeb(const double ts, const double rkn, const double ws, const double tak, const double g,
	const double qp, const double densa, const double cp, const double hneu, const double pr, const double ea,
	const double tk, const double dens, const double cs, const double rs, const double tavek, const double qsn,
	const double qli, const double fc, const double sbc, const double qnetob, const int iradfl,
	const ArrayXd &mtime, const ArrayXd &param, const int iTsMethod, const double w, const double dt,
	const double ub, const double refDepth)
{
	//      function to evaluate the surface energy balance for use in solving for
	//      surface temperature
	double LanS, LanG, LanE_Ze, LanE_De, LanE_Ze2, LanE_De2, Ze2;
	double zs, rkin, qcs, fKappaS, fKappaG, d1, dlf, ze, result, de2;

	// Constant data set
	const double rhow  = 1000.0;		//  Density of Water (1000 kg/m^3)
	const double w1day = 0.261799;		//  Daily frequency (2pi/24 hr 0.261799 radians/hr)

	//  Parameters
	const double es    = param(3-1);	//  emmissivity of snow (nominally 0.99)
	const double cg    = param(4-1);	//  Ground heat capacity (nominally 2.09 KJ/kg/C)
	const double z     = param(5-1);	//  Nominal meas. height for air temp. and humidity (2m)
	const double rho   = param(7-1);	//  Snow Density (Nominally 450 kg/m^3)
	const double rhog  = param(8-1);	//  Soil Density (nominally 1700 kg/m^3)
	double de          = param(11-1);	//  Thermally active depth of soil (0.1 m)
	LanS               = param(15-1);   // the thermal conductivity of fresh (dry) snow (0.0576 kJ/m/k/hr)
	LanG               = param(16-1);   // the thermal conductivity of soil (9.68 kJ/m/k/hr)

	const double wlf   = param(17-1);	//  Low frequency fluctuation in deep snow/soil layer (1/4 w1 = 0.0654 radian/hr)
	const double rd1   = param(18-1);	//  Amplitude correction coefficient of heat conduction (1)
	const double fstab = param(19-1);	//  Stability correction control parameter 0 = no corrections, 1 = full corrections

	zs = w*rhow/rho;

	rkin = rkinst(rkn, ws, tak, ts, z, g, fstab);

	// To make the UEB2 work under Unix/linux system the fancy stuff like "Select case" shall not be used
	//select case(iTsMethod)		Four choice of the surface temperature modeling
	//case (1)						1. Old snow, use calibrated snow surface heat conductance
	if (iTsMethod == 1) {
			qcs = dens*cs*rs*(ts - tavek);	// 2. Revised scheme LanE/Ze of the snow surface
											// 3. Force restore approach
		//case (2)                                 4. Modified force restore approach.
	} else if (iTsMethod == 2) {
		fKappaS = LanS/(rho*cs);
		fKappaG = LanG/(rhog*cg);

		d1 = sqrt(2*fKappaS/w1day);
		if (zs >= rd1*d1) {

			LanE_Ze = LanS/(rd1*d1);
		} else {
			 // call the subroutine to update the heat conductivity. LanE()
			LanE_Ze = LanE(LanS, LanG, zs, rho, rhog, cs, cg, rd1, ze, w1day);
			LanE_Ze = LanE_Ze/ze;
		}

			qcs = LanE_Ze*(ts - tavek);
	} else if (iTsMethod == 3) { //case (3)
			fKappaS = LanS/(rho*cs);
			fKappaG = LanG/(rhog*cg);

			d1 = sqrt(2*fKappaS/w1day);

			if (zs >= rd1*d1) {
				LanE_Ze = LanS/(rd1*d1);
				ze = rd1*d1;
	        } else {
				// call the subroutine to update the heat conductivity. LanE()
				LanE_Ze = LanE(LanS, LanG, zs, rho, rhog, cs, cg, rd1, ze, w1day);
				LanE_Ze = LanE_Ze/ze;
	        }

			de = ze/rd1;
			LanE_De = LanE_Ze/de*ze;
			qcs = LanE_De*(ts - ts_save::Ts_old)/(w1day*dt) + LanE_Ze*(ts - tavek);


		} else {	//case (4)   //change to all default cases. If not for model comparison
			fKappaS = LanS/(rho*cs);
			fKappaG = LanG/(rhog*cg);

			d1  = sqrt(2*fKappaS/w1day);
			dlf = sqrt(2*fKappaG/wlf);

			if (zs >= rd1*d1) {
				LanE_Ze = LanS/(rd1*d1);
				ze = rd1*d1;
			} else {
				// call the subroutine to update the heat conductivity. LanE()
				LanE_Ze = LanE(LanS, LanG, zs, rho, rhog, cs, cg, rd1, ze, w1day);
				LanE_Ze = LanE_Ze/ze;
			}

			if (zs >= rd1*dlf) {
				LanE_Ze2 = LanS/(rd1*dlf);
				Ze2 = rd1*dlf;
			} else {
				// call the subroutine to update the heat conductivity. LanE()
				LanE_Ze2 = LanE(LanS, LanG, zs, rho, rhog, cs, cg, rd1, Ze2, wlf);
				LanE_Ze2 = LanE_Ze2/Ze2;
			}

			de = ze/rd1;
			LanE_De = LanE_Ze/de*ze;
			de2 = Ze2/rd1;
			LanE_De2 = LanE_Ze2/de2*Ze2;

			if (ub <= 0.0 || refDepth <= 0.0) {
				qcs = LanE_De*(ts - ts_save::Ts_old)/(w1day*dt) + LanE_Ze*(ts - ts_save::Ts_Ave) +
					LanE_De2*(ts_save::Ts_Ave - ts_save::Tave_ave);
	        } else if (refDepth  >  rd1*d1) {
				qcs = LanE_De*(ts - ts_save::Ts_old)/(w1day*dt)+LanE_Ze*(ts - ts_save::Ts_Ave) +
					LanE_De2*(ts_save::Ts_Ave - ts_save::Tave_ave);
			} else {
                 qcs = LanE_Ze*ze*(ts - tk)/refDepth;
			}


			//End select
		}


       	result = qp + rkin*densa*cp*(tak - ts) + (hneu*densa*0.622*rkin)/pr*(ea - svp(ts - tk)) - qcs;
	    if (iradfl == 0) {
			result += qsn + qli - (1.0 - fc)*es*sbc*pow(ts, 4);
	    } else {
			result += qnetob;
	    }

	return result;
}


// *********************************************************************
double QcEst(const double ts,  const double rkn, const double ws, const double tak,
	const double g, const double qp,
	const double densa, const double cp, const double hneu,
	 const double pr, const double ea, const double tk, const double dens,
	 const double cs, const double rs, const double tavek, const double qsn,
	 const double qli, const double fc, const double sbc,
	 const double qnetob, const int iradfl, const ArrayXd &mtime,
	 const ArrayXd &param, const int iTsMethod, const double w, const double dt)
{
	//       function to estimate surface heat conduction for use in solving for
	//       surface temperature
	double result, rkin;
	const double es    = param(3-1);	//  emmissivity of snow (nominally 0.99)
	const double z     = param(5-1);	//  Nominal meas. height for air temp. and humidity (2m)
 	const double fstab = param(19-1);	//  Stability correction control parameter 0 = no corrections, 1 = full corrections

	rkin = rkinst(rkn, ws, tak, ts, z, g, fstab);
	result = qp + rkin*densa*cp*(tak - ts) + (hneu*densa*0.622*rkin)/pr*(ea - svp(ts - tk));
	if (iradfl == 0) {
		result += qsn + qli-(1.0 - fc)*es*sbc*pow(ts, 4);
	} else {
		result += qnetob;
	}

	return result;
}

/*
**************************** SRFTMPO () ***************************
C     Computes the surface temperature of snow
      FUNCTION SRFTMPO(QSI,A,QLI,QPIN,EA,TA,TAVE,TK,PR,RA,CP,RHO,RKN,
     *   HNEU,ES,SBC,CS,RS,W,qnetob,IRADFL,WS,Z,G,FC,fstab)
      dimension tint(0:10,2)
      DATA ncall,tol/0,0.05/
      ncall = ncall+1
      NITER = 10
      TSTAR = TA
      qp=qpin    ! store input variable locally without changing global value
      if (w <= 0. and. qp > 0.)qp=0.   ! ignore the effect of precip advected
C      energy on the calculation of surface temperature when there is no snow.
C      Without this ridiculously high temperatures can result as the model
C      tries to balance outgoing radiation with precip advected energy.
      TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *        RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab,mtime) !YJS pass modeling time
C      The first calculation gets the solution on the right side of ta
C      since it uses neutral stability and the linearization (equivalent
C      to Newton-Rhapson) will move in the direction of a solution in the case
C      of a well behaved function.  See Notebook 8/3/94 for further elaboration.
C      Use this to place bounds on surface temperature for iteration.
      if (tsurf > ta) {
        tlb=ta
        tub=ta+30.   !  Max upper bound 30 C warmer than surface
        if (tsurf > tub)tsurf=(tlb+tub)*.5  ! sometimes tsurf is outside these bounds
      else
        tlb=ta-30.
        tub=ta
        if (tsurf < tlb)tsurf=(tlb+tub)*.5
      endif
      tint(0,1)=tstar
      tint(0,2)=tsurf
C   Now iterate
      tstar=tsurf
      DO 10 I=1,NITER
      TSURF = SNOTMP(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *        RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab,mtime) !YJS pass modeling time

C      cerr << ta,tstar,tsurf
      tint(i,1)=tstar
      tint(i,2)=tsurf
      if (tlb <= tsurf && tsurf .le. tub) {
         if (tsurf > tstar) {
           tlb=tstar   ! increasing so can increase lower bound.
         else
           tub=tstar
         endif
       else if (tsurf  >  tub) {   ! upper bound overshot
         tlb=tstar                   ! increase lower bound and
         tsurf=(tstar+tub)*.5        ! guess at solution halfway
       else    ! tsurf < tlb  here, i.e. lower bound overshot
         tub=tstar
         tsurf=(tstar+tlb)*.5
       endif
C    Check for convergence
       if (ABS(TSURF-TSTAR) <  tol) THEN
          GO TO 20
       ELSE
          TSTAR = TSURF
       ENDIF
 10    CONTINUE
C       cerr <<  'slipped through the loop in SRFTMP()',ncall
C       do 15 i = 0,niter
C 15      cerr << tint(i,1),tint(i,2)
C       cerr << "tsurf",tsurf
C    Newton rhapson not converging so use bisection
       f1=sntmpb(tlb,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
       f2=sntmpb(tub,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
       tsurf=(tlb+tub)*.5
       if (f1*f2  > 0.) {
         cerr << 'SRFTMP has failed to find a solution',ncall
     *     ,tlb,tub,tsurf
       else
         nib=(alog(tub-tlb)-alog(tol))/alog(2.)
         do 16 i=1,nib
C         cerr << nib,tlb,tub,tsurf
         f=sntmpb(tsurf,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)
         if (f*f1  > 0.) {
           tlb=tsurf
           f1=f
         else
           tub=tsurf
C          f2=f
         endif
         tsurf=(tlb+tub)*.5
 16      continue
C         cerr << "Bisection",nib,tsurf
       endif
 20     if (W > 0..AND.TSURF > 0.) THEN
          SRFTMPO = 0.0
       ELSE
          SRFTMPO = TSURF
       ENDIF
C       cerr << "Final",srftmp
       RETURN
       END

**************************** SNOTMP () ******************************
C     Function to compute surface temperature using Penman/surface
C     resistance analogy and equilibrium approach

      function snotmp(tstar,qsi,a,qli,qp,ea,ta,tave,tk,pr,ra,cp,rho,
     *        rkn,hneu,es,sbc,cs,rs,qnetob,iradfl,ws,z,g,fc,fstab,mtime)

	real LanmdaE, LanS, LanG, LanE
      qsn = qsi*(1.0-a)
C     To ensure all temperatures in kelvin
      tak = ta+tk
      tstark = tstar+tk
      tavek  = tave+tk
      densa  = pr/(ra*tak)
c     densa in kg/m3 if pr in pa
c      cp1 = cp/1000.0
c     cp1 in kj/kg/ok


      DENS = RHO
      RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab)

	LanmdaE=LanE(LanS, LanG, Zs, DENS, rhog, cs, cg, r, ze, w1day)

      UPPER = QP+(DENSA*CP*TAK)*RKIN
     *        -(HNEU*DENSA*0.622*RKIN)/PR*(SVPI(TSTAR)
     *        -EA-DELTA(TSTAR)*TSTARK)
     *        +LanmdaE*TAVEK                     !change to combined layer (snow and soil)
      DEN = LanmdaE+(DENSA*CP)*RKIN+(DELTA(TSTAR)*
     *       HNEU*DENSA*.622*RKIN)/PR
      if (IRADFL.EQ.0) THEN
         UPPER = UPPER + QSN+QLI+3.0*(1.-FC)*ES*SBC*TSTARK**4
         DEN = DEN+4.0*(1.-FC)*ES*SBC*TSTARK**3
      ELSE
         UPPER = UPPER + qnetob
      ENDIF

      SNOTMP = UPPER/DEN-TK

      RETURN
      END
C*************************** sntmpb () ******************************
C     Function to compute surface temperature using bisection

      FUNCTION sntmpb(TSTAR,QSI,A,QLI,QP,EA,TA,TAVE,TK,PR,RA,CP,RHO,
     *         RKN,HNEU,ES,SBC,CS,RS,qnetob,IRADFL,WS,Z,G,FC,fstab)

      QSN = QSI*(1.0-A)
C     To ensure all temperatures in kelvin
      TAK = TA+TK
      TSTARK = TSTAR+TK
      TAVEK  = TAVE+TK
      DENSA  = PR/(RA*TAK)
C     DENSA in kg/m3 if PR in Pa
C      CP1 = CP/1000.0
C     CP1 in kj/kg/oK
      DENS = RHO
      RKIN=RKINST(RKN,WS,TAK,TSTARK,Z,G,fstab)

      UPPER = QP+(DENSA*CP*TAK)*RKIN
     *        -(HNEU*DENSA*0.622*RKIN)/PR*(SVPI(TSTAR)
     *        -EA)
     *        +DENS*CS*RS*TAVEK
      DEN = DENS*CS*RS+(DENSA*CP)*RKIN
      if (IRADFL.EQ.0) THEN
         UPPER = UPPER + QSN+QLI-(1.-FC)*ES*SBC*TSTARK**4
      ELSE
         UPPER = UPPER + qnetob
      ENDIF
      sntmpb=upper-den*tstark
      RETURN
      END
*/
// ******************************* RKINST() **********************
double rkinst(const double rkn, const double ws, const double ta,
	const double ts, const double z, const double g, const double fstab)
{
	//     function to calculate no neutral turbulent transfer coefficient using the
	//     richardson number correction.
	double result, rich;
	if (ws <= 0.0) {
		result = 0.0;    //  No wind so no sensible or latent heat fluxes.
	} else {
		rich = g*(ta-  ts)*z/(ws*ws*ta);    // ta must be in K
		if (rich >= 0.0) {
			result = rkn/(1.0 + 10.0*rich);
			// PhiMH=1+10*rich    Change back
		} else {
			//          result=rkn*(1-10.*rich)
			//PhiMH=(1-16*rich)**(-.75)
			result = rkn*min(3.0, pow(1.0-16.0*rich, 0.75));
		}
	}
	//  Linear damping of stability correction through parameter fstab
	result = rkn+fstab*(result-rkn);
	// if (PhiMH != 0.0) result=rkn/PhiMH Change back to the original design purpose

	return result;
}


// ******************************* DELTA () ****************************
//     Function to compute gradient of saturated vapour pressure,
//     temperature function over ice
//     Uses Lowe (1977) polynomial
double delta(const double t)
{
	double result;
	double a[7] = { 0.5030305237, 0.0377325502, 0.001267995369,
				2.477563108e-5, 3.005693132e-7, 2.158542548e-9,
				7.131097725e-12 };

	result = a[0] + t*(a[1] + t*(a[2] + t*(a[3] + t*(a[4] + t*(a[5] + t*a[6])))));
	result *= 100.0;    // convert from mb to Pa

	return result;
}

// ************************** PARTSNOW () ************************
//     Partitioning of precipitation into rain and snow

double partsnow(const double p, const double ta, const double tr, const double ts)
{
	double result;
	if (ta < ts) {
		result = p;
	} else if (ta > tr) {
		result = 0.0;
	} else {
		result = p*(tr - ta)/(tr - ts);
	}

	return result;
}


// *************************** QPF () ***************************
//     Calculates the heat advected to the snowpack due to rain
double qpf(const double pr, const double ta, const double to, const double ps,
	const double rhow, const double hf, const double cw, const double cs)
{
	double tsnow, train;
	if (ta > to) {
		train = ta;
		tsnow = to;
	} else {
		train = to;
		tsnow = ta;
	}

	return pr*rhow*(hf + cw*(train - to)) + ps*rhow*cs*(tsnow - to);
}


// ************************ TURBFLUX () ***************************
//     Calculates the turbulent heat fluxes (sensible and latent
//     heat fluxes) and condensation/sublimation.
int turbflux(const double pr, const double ra, const double ta, const double tk,
	const double ts, const double z, const double g, const double cp,
	const double rkn, const double ws, double &ea, const double rhow,
	const double hneu, double &qh, double &qe, double &e, const double fstab)
{

	double rhoa, tak, tsk, rkin, es;

	tak = ta + tk;
	tsk = ts + tk;
	rkin = rkinst(rkn,ws,tak,tsk,z,g,fstab);
	rhoa = pr/(ra*(tak));
	//     rhoa in kg/m3
	qh = rhoa*(ta-ts)*cp*rkin;
	es = svpi(ts);
	qe = 0.622*hneu/(ra*(tak))*rkin*(ea - es);
	e = -qe/(rhow*hneu);	//     e in  m/hr

	return 0;
}


// *********************** FMELT () *************************
//     Calculates the melt rate and melt outflow
double fmelt(const double ub, const double rhow, const double w, const double hf,
	const double lc, const double rid, const double ks, const double pRain)
{
	double result, uu, ss;

	uu = 0.0;
	if (ub < 0.0) {
		result=0.0;
	} else if (w <= 0.0) {
		result = pRain;
		if (pRain <= 0.) {
			result = 0.0;
		}
	} else {
		uu = ub/(rhow*w*hf);
		//                              liquid fraction
		if (uu > 0.99) {
			uu = 0.99;
			//                              TO ENSURE HIGH MELT RATE WHEN ALL LIQUID
		}
		if ((uu/(1.0 - uu)) <= lc) {
			ss = 0.0;
		} else {
			ss = (uu/((1.0 - uu)*rid) - lc/rid)/(1.0 - lc/rid);
		}
		result = ks*pow(ss, 3);
	}
	if (result < 0.0) {
		cerr << "fmelt is negative\n";
		exit(EXIT_FAILURE);
	}

	return result;
}


// *************************** SVP () *****************************
//     Calculates the vapour pressure at a specified temperature over water or ice
//     depending upon temperature.  Temperature is celsius here.
double svp(const double t)
{
	double result;

	if (t >= 0.0) {
		result = svpw(t);
	} else {
		result = svpi(t);
	}

	return result;
}

// *************************** SVPI () *****************************
//     Calculates the vapour pressure at a specified temperature over ice.
//     using polynomial from Lowe (1977).
double svpi(const double t)
{
	double result;

	result = 6.109177956 + t * (0.503469897 + t*(0.01886013408 +
		t * (0.0004176223716 + t * (5.82472028e-06 + t *
		(4.838803174e-08 + t * 1.838826904e-10)))));
	result *= 100.0;   // convert from mb to Pa

	return result;
}

// *************************** SVPW () *****************************
//     Calculates the vapour pressure at a specified temperature over water
//     using polynomial from Lowe (1977).
double svpw(const double t)
{
	double result;

	result = 6.107799961 + t * (0.4436518521 + t * (0.01428945805 +
		t * (0.0002650648471 + t * (3.031240936e-06 + t *
		(2.034080948e-08 + t * 6.136820929e-11)))));
	result *= 100.0;   // convert from mb to Pa
	return result;
}

// **************************** PREHELP () ***************************
//      Routine to correct energy and mass fluxes when
//      numerical overshoots dictate that W was changed in
//      the calling routine - either because W went negative
//      or due to the liquid fraction being held constant.

int prehelp(const double w1, const double w, const double dt, double &fm,
	const double fm1, const double fac, const double ps, const double pRain,
	double &e, const double rhow, const double hf, double &q, double &qm,
	double &mr, double &qe, const double hneu)
{
	double qother;
	//   The next statement calculates the value of FM given
	//    the W and w1 values
	fm = (w1 - w)/dt*fac-fm1;
	//   The next statement calculates the changed MR and E due to the new FM.
	//    FM was = PRAIN+PS-MR-E
	//    Here the difference is absorbed in MR first and then E if mr < 0.

	mr = max(0.0, (ps + pRain - fm - e));
	e = ps + pRain - fm - mr;
	//    Because melt rate changes the advected energy also changes.  Here
	//     advected and melt energy are separated,
	qother = q + qm - qe;
	//     then the part due to melt recalculated
	qm = mr*rhow*hf;
	//     then the part due to evaporation recalculated
	qe = -e*rhow*hneu;
	//     Energy recombined
	q = qother - qm + qe;

	return 0;
}

// ******************************** ALBEDO () *****************************
//     Function to calculate Albedo
//     BATS Albedo Model (Dickinson et. al P.21)
double albedo(const double tausn, const double coszen, const double d,
	const double aep, const double abg, const double avo, const double airo)
{
	double b, cs, cn, fage, result, fzen, rr, avd, avis, aird, anir;
		b = 2.0;
		cs = 0.2;
		//      avo = 0.95
		cn = 0.5;
		//      airo = 0.65

		fage = tausn/(1.0 + tausn);

		if (coszen < 0.5) {
			fzen = 1.0/b*((b + 1.0)/(1.0 + 2.0*b*coszen) - 1.0);
		} else {
			fzen = 0.0;
		}
		avd = (1.0 - cs*fage)*avo;
		avis = avd + 0.4*fzen*(1.0 - avd);
		aird = (1.0 - cn*fage)*airo;
		anir = aird + 0.4*fzen*(1.0 - aird);
		result = (avis + anir)/2.0;
		if (d < aep) {   // need to transition albedo to a bare ground value
			rr = (1.-d/aep)*exp(-d*0.5/aep);
			result =rr*abg + (1.0 - rr)*result;
		}

		return result;
}

// ******************************** AGESN () *****************************
//     Function to calculate Dimensionless age of snow for use in
//     BATS Albedo Model (Dickinson et. al P.21)

int agesn(double &tausn, const double dt, const double ps,
	const double tsurf, const double tk, const double dnews)
{
	double tsk, r1, r2, r3;

	tsk = tsurf + tk;  // express surface temperature in kelvin
	r1 = exp(5000.0*(1.0/tk - 1.0/tsk));
	r2 = pow(r1, 10);
	if (r2 > 1.0) {
		r2 = 1.0;
	}
	r3 = 0.3;
	//   Dickinson p 23 gives decay as DT*(R1+R2+R3)/tau0  with
	//    tau0 = 10**6 sec.  Here 0.0036 = 3600 s/hr * 10**-6 s**-1
	//    since dt is in hours.
	//      tausn = max((tausn+0.0036*(R1+R2+R3)*DT)*
	//     &            (1.0 - 100.0*PS*DT),0.)
    tausn = tausn + 0.0036*(r1 + r2 + r3)*dt;
	if (ps > 0.0) {
		if (dnews  >  0.0) {
			tausn = max(tausn*(1.0 - ps*dt/dnews), 0.0);
		} else {
			tausn = 0.0;
		}
	}

	return 0;
}


// ***************************** LanE () *************************************
//	Function to get the LanE which is the thermal conductivity by ze

double LanE(const double lans, const double lang, const double zs, const double rho,
	const double rhog, const double cs, const double cg, const double r, double &ze, const double w1day)
{
	double fkappas, fkappag, d1, d2, result;

	fkappas = lans/(rho*cs);
	fkappag = lang/(rhog*cg);
	d1 = sqrt(2*fkappas/w1day);
	d2 = sqrt(2*fkappag/w1day);
	result = min(r, zs/d1)*d1/lans + max(0.0, (r - zs/d1))*d2/lang;
	ze = min(r, zs/d1)*d1 + max(0., (r - zs/d1))*d2;

	if (result != 0.0)
		result = 1.0/result*ze;

	return result;
}


// *************************** refDep() **************************************
//     function to calculate the refreezing depth
double refDep(const double flans, const double a, const double b, const double hf,
	const double rhom, const double dt, const double x1 )
{
	double result, temp;

	temp = flans*flans + 2*b*(-a*flans*dt/(hf*rhom) + flans*x1 + 0.5*b*x1*x1);

	if (temp < 0.0 || a > 0.0 || b == 0.0) {
		result = 0.0;
	} else {
		result = (-flans + sqrt(temp))/b;
	}

	return result;
}

// ************************** Grad () ********************************
// Linear simple function to calculate the gradient of Q and T

int Grad(const double qc1, const double qc2, const double t1, const double t2, double &a, double &b)
{
	if ((t2 - t1) != 0.0) {
		b = (qc2 - qc1)/(t1 - t2);
		a = qc1 + b*t1;
	}

	return 0;
}

// ********************* Caculate the daily average  *********************
//  Calculate the daily average value
//  n number of records, a minimum value -100 or some
double daily_ave(const ArrayXd &backup, const int n, const double a)
{
	int i;
	double sum = 0.0, result;
	double count = 0.0;
	for (i = 0; i < n; i++) {
		if (backup(i) > a) {
			sum += backup(i);
			count += 1.0;
		}
	}
	if (count != 0.0) {
		result = sum/count;
	} else {
		result = 0.0;
	}

	return result;
}

