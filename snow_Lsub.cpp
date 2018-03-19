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

//   Lumped snowmelt module can be called outside  (With depletion curve)
//   input of snow module include
// year,month,day,hour,DT,nStep -- year, month, day, hour,  time step (hour), number of time step
// INPT  -- input forcing
// 	  inpt(1,1)=t          ! temperature to subwatershed (C)
// 	  inpt(2,1)=rain(j)    !precipitation
// 	  inpt(3,1)=ws(j)      !wind speed  (m/s)
// 	  inpt(4,1)=rh         !relative humidity
// 	  inpt(5,1)=qsi(j)     !short wave  (kJ/m^2/h)
// 	  inpt(6,1)=qnet(j)    !net radiation  (kJ/m^2/h)
// 	  inpt(7,1)=trange(j)  !temperature range  (C, it will be updated to zeinth
//           angle while calling the UEB point snowmelt module)
// SITEV -- site variables
//        site variables (1-5)
//        sitev(1)=0                     !forest cover fraction
//        sitev(2)=1.0                   !drift factor  (No detailed information give 1)
//        sitev(3)=pa                    !air pressure
//        sitev(4)=3.6                   !ground heat flux  Assign a constant 1 W/m^2
//        sitev(5)=0.1                   !albedo extinction parameter (veghght)

// STATEV	1:5 point state variables (U, W, dimensionless age, refdepth, totalrefdepth)
//        plus six depletion curve lumped parameters 6:11 (Wmax, areafrac, meltflag,Wtgt,Wtgtmax,Aftgt)

// PARAM  --  snowmelt module parameters (as above)
// iflag  -- flags to call point snowmelt module
// 	   iflag(1) = irad   ! on input this has value
//         	0 is no measurements - radiation estimated from diurnal temperature range
// 			1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
// 			2 is incoming shortwave and longwave radiation read from file (measured)
// 			3 is net radiation read from file (measured)
//          internally iflag(1) is passed to snowueb with the meaning
//          0=radiation is shortwave in col 5 and longwave in col 6, else = net radiation in column 7
//        iflag(2) = 0      ! no(/yes 1) printing
//        iflag(3) = 9      ! unit to which to print
//        iflag(4) = 1      ! how albedo calculations are done
// 	   iflag(5) = 4      ! model option for surface temperature approximation (modified
//            force-restorce approach (Luce, 2000; Luce and Tarboton, 2001; You et al.  ;
//            You, 2002))
// dtbar  -- monthly mean daily temperature range
// cump,cume,cummr  -- cumulative precipitation (with df), cumulative evaporation, cumulative melt rate
// outv  -- output variables (from point snowmelt modeling (UEB) 13 variables, see snow menu (Tarboton, 1995,
//      Tarboton and Luce (1996))
// tsbackup,tavebackup  --  temperature storage of modeled snow surface temperature, snow average
//      temperature in last 24 hours (48 means 0.5 hour time step is possible)( for modified
//      force restore approach
// dfc  -- depletion curve (wa/wmax, afraction)
// baout,out  -- (basin output and point snowmelt output)

int snowLSub(int &year, int &month, int &day, double &hour, const int dt,
	const int nStep, ArrayXXd &inpt, ArrayXd &sitev, ArrayXd &statev,
	const ArrayXd &param, ArrayXi &iflag, const ArrayXd &dtbar, const int nstepday,
    double &cump, double &cume, double &cummr, ArrayXd &outv,
    ArrayXd &tsbackup, ArrayXd &tavebackup,
    const int ndepletionpoints, double **dfc, const int modelelement, const int jj)
{
	double lat;
	ArrayXd mtime(5);      // YJS pass to reflect the change of Snow
	// real dfc(ndepletionpoints,2)

	// Constant data set
	const double io    = 4914.0;		//  solar constant  kJ/m^2/hr
    const double tk    = 273.15;		//  Temperature to convert C to K (273.15)
    const double sbc   = 2.0747e-7;		//  Stefan boltzman constant (2.0747e-7 KJ/m^2/hr/K)
    const double hf    = 333.5;			//  Heat of fusion (333.5 KJ/kg)
    const double cs    = 2.09;			//  Ice heat capacity (2.09 KJ/kg/C)
    const double rhow  = 1000.0;		//  Density of Water (1000 kg/m^3)
    const double bca   = 0.8; 			//  A in Bristow Campbell formula for atmospheric transmittance
    const double bcc   = 2.4; 			//  C in Bristow Campbell formula for atmospheric transmittance

	double afrac, baw, wmax, meltflag, wtgt, wtgtmax, aftgt, slope, azi;
	double afnew, woldsca, wnew, wnewsca, scaw, basub, cg, rhog, de, atff, hri0;
	double afracnew, baswin, rh, qsiobs, qnetob, trange, ta, hri, coszen;
	double atfimplied, cf, iradfl, deltaw, dw1, dw2, r, wref;
	int irad, i, iit;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> snowLSub(" << ncalls << ")" << std::endl;
    }
	caller = "snowLSub";
#endif
	// Unravel state variables
	afrac    = statev(6);
	baw      = statev(1);
	wmax     = statev(5);
	meltflag = statev(7);
	wtgt     = statev(8);     // target snow water equivalent for snowfall hysteresis loops during melt
	wtgtmax  = statev(9);
	aftgt    = statev(10);
	irad     = iflag(0);
	slope    = sitev(5);
	azi      = sitev(6);
	lat      = sitev(7);
	// debugging code to stop at a specific date
	//       if ((year .eq. 1960) && (month .eq. 10) &&
	//      + (day .eq. 27) && ! (hour .gt. 8) &&
	//     + (modelelement .eq. 110))then
	// 	   year=year
	// 	endif
double afracStart = afrac;
    for (iit = 1; iit <= nStep; iit++) {    // time loop
		if (afrac > 0.0 && baw > 0.0) {
			//   See discussion in "depletion curve logic.doc".  Revised to avoid spikes this was causing.
			//     The idea here is to use a weighted average of the snow water equivalent depths over the old and
			//     new snow covered area as the representative snow water equivalent to give the point model while
			//     on an excursion.
			if (meltflag == 0.0 && aftgt > 0.0 && afrac > aftgt && baw > wtgt) {  //  here we are on an excursion
				// 			scaw=wtgt/aftgt+baw-wtgt   ! dgt 7/25/05.  initial approach superceded by the following
				afnew = afrac - aftgt;
				woldsca = wtgt/aftgt;
				wnew = baw - wtgt;
				wnewsca = wnew/afnew;
				scaw = (wnewsca*wnew + woldsca*wtgt)/baw;
			} else {
				scaw = baw/afrac;
			}
		} else { // initialize for no snow
			scaw     = 0.0;
			afrac    = 0.0;
			baw      = 0.0;
			wmax     = 0.000001;   //  dgt 7/25/05 to avoid nan when first used
			wtgtmax  = 0.0;
			meltflag = 1.0;
		}
		statev(1) = scaw;

		//    DGT 7/26/05.  Statev(1) in the lumped model stores the average temperature except in the case
		//    when liquid water is present in which case it stores the liquid fraction.
		//    The assumption is that these quantities remain constant during depletion curve adjustments to W
		//    and are used to reinstate a physically realistic U.
		cg   = param(4-1);    //  Ground heat capacity (nominally 2.09 KJ/kg/C)
		rhog = param(8-1);   //  Soil Density (nominally 1700 kg/m^3)
		de   = param(11-1);    //  Thermally active depth of soil (0.1 m)

		if (statev(0) < 0.0) {
			//  Here energy content state variable was temperature and needs to be converted to energy
			statev(0) = statev(0)*(rhog*de*cg + rhow*scaw*cs);
		} else {
			if (scaw > 0.0) {
			//  Here energy content state variable was liquid fraction
				statev(0) = statev(0)*rhow*scaw*hf;
			} else {
				//  Here energy content variable was temperature.  scaw is 0 so energy is only of soil
				statev(0) = statev(0)*rhog*de*cg;
			}
		}

		//    ***************** This is start of point model
		//
		//    From the  input list
		rh     = inpt(4-1,0);
		qsiobs = inpt(5-1,0);
		qnetob = inpt(6-1,0);
		trange = inpt(7-1,0);
		ta     = inpt(1-1,0);
		hyri(year, month, day, hour, dt, slope, azi, lat, hri, coszen);
		inpt(7-1,0) = coszen;
		if (irad <= 1) {
			atf(atff, trange, month, dtbar, bca, bcc);
			if (irad == 0) {  // need to estimate radiation based on air temp.
				inpt(5-1,0) = atff*io*hri;
			} else {
				//       Need to call HYRI for horizontal surface to perform horizontal
				//       measurement adjustment
				hyri(year, month, day, hour, dt, 0.0, azi, lat, hri0, coszen);
				//       If HRI0 is 0 the sun should have set so QSIOBS should be 0.  If it is
				//       not it indicates a potential measurement problem. i.e. moonshine
				if (hri0 > 0.0) {
					atfimplied = min(qsiobs/(hri0*io), 0.9); // To avoid unreasonably large radiation when hri0 is small
					inpt(5-1,0) = atfimplied*hri*io;
				} else {
					inpt(5-1,0) = qsiobs;
					if (qsiobs != 0.0) {
						cerr << "Warning, you have nonzero nightime";
						cerr << " incident radiation of " << qsiobs;
						cerr << " at date " << year << "/" << month << "/" << day << ", " << hour << '\n';
					}
				}
			}
			//         cloud cover fraction dgt 10/13/94
			cf = 1.0 - atff/bca;
			qlif(inpt(6-1,0), ta, rh, tk, sbc, cf);
			iradfl = 0.0;
		} else {
			iradfl = 1.0;
			inpt(5-1,0) = qnetob;
		}
		iflag(1-1) = iradfl;
		updatetime(year, month, day, hour, dt);

		if (iflag(2-1) == 1) {
			cout << "iflag(2) " << iflag(1) << " iflag(3) " << iflag(3-1) << endl;
			//for (i = 0; i < 7; ++i) {
                //cout <<
		   //write(iflag(3),*)year,month,day,hour,atff,hri,(inpt(i-1,0),i=1,7)
		   //}
		}

		//     set modeling time
		mtime(0) = (double)year;
		mtime(1) = (double)month;
		mtime(2) = (double)day;
		mtime(3) = (double)hour;
		mtime(4) = modelelement;   // This used to know where we are for debugging

		//    cumulative variables used for single time step over snow covered area
		cump  = 0.0;
		cume  = 0.0;
		cummr = 0.0;
		snowueb2(dt, 1, inpt, sitev, statev, tsbackup, tavebackup, nstepday, param, iflag,
			cump, cume, cummr, outv, mtime, modelelement, jj);   //Add a pass varible to snow mtime

		//  ************************ End of point model

		//      DeltaW = Change in water equivalent
		deltaw = statev(2-1) - scaw;
		if (outv(2-1) > 0.0 && deltaw > 0.0) {   // Here an increase in snow water equivalent due to snowfall
			if (meltflag == 1) {
				meltflag = 0.0;   //  start new excursion
				wtgt = baw;
				aftgt = afrac;
				wtgtmax = wtgt;  // establish new upper limit of new excursion
			}
			baw = baw + deltaw;
			wtgtmax = max(baw, wtgtmax);
			if (baw > wmax) {
				wmax = baw;
				meltflag = 1.0;  // end up on main curve if wmax increases
			} else if (baw > 3.0*wtgt) {   //  DGT 7/25/05   If new snow dominates old snow we
				//   reset the process to be on the depletion curve.
				//   Here domination is arbitrarily defined as greater than 3 times.
				meltflag = 1.0;
				wmax = baw;
			}
			afracnew = 1.0;
			baswin = cummr;    // basin average surface water input
			basub = cume;   // basin average sublimation
		} else if ( deltaw > 0) {  // Here a rare case when W increases but there is no snowfall.
			//    This may be due to rainfall absorbing and freezing into the snow or condensation.
			//     This is assumed to only be effective over snow covered area fraction
			baw += deltaw*afrac;
			wtgtmax = max(baw, wtgtmax);
			if (baw > wmax) {
				wmax = baw;
				meltflag = 1.0;   // end up on main curve if wmax increases
			}
			//    DGT 7/25/05  Logically the dominates check above should be done here too, but it is not
			//    because if Wmax is reset Afrac should go to 1 to put the process at the top of a depletion curve
			//    and we are not increasing Afrac because there is no snowfall.  The implication is that baw may
			//    go greater than 3 * Wtgt and the depletion curve only get reset next time snow falls.
			//    The 3 is arbitrary so this does not seem like a serious shortcoming.
			afracnew = afrac;   // (emphasize does not change)
			//          One could increase Wmax here to keep on main depletion curve.  However I decided to handle this issue below with the Min statements so that Af never increases when there is melt.  Excursions below the depletion curve will therefore return back to the same Wa before further reductions in Af occur.
			baswin = cummr*afrac + cump*(1.0 - afrac);  //  All precipitation is surface water input over snow free area
			basub = cume*afrac;  // Sublimation is only counted over snow covered area.  Calling model has to handle ET from soil
		} else {   // Here (DeltaW < 0) so there is melt
			if (meltflag == 0 && (baw - deltaw*afrac) < wtgt) {
				//    Here reduction passes target so do in two steps to avoid overshooting and large steps towards end of season
				dw1 = -(baw - wtgt)/afrac;
				dw2 = -(deltaw - dw1);
				//              now logically in two steps
				//                 baw=baw+DW1*Afrac+Aftgt*DW2
				//              which simplifies to
				baw = wtgt + aftgt*dw2;
			} else {
				baw += deltaw*afrac;
			}
			if (baw < 0) {
				//       There is a mass balance issue here because baw is insufficient to supply M and E that
				//       result in the reduction DeltaW.
				r = 1.0 - (baw/afrac*deltaw);   // This assumes the shortfall will be evenly spread between M and E
				baw = 0.0;
				cummr = cummr*r;
				cume = cume*r;
				//        DGT 7/26/05  Since snow is disappearing here we should make sure energy behaves properly
				//        Set energy content to 0.  There was snow so any energy was in liquid fraction not temperature greater
				//        than 0 and that is now gone.
				if (statev(0) >= 0.0) {
					statev(0) = 0.0;
					// 			else
					// 				statev(1)=statev(1)    ! redundant but debugging to see if this ever occurs
				}
			}
			baswin = cummr*afrac + cump*(1.0 - afrac);
			basub = cume*afrac;
			if (baw < wtgt)
				meltflag = 1;  // Back on main melt curve
			if (meltflag == 1) {
				afracnew = min(ablat(baw/wmax, ndepletionpoints, dfc), afrac);   // Regular depletion curve.
				// Make sure that Af does not increase in this step for special case where W may have increased
				// to be below depletion curve by absorbing of precip or snow.  Assume Af remains the same until W returns to original Af
			} else {
				wref = wmax - (wtgtmax - baw)/(wtgtmax - wtgt)*(wmax - wtgt);  // See depletion curve logic.doc for derivation and picture
				afracnew = min(ablat(wref/wmax, ndepletionpoints, dfc), afrac);
			}
		}
		//   DGT 7/25/05.  Save in statev(1) the energy state variable that persists through depletion curve
		//    baw and scaw adjustments
		//
		if (statev(0) <= 0.0) {
			//  Here convert energy to temperature - using statev(2) output from the point model prior to
			//  depletion curve adjustments
			statev(0) = statev(0)/(rhog*de*cg + rhow*statev(1)*cs);
		} else {
			if (baw > 0.0) {
				//   Here snow water equivalent is positive and energy content positive so
				//   energy content variable is the liquid fraction
				//  Note - this is a strictly greater than.  A tolerance could be considered
				//  Note also that the test is on baw because depletion curve adjustments may result in this
				//  being 0 when statev(2) is not so then the temperature inter
				if (statev(1) <= 0.0) {
					if (baw > 1.0e-6) {
						//   Due to arithmetic precision in the division and multiplication by afrac
						//    baw > 0 can occur when statev(2) = 0.  Only consider this an error if > 1e-6
						cerr << "error, baw > 0, statev(2) <= 0.0\n";
						baw = 0.0;
						statev(0) = 0.0;
					} else {
						baw = 0.0;
						//   Here energy content is positive and snow water equivalent is 0 so
						//   energy content state variable is the average soil temperature
						statev(0) = statev(0)/(rhog*de*cg);
					}
				} else {
					statev(0) = statev(0)/(rhow*statev(1)*hf);
				}
			} else {
				if (statev(1) > 0.0) {
					//   Here baw was 0 but statev(2) is positive.  This can occur due to depletion curve adjustments
					//   Treat this as a case where snow water equivalent is 0 but disregard the previous energy content
					//   because it is greater than 0 and reflects liquid water fraction which has been removed.
					statev(0) = 0.0;
					cerr << "Liquid fraction energy content set to 0\n";
					cerr << "because depletion curve set W to 0\n";
				} else {
					//   Here energy content is positive and snow water equivalent is 0 so
					//   energy content state variable is the average soil temperature
					statev(0) = statev(0)/(rhog*de*cg);
				}
			}
		}
		//    Output variables
		outv(15-1) = baw;
		outv(16-1) = afrac;
		if (afracnew == 1)
			outv(16-1) = afracnew;  // output old Afrac unless in increased because old Afrac was used in mass balances
		afrac = afracnew;
		outv(17-1) = meltflag;
		outv(18-1) = wtgt;
		outv(19-1) = wtgtmax;
		outv(20-1) = wmax;
		outv(21-1) = aftgt;
		outv(22-1) = baswin;
		outv(23-1) = basub;
		if (iflag(2-1) == 1) {
			cout << iflag(3-1) << " ";
			for (i = 15; i <= 23; i++) {
				cout << outv(i-1) << " ";
			}
		}
	}// end of time step loop

	iflag(1-1)   = irad;  // reinstate iflag(1) for next call
	statev(2-1)  = baw;
	statev(6-1)  = wmax;
	statev(7-1)  = afrac;
	statev(8-1)  = meltflag;
	statev(9-1)  = wtgt;
	statev(10-1) = wtgtmax;
	statev(11-1) = aftgt;
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving snowLSub(" << ncalls << ")" << "\n" << endl;
    }
    ncalls++;
#endif

	return 0;
}


// *******************************************************************

double ablat(const double DimBaw, const int ndfc, double **dfc)
{
	//       Function to interpolate dimensionless depletion curve
	double Afrac;
	int i;
	//     DimBaw is dimensionless Basin Average Water Equivalent corresponding to column (1) of dfc
	//     ndfc is number of points rows in dfc
	//     dfc is dimensionless (0-1) basin average snow water equivalent (col 1) and snow covered area fraction (0-1) col 2

	if (DimBaw >= 1.0) {
		Afrac = 1.0;
	} else {
		if (DimBaw <= 0.0) {
			Afrac = 0.0;
		} else {
			for (i = 1; i <= ndfc; i++) {
				if (i == 1) {
					Afrac = dfc[0][1];   // DGT 7/25/05 to avoid Afrac ever being interpolated as 0 and avoid divide by 0 errors
				} else {
					if (dfc[i-1][0] >= DimBaw) {
						Afrac = dfc[i-1-1][2-1] + (dfc[i-1][2-1] - dfc[i-1-1][2-1])*(DimBaw - dfc[i-1-1][1-1])/(dfc[i-1][1-1] - dfc[i-1-1][1-1]);   // assumes points not repeated in the first column of dfc
						goto L1000;
					}
				}
			}
		}
	}

L1000: if (Afrac >= 1.0)
		Afrac = 1.0;

	return Afrac;
}
