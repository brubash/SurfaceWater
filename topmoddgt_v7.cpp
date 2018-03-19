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

using namespace Eigen;
using namespace std;

// warning, Ross has subverted the meaning of qinst for nooksack

// This version, V2, has been corrected for -ve seepage under saturated conditions.

// *******************************************************************
//
//  Program to use discrete distribution Z-based TOPMODEL
//
// ********************************************************************
// RPI 22/3/2004 - a lot of these variables in V6 were automatic arrays, i.e.
// the space allocated to them disappeared when the subroutine was exited.
// This did not matter in V6 as they were not required outside this routine.
// However, in V7 they have to retain their values between calls and so must
// be made explicit shape arrays. One way to do this is to put them in the
// calling sequence. Also there is scope for avoiding repeated calculations
// at the start of the routine - this will become a big problem on long runs
// as many unnecessary calculations will end up being made. To overcome this problem

int intercept(double &cv, const double rit, double &ae, const double delt, const double x1, const double cr, double &r1, double &ud);
int soil(const double sr, const double r2, double &r3, double &rie, double &rd, double &srn,
	const double zr, const double dth1, const double dth2, const double szf,
	const double ak0fzrdt, const double c, const double soilc,
	const double psif, const bool global, const double ad_cap, double &ad_now);

int topmod(double **si, const ArrayXXd &Sp, const int isub, const int *Nka, const double Lambda,
	double **atb, double **pka, const int *nd, double **cl, double **pd, const double units,
	Array<int,Dynamic,1> &irr, const bool modwrt, const int ipsub, const int ipatb, const int stim,
	const double r, const double pet, const long int interval, const double art_drainage, const double rate_irrig,
	const int imonth, const int ndata, const int mps, const int mpe, double &qinst_out, double &dr_out,
	const int ndump, int *ntdh, const int istep, const int maxC, double &zbm, const int maxA,
	const int maxSlp, const int maxInt, double &sumr, double &sumq, double &sumae,
	double &s0, double &q0, double &sr, double &cv, double &aciem,
	double &acsem, double &sumpe, double &sumie, double &sumqb, double &sumce,
	double &sumsle, double &sumr1, double &qb, ArrayXd &qinst, ArrayXd &dr, double &sumqv,
	double &sumse, double &zbar, const double zbar_new, double **tdh, double &zr, double &ak0fzrdt, double &logoqm,
	double &qvmin, double &dth, double &sumad, double &evap_mm, double &qlat_mm,
	const int ipflag, ArrayXd &rirr, const int js)
{
	//         reversed the  order of DR and QINST back to what they were 25/9/98
	//   DGT had problems with qinst flows not preserving mass balance so resorted
	//   to routing the time step averages DR - simply by interchanging arguments
	//   because the DR passed back to CALCTS was not used for anything.

	//     2            r,PET,TIMES,NDATA,MPS,MPE,NGAUGE,DR,QINST,MP1)
	//   Understanding of Arguments DGT March 1998
	//   si - Array of subbasin initil conditions - array with column for each subbasin
	//       3 entries in each column SRO, ZBAR0, CV0
	//   Sp - Array of subbasin parameters - array with column for each subbasin
	//       11 entries
	//        1 - area - Subbasin area
	//        2 - szf   - exponential decay parameter
	//        3 - AKO  - vertical surface conductivity
	//        4 - DTH1 - difference between saturation and field capacity
	//        5 - DTH2 - difference between field capacity and wilting point
	//        6 - soilc - This is the maximum storage capacity (m) of the root zone.
	//            P-E enters the root zone and only exits if soilc is satisfied
	//        7 - C  Brooks Corey / Clapp Hornberger parameter
	//        8 - PSIf - Green - Ampt wetting front suction parameter
	//        9 - chv - Channel velocity
	//        10 - CC  Interception capacity
	//        11 - CR  Drying et enhancement ratio for drying rate of intercepted water.
	//
	//    isub - active subbasin number
	//    NKA - Array of number of ln(a/tan b) classes for each subbasin
	//    Lambda - Lambda for active subbasin
	//    ATB - two dimensional array with columns the ln(a/tan b) distribution for
	//          each subbasin
	//    PKA - two dimensional array with columns the proportion of basin corresponding
	//          to ln(a/tan b) values
	//    nd - array of number of flow distances for each subbasin
	//    CL - two dimensional array with columns giving flow distances for each subbasin
	//    PD - two dimensional array with columns giving cumulative flow distance
	//         probabilities for each subbasin
	//    units - a units conversion parameter
	//    IRR - A 2 dimensional output array of various diagnostic quantities
	//    modwrt - logical variable to control I/O
	//    IPSUB - Subbasin that output will be for - no longer used
	//    IPATB - not used
	//    STIM - Base time for Tideda output
	//    r - rainfall data - row for each raingage
	//    PET - potential et array
	//    interval - time step associated with each input value (seconds).  Assumed all equal -
	//        only first one is used.
	//    ndata - number of time steps
	//    ngauge - raingauge number being used for this subbasin
	//    dr - array that is direct runoff returned for this subbasin
	//    qinst - time series of outflow rates returned
	//    mp1 -not used

	static double sumi, sume;
	static string text;
	static double bal;
	static int tim;

	// local topmodel variables
	static int it, i, ind;
	static double sum, szf, soilc, dth1, dth2, ak0;
	static double te, chv, tdh1, tdh2, ae, rit, qv;
	static double sr0, zbar0;
	static double rof, roff, qlat, s1;
	static double rofex, rofrex, acie;
	static double dt, temp, zbarn;
	static double acse, qlatin;

	static ArrayXd srz(maxA), suz(maxA);

	const double i4max = 2147000000.0;
	// rpi 14/7/03 - turn scalars into arrays for removal of time loop
	// rpi 14/7/03 - the array dimension is maxslp
	// real*8 tdh(max_ntdh,maxslp)

	// landcare
	static int ymd, hms;
	static Array<int,Dynamic,1> surfrunoff(maxA);
	static ArrayXd soilwetness(maxA);
	//
	static int idebug, jt;
	static int nmin, jtdh, itdh;
	static int isr0msg=0, icv0msg=0;
	static double area, c, psif, cc, cr, r1, ud, r2;
	static double cv0;
	static double r3, rie, rd, srn, atbsat, atbinf, riei, rdi, rsei, ets, r3i, asatl;
	static double ainfl, amidl;
	static double puninf, pinf, psat, zi, sri, r3itemp, rdtemp, srni, balcv;
	static double balsat, balsoil;
	static double fi, crop_coeff, fsprinkler; //impervious fraction, crop coefficient is kc

	//  Artificial drainage variables
	static double adi, ad_cap, ad_now, ad_now_i;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> topmod(" << ncalls << ")" << std::endl;
    }
    caller = "topmod";
#endif

	// ***********************************************************************
	// Calculate the number of items to be put into the model's tideda file
	//      OPEN(UNIT=30,FILE='testclpd.txt',STATUS='unknown')
	idebug = 0; // put this > 0 to write water balance data to toperror.txt
	//   unravelling basin properties
	area       = Sp(1-1,isub-1);
	szf        = Sp(2-1,isub-1);
	ak0        = Sp(3-1,isub-1);
	dth1       = Sp(4-1,isub-1);
	dth2       = Sp(5-1,isub-1);
	soilc      = Sp(6-1,isub-1)*(dth1 + dth2);
	c          = Sp(7-1,isub-1);
	psif       = Sp(8-1,isub-1);
	chv        = Sp(9-1,isub-1);
	cc         = Sp(10-1,isub-1);
	cr         = Sp(11-1,isub-1);   //  Ratio by which evaporation from interception
	fi         = Sp(15-1,isub-1); // raw 12-jan-2005 fraction impervious
	fsprinkler = Sp(21-1,isub-1);
	crop_coeff = Sp(25+imonth-1,isub-1);
	ccompr2::t0 = Sp(25+12,isub-1);
	//    Derived parameters
	//      T0 = AK0/szf RAW special for nooksack 21-March-2005
	te = ccompr2::t0;   //  No spatial variability in Hyd conductivity yet

	// The following steps should only be executed at the start of a data sequence
	if (istep == 0) {
		//    is enhanced over the reference potential rate which is for unlimited
		//    transpiration from the vegetation
		dt  = interval/3600.0;  //  Time step in hours - Assumed equal interval time steps
		tim = stim;
		q0  = dt*te*exp(-Lambda);
		dth = dth1 + dth2;
		// No storage! need to prevent divide by zero
		if (soilc <= 0.0 && dth <= 0.0) {
			dth = 1.0;
		}
		zr = soilc/dth;
		ak0fzrdt = dt* ak0*exp(-szf*zr); // parameter used in infiltration excess
		//   Parmeters for analytic integration of sat zone store - See notebook analysis
		//   DGT April 1,1998
		logoqm = -log(szf/dth1) - log(1.0e-9*dt);
		//   The 1d-9 above is a m/hr tolerance on numerical estimation of changes in zbar
		qvmin = 1.0e-7*dth1/szf;

		//   The 1d-7 is a tolerance relating the analytic solution with qv not equal to
		//    zero to the solution with qv = 0.
		// Calculate Time delay histogram given chv (Distance/Time Step)
		// By always calculating this each time we enter this subroutine we
		// save a lot of space but at the price of time.
		ntdh[isub-1] = (int)(cl[nd[isub-1]-1][isub-1]/(chv*dt))+1;
		if (ntdh[isub-1] == 1) {
			tdh[0][isub-1] = 1.0;
		} else {
			if (ntdh[isub-1] > MAX_NTDH)
				ntdh[isub-1] = MAX_NTDH;
			tdh1 = 0.0;
			sum = 0.0;
			for (itdh = 1; itdh <= ntdh[isub-1]-1; itdh++) {
				for (ind = 2; ind <= nd[isub-1]; ind++) {
					if (chv*dt*itdh <= cl[ind-1][isub-1]) {
						tdh2 = pd[ind-2][isub-1] + (pd[ind-1][isub-1] - pd[ind-2][isub-1])*(itdh*chv*dt - cl[ind-2][isub-1]) / (cl[ind-1][isub-1] - cl[ind-2][isub-1]);
						tdh[itdh-1][isub-1] = tdh2 - tdh1;
						sum += tdh[itdh-1][isub-1];
						tdh1 = tdh2;
						break;	// continue itdh loop
					}
				}
			}
			tdh[ntdh[isub-1]-1][isub-1] = 0.0;
			for (itdh = 1; itdh <= ntdh[isub-1]; itdh++) {
				tdh[itdh-1][isub-1] = tdh[itdh-1][isub-1]/sum;
			}
			// rpi 960306 removed this to prevent all routing being dumped into the
			// last tdh position by default. Instead spread residual error
			// proportionately across all values
			//      TDH(ntdh)=1-sum
		}
		//    Unravelling initial conditions
		sr0 = si[0][isub-1];
		if (sr0 > soilc) {
			if (isr0msg <= 10) {
				cerr << "Initial Sr0 " << fixed << setw(12) << setprecision(5) << sr0;
				cerr << " too big, reduced to soilc " << scientific << setw(12) << setprecision(6) << soilc << '\n';
			}
			isr0msg++;
			sr0 = soilc;
		}
		zbar0 = si[1][isub-1];
		cv0 = si[2][isub-1];
		if (cv0 > cc) {
			if (icv0msg <= 10) {
				cerr << "Initial cv0 " << fixed << setw(8) << setprecision(5) << cv0;
				cerr << " too big, reduced to cc " << scientific << setw(12) << setprecision(6) << cc << '\n';
			}
			icv0msg++;
			cv0 = cc;
		}
		irr  = 0;		// array assignment
	    rirr = 0.0;		// array assignment
		zbar = zbar0;
		sr   = sr0;
		cv   = cv0;
		//  initial storage for mass balance check
		s0 = -(zbar)*dth1 + sr + cv;

		// reset maximum contributing area variables for new storm period
		// rpi reintroduced the next  1 lines for infil. excess mechanism 13/1/94
		//      write(30,9955) ntdh,cl,pd,tdh
		// 9955 format(1x,'ntdh,cl,pd,tdh',i4/(1x,5g14.3))
		aciem  = 0.0;    // max area contributing by infiltration excess
		acsem  = 0.0;    // max area contributing by saturation excess
		zbm    = zbar;    // max zbar
		sumr   = 0.0;    // sum of rainfall
		sumae  = 0.0;   // sum of actual evapotranspiration
		sumpe  = 0.0;   // sum of potential evapotranspiration
		sumq   = 0.0;    // sum of all basin outflow - before routing
		sumie  = 0.0;   // sum of infiltration excess runoff - before routing
		sumse  = 0.0;   // sum of saturation excess runoff - before routing
		sumqb  = 0.0;   // sum of baseflow - before routing
		// sumqm = 0.0;   ! sum of measured flows
		sumce  = 0.0;   // sum of canopy actual e/t
		sumsle = 0.0;  // sum of soil actual e/t (ie when r3<0)
		sumr1  = 0.0;   // sum of throughfall
		// sumr2 = 0.0;   ! sum of net forcing at soil surface
		sumqv  = 0.0;   // sum of soil drainage to sat zone
		// suminf = 0.0;  !sum of infilt at soil surface (ie when r3>0)
	    temp   = -szf*zbar;
	    if ( temp < -100.0 ) {
			qb = 0.0;
	    } else {
			qb = q0*exp(temp);
	    }
		//  initialise output array
		qinst = 0.0;	// Array assignment
		dr = 0.0;		// Array assignment
      	qinst(0) = qb/units*area/interval;  // dgt converts m/ts to mm^3/sec
		qinst_out = qinst(0); // rpi 16/7/2003 - needed to get at qinst
		dr_out = dr(0);   // dgt 6/10/05 - to capture dr out rather than use irr which suffers from rounding
		// and does not occur if printing is off
		for (it = 1; it <= min(ntdh[isub-1]-1, ndata); it++) {   // dgt 3/16/98  just in case ndth > ndata
			// added by rpi to transfer out initial baseflow for indepth to use
			temp = 1.0;
			for (jt = 1; jt <= it; jt++) {	// changed to go to it so as not to double count first step
				temp -= tdh[jt-1][isub-1];
			}
        	temp = max(0.0, temp);  // dgt 5/26/98   guard against small negative
			qinst(it) = qinst(0)*temp;
			// dr(1,1) gets over-written below so value has to go into it+1 (fortran indexing, it in C indexing)
			dr(it) = qinst(0)*temp;     // mass conservative time step integrated flow mm^3/ts
		}
		goto L4992;
	} // time0


	// ***********************************************************************

	//  Start loop on time steps
	//      DO 50 IT=1,NDATA


	//  Convert inputs to m per time step
	it = istep;
	ae = pet*units*crop_coeff;       // 12-jan-2005 crop_coeff is crop factor (kc)
	sumpe += ae;
	rit = (r + rate_irrig)*units;   //(isub)
	ad_cap = art_drainage*dt;   // dgt 6/28/05   artificial drainage capacity in units per time step (m/day)
	qv = 0.0;
	sumi = 0.0;  // rpi 22/3/2004 removed the subscript as unnecessary
	//	srav=0.0;
	//
	//   interception layer   dgt 3/98
	// 12-jan-2005 of the irrig water, a fraction 1-fsprinkler bypasses the canopy and reaches the ground directly
	//so don't apply all the water to the canopy top. add some of it to the thrufall instead
	intercept(cv, rit - (1.0 - fsprinkler)*rate_irrig*units, ae, dt, cc, cr, r1, ud);  // dgt 8/18/05  multiplied rate_irrig by units
	r1 += (1.0 - fsprinkler)*rate_irrig*units;
	sumce += ae;
	sumr1 += r1;
	//	write(21,*)'it,balcv,rit,cv,cvs,r1,ae', it,balcv,rit,cv,cvs,r1,ae
	//  cv is interception state variable [m] (input and returned)
	//  RIT is precipitation amount [m] (input)
	//  AE is reference evaporation amount [m] on input, and upon return is net evap
	//     from interception
	//  dt is time step [hr] (input)
	//  CC is interception capacity parameter [m] (input)
	//  CR is interception evaporation enhancement factor parameter (input)
	//  r1 net precipitation (throughfall + stemflow)
	//  ud is unsatisfied transpiration demand that goes through to soil.
    r2 = r1 - ud;   // Net forcing on the soil zone.
	//      sumr2=sumr2+r2
    // thinks some ET is missing from the above
	// since r2=r1-ud, then the min(r1,ud) was already removed before soil
	sumsle += min(r1, ud);
	sume = ae + min(r1, ud); // removed the subscript as unnecessary
	//   basinwide soil zone calculations
	soil(sr, r2, r3, rie, rd, srn, zr, dth1, dth2, szf, ak0fzrdt, c, soilc, psif, true, ad_cap, ad_now);
	//   calculate saturation thresholds
	atbsat = szf*zbar + Lambda;   //  Saturation threshold
	//       Points with ln(a/tan b) greater than this are saturated
	atbinf = atbsat - szf*zr;   //  points with ln(a/tan b) greater than this
	//       are 'influenced' by shallow water table

	nmin = 2;

	//    Start loop on a/TtanB increments at 2 because each increment represents a bin from i-1 to i
	//    This code assumes ln(a/tan b) values are in increasing order and NOT repeated
	acie = 0.0;      //   infiltration excess area
	acse = 0.0;      //   saturation excess area
	rofex = 0.0;     //   runoff saturation excess
	rofrex = 0.0;    //   runoff infiltration excess
	// raw 12-jan-2005 impervious area
	rof = max(r2, 0.0)*fi;
	rofrex = rofrex+rof;
	sumie += rof;
	acie  += fi;
	sumad = 0.0;   //  dgt 6/28/05 this is a variable to accumulate artificial drainage across classes

	for (i = nmin; i <= Nka[isub-1]; i++) {
		//       Check for saturation
        if (atb[i-2][isub-1]  >=  atbsat) {
			//          Saturated
			riei = 0.0;   //  No infiltration excess runoff
			if (r2 <= 0.0) {   //  Potential ET greater than precip
				//  To preserve mass balance since the soil is saturated any change in sr needs
				//  to be added to saturated zone store.  Do this as a component (which may be negative)
				//  of the soil zone drainage.
				rdi = r2 + sr - srn;   // ET demand taken from sat zone
				rsei = 0.0;   //  No runoff
				ets = -r2;   //  ET from soil satisfies all ET demand
				r3i = r2;  // RPI 3/5/01
			} else {
				rsei = r2;
				//		  rdi=sr-srn   ! CX
				// RPI and RAW 24/4/01 prevent -ve seepage in this case, i.e., there is
				// rainfall but no infiltration because saturated conditions
				rdi = max(0.0, sr - srn);   // For mass balance change in soil moisture taken from sat zone

				r3i = 0.0;   // No infiltration
				ets = 0.0;
				//      no actual evap this case because it is already taken care of in interception
			}
			//  DGT 6/28/05 Artificial drainage from saturated soil store is at capacity
			adi = ad_cap;
			rdi = rdi - ad_cap;   //  Keep mass balance right
			//      Accumulate saturated area
			acse += (1.0 - fi)*pka[i-1][isub-1];
		} else if (atb[i-1][isub-1] <= atbinf) {
			//        Not influenced by water table at all
			//        Use deep water table values
			r3i = r3;
			riei = rie;
			rdi = rd;
			rsei = 0.0;
			ets = 0.0;
			if (r3i < 0.0)
				ets = -r3i;   //  ET
			//  DGT 6/28/05 Artificial drainage at uninfluenced soil value
			adi = ad_now;   // No need to adjust rd because soil did it
		} else {

			//      Here part of the zone is either saturated or influenced by the water table

			//        compute ln(a/tan b) limits
			asatl = min(atbsat, atb[i-1][isub-1]);
			ainfl = max(atbinf, atb[i-2][isub-1]);
			amidl = (asatl + ainfl)*0.5;   //  Use midpoint a/tan b as representative
			//        compute proportions in each category
			temp = 1.0/(atb[i-1][isub-1] - atb[i-2][isub-1]);
			puninf = (ainfl - atb[i-2][isub-1])*temp;
			pinf = (asatl - ainfl)*temp;
			psat = (atb[i-1][isub-1] - asatl)*temp;

			//       Saturated part

			if (r2 <= 0.0) {   //  Potential ET greater than precip
				//  To preserve mass balance since the soil is saturated any change in sr needs
				//  to be added to saturated zone store.  Do this as a component (which may be negative)
				//  of the soil zone drainage.
				rdi = psat *(r2 + sr - srn);   // ET demand in excess of soil zone demand taken from sat zone
				rsei = 0.0;   //  No runoff
				ets = -psat*r2;   //  ET from soil satisfies all ET demand
				r3i = r2*psat;
			} else {
				rsei = psat*r2;
				// RPI and RAW 24/4/01 prevent -ve seepage in this case, i.e., there is
				// rainfall but no infiltration because saturated conditions
				rdi = max(0.0, psat*(sr - srn));
				//	      rdi = psat*(sr-srn)  !CX
				ets = 0.0;
				r3i = 0.0;
				//      no actual evap this case because it is already taken care of in interception
			}
			acse += psat*(1.0 - fi)*pka[i-1][isub-1];
			//  DGT 6/28/05   Artificial drainage
			adi = ad_cap*psat;
			rdi -= ad_cap*psat;

			//    Influenced part
			//    Calculate influenced category zi
			zi = zbar + (Lambda - amidl)/szf;
			sri = sr + (zr - zi)*(1.0 - sr/soilc)*dth; // local enhancement of soil moisture
			soil(sri, r2, r3itemp, riei, rdtemp, srni, zr, dth1, dth2, szf, ak0fzrdt, c, soilc, psif, false, ad_cap, ad_now_i);  // RPI 6/5/01 introduced r3ie
			//     6/28/05   DGT introduced last two arguments for artificial drainage

			//          The main use of the above call is to get
			//           local infiltration or ET.  This is r3i, which is then used
			//           with the global storage change to infer local sat zone addition
			//           or withdrawal.  The concept of flux to a sat zone is suspended
			//           here because the sat zone and soil zone have overlapped.  This
			//           is for mass balance.
			//           The rd returned from soil is not used.  Instead mass balance
			//           is used to infer the local sat zone storage reallocation.
			// RPI and RAW 24/4/01 prevent -ve seepage in this case, i.e., there is
			// rainfall but no infiltration because saturated conditions
			rdtemp = pinf*(r3itemp + sr - srn);
			if (r2 > 0.0)
				rdtemp = max(0.0, rdtemp);
			//          rdi= rdi + pinf*(r3i+sr-srn)
			rdi += rdtemp;
			if (r3itemp < 0.0)
				ets -= pinf* r3itemp;
			r3i += r3itemp*pinf;  // RPI 3/5/01 get weighted contribution
			riei = pinf * riei;
			//  DGT 6/28/05 artificial drainage
			adi = ad_now_i*pinf;

			//    Uninfluenced part
			//        Use deep water table values
			riei += puninf*rie;
			rdi += puninf*rd;
			if (r3 < 0.0)
				ets = ets -puninf*r3;   //  ET
			r3i += r3*puninf;  // RPI 3/5/01 get weighted contribution
		}
		//   summing contributions
		qv     += rdi*(1.0 - fi)*pka[i-1][isub-1];
		sumi   += r3i*(1-fi)*pka[i-1][isub-1]; // rpi 22/3/2004 removed the subscript as unnecessary
		sume   += ets*(1-fi)*pka[i-1][isub-1];  // rpi 22/3/2004 removed the subscript as unnecessary
		sumsle += ets*(1.0 - fi)*pka[i-1][isub-1];
		if (rsei >= 0.0) {
			rof = rsei*(1.0 - fi)*pka[i-1][isub-1];
			rofex += rof;
			sumse += rof;
		} else {
			cerr << "saturation excess less than 0 error, rsei&psat=" << rsei << " " << psat << '\n';
			exit(EXIT_FAILURE);
		}

			//  case of infiltration excess
		if (riei >= 0.0) {
			rof = riei*(1-fi)*pka[i-1][isub-1];
			rofrex += rof;
			sumie  += rof;
			// rpi reintroduced the next  1 lines for infil. excess mechanism 13/1/94
			if (riei > 0.0)
				acie += acie  + (1.0 - fi)*pka[i-1][isub-1];
		} else {
			cerr << "infiltration excess less than 0 error\n";
		}
		//  dgt 6/28/05   artificial drainage
		sumad += adi*(1.0-fi)*pka[i-1][isub-1];

		if ( (modwrt) && ((it >= mps) && (it <= mpe))) {
			// 17/4/98 we don't have z(i) anymore, so just calculate at midpoint
			lundatFile << " " << dec << setw(4) << isub;
			lundatFile << dec << setw(4) << it;
			lundatFile << dec << setw(4) << i;
			lundatFile << fixed << setw(13) << setprecision(3) << atb[i-1][isub-1];
			lundatFile << fixed << setw(13) << setprecision(3) << zbar + (Lambda-(atb[i-2][isub-1] + atb[i-1][isub-1])*0.5)/szf;
			lundatFile << fixed << setw(13) << setprecision(3) << rsei;
			lundatFile << fixed << setw(13) << setprecision(3) << riei << '\n';
		}
		// Landcare, surfrunoff in um/timestep (rsei is in m),
		// soilwetness is saturated depth as a fraction of soil depth = (ZR-ZI)/ZR = 1-ZI/ZR
		surfrunoff(i-1) = min((rsei + riei)*1.0e6, i4max) + 0.5;
		soilwetness(i-1) = min(1.0, max(0.0, 1.0 - (zbar + (Lambda - (atb[i-2][isub-1] + atb[i-1][isub-1])*0.5)/szf)/zr));

		//      if (isub  ==  87 && ipflag  ==  1 && istep  ==  245)then
		//	  write(3284,*)i,r3i,ets
		//	endif
	}
	if (modwrt && mpe == -1) {
		if (it == 1) {  //write a header for this subcatchment
			lundatFile << isub << " " << ndata;
			for (i = 0; i < maxSlp; i++) {
				lundatFile << " " << Nka[i];
			}
			lundatFile << '\n';
		}
		td81micdh(ymd, hms, tim);
		lundatFile << dec << setw(8) << ymd;
		for (i = 0; i < Nka[isub-1]; i++) {
			lundatFile << dec << setw(8) << surfrunoff[i];
		}
		for (i = 0; i < Nka[isub-1]; i++) {
			lundatFile << fixed << setw(6) << setprecision(2) << soilwetness[i];
		}
		if (isub == 1)
			tim += interval;
	}
	// RPI and RAW 24/4/01 have prevented -ve seepage so now
	// we must calculate the correct soil moisture value to conserve mass
	srn = sr + sumi - qv - sumad; //  DGT 6/28/05 added artificial drainage effect

	//      if (ipflag  ==  1)then
	//	   write(3284,'(1x,i5,7e14.5)')istep,r2,sumi,sume,sr,srn,qv,sumad
	//	endif
	//	srn=sr+sumi-qv ! RPI 22/3/2004 removed the subscript as unnecessary

	//    update soil zone store for next iteration

	sr = srn;
	roff = rofex + rofrex + sumad;   //  dgt 6/28/05  adding artificial drainage to runoff
	//	roff=rofex+rofrex   //  total runoff
	//      if ( modwrt && (isub  ==  ipsub) ) then
	if ( modwrt ) {
		  irr(6) = min( roff*1.0e6, i4max ) + 0.5; //um/timestep
		  rirr(6) = roff*1000.0;   // units for rirr are mm wherever possible
	}

	//   What follows is an analytic integral of
	//     d(DTH1*ZBAR)/dt = -QV + Q0 exp(-szf * Zbar)
	//    over a unit time step (because QV and Q0 are both per time step
	//    Then base flow qb is obtained from the implied storage difference
	//    In this derivation zbar can be negative - therefore it makes no
	//    physical sense to try interpret it in terms of the soil zone issues.
	//    dth1 is simply a multiplying factor that scales saturated zone storage.

    //  QB=Q0*DEXP(-szf*ZBAR)
    //  ZBAR=ZBAR-(QV-QB)/DTH1

	//      TEMP=DEXP( MIN(szf*ZBAR, 100D0 ) )
	//      if ( QV > Q0*1D-6 ) THEN
	//	    TEMP1=DEXP( DMAX1(-QV*szf/DTH1, -100D0 ) )
	//	    ZBARN = LOG( Q0/QV + TEMP1*(TEMP-Q0/QV) ) / szf
	//      ELSE
	//	    ZBARN = LOG( TEMP + Q0*szf/DTH1 ) / szf
	//      ENDIF
	//   DGT revised following careful analysis in notebook 4-1-98
#ifdef ZBAR_IN
    zbarn = zbar_new;
#else
	if (max(zbar*szf, zbar*szf-qv/dth1) > logoqm) {
	//  Forward difference
		zbarn = zbar + (-qv + q0*exp(-szf*zbar))/dth1;
	} else {
		if (fabs(qv) > qvmin) {
			//   Use non zero qv solution
			zbarn = log(q0/qv + exp(-qv*szf/dth1)*(exp(szf*zbar) - q0/qv))/szf;
		} else {
			//    Use zero qv solution
			zbarn = log(exp(szf*zbar) + q0*szf/dth1)/szf;
		}
	}
#endif
	sumqv += qv;
	qb = qv + (zbarn - zbar)*dth1;
	//       write(51,592) it,isub,zbar,zbarn,dth1,qb,qv
	// 592   format(1x,2i3,5f12.6)
	zbar = zbarn;

	sumq += qb + roff;
	sumr += rit;
	sumqb += qb;
	qlat_mm = (qb + roff)/units;
	//  dgt changed below to get units right
	//  convert runoff back to input data units
	qlat = (qb + roff)/units*area/interval;   // mm^3/sec
	// calculate instantaneous flow for reach routing - rpi
	//      qlatin=(q0*dexp(-szf*zbarn)+roff)/units*area/interval   ! mm^3/sec
	qlatin = (qb)/units*area/interval;   // mm^3/sec
	sumad = sumad/units*area/interval;   // dgt 6/28/05  - put in vol units

	//  use time delay histogram to route runoff (assumed known)
	// rpi 31/7/2002 for neq=0 purposes we have to let this convolution go
	// beyond m values to m+ntdh values
	for (itdh = 1; itdh <= ntdh[isub-1]; itdh++) {
		dr(itdh-1) = dr(itdh-1) + qlat*tdh[itdh-1][isub-1];
		qinst(itdh-1) += qlatin*tdh[itdh-1][isub-1];
	}
	qinst_out = qinst(0); // rpi 16/7/2003 - needed to get at qinst
	dr_out = dr(0);  // dgt 6/10/05 to get at dr

	sumae += sume;

	// rpi reintroduced the next  1 lines for infil. excess mechanism 13/1/94
	if (acie > aciem)
		aciem = acie;
	if (acse > acsem)
		acsem = acse;
	if (zbm > zbar)
		zbm = zbar;

	//  End of time step loop


	//   DGT 8/17/05  Add basin_evap so as not to have to use irr for this
      evap_mm = sume*1.0e3;   // This is in mm
	//      if ( MODWRT && (isub  ==  IPSUB)) THEN
	if ( modwrt ) {
		if (isr0msg > 10 && (it == 1 && isub == 1))
			cerr << " ***** message about sr0 exceeding soilc exceeded a total of " << dec << setw(6) << isr0msg << " times\n";
		if (icv0msg > 10 && (it == 1 && isub == 1))
			cerr << " ***** message about cv0 exceeding soilc exceeded a total of " << dec << setw(6) << icv0msg << " times\n";
		irr(1-1)  = it;
		// fluxes to go into irr as um/hr
		// need 1d3/area for irr(2,) since dr is in mm^3/int, not m/int like qv
		// dr is now in mm^3/sec therefore need interval which is time step in seconds
		irr(1-1)  = min( dr(0)*1.0e3/(area)*interval, i4max ) + 0.5;  // rpi 16/7/2003 changed it to 1 in dr() um/timestep
	    rirr(2-1) = dr(0)*interval/area;   // in mm
		// added the following 3 lines - see above
	    irr(3-1)   = min( qb*1.0e6, i4max ) + 0.5; //um/timestep
		rirr(3-1)  = qb*1000.0;  // mm
	    irr(4-1)   = min( qv*1.0e6, i4max) + 0.5; //um/timestep
		rirr(4-1)  = qv*1000.0; // mm
	    irr(5-1)   = min( rofex*1.0e6, i4max ) + 0.5; //um/timestep
	    rirr(5-1)  = rofex*1000.0; // mm
	    irr(6-1)   = min( rofrex*1.0e6, i4max ) + 0.5; //um/timestep
		rirr(6-1)  = rofrex*1000.0; //mm
	    irr(8-1)   = acie*100.0; //percentage
		rirr(8-1)  = acie*100.0;
	    irr(9-1)   = acse*100.0; //percentage
		rirr(9-1)  = acse*100.0;
	    irr(11-1)  = min( zbar*1.0e6, i4max) + 0.5; //micrometre
		rirr(11-1) = zbar*1000.0;  //mm
		// revised next two lines, since cv and sr replace srz and suz 17/4/98
	    irr(10-1)  = min( cv*1.0e6, i4max) + 0.5; //micrometre
		rirr(10-1) = cv*1000.0;  //mm
	    irr(12-1)  = min( sr*1.0e6, i4max) + 0.5; //micrometre
		rirr(12-1) = sr*1000.0; //mm
		irr(13-1)  = min(pet*units*1.0e6, i4max) + 0.5;   // dgt added pet !um/timestep
		rirr(13-1) = pet;   // mm
		irr(14-1)  = min(sume*1.0e6, i4max)+0.5;   // dgt added actual et ! rpi 22/3/2004 removed the subscript as unnecessary
		rirr(14-1) = sume*1000.0; // mm
		s1 = -zbar*dth1 + sr + cv;
		bal = sumr - sumq - sumae + (s0 - s1);
	    irr(15-1) = min(bal*1.0e9, i4max) + 0.5; //dgt mass balance discrepancy
		//            in units of m x e-9, i.e. 1000 ths of mm.
		rirr(15-1)=bal*1000.0; //mm

		if ( idebug > 0 ) {
			cerr << "it=" << it << "water balance: bal=" << bal << " sumr=" << sumr << '\n';
			cerr << "isub,bal,sumr , s1, s0, sumq , sumae \n";
			cerr << isub << " " << bal << " " << sumr << " " << s1 << " " << s0 << " " << sumq << " " << sumae << '\n';
			balcv = sumr - ((cv - cv0) + sumr1 + sumce);
			cerr << "isub,balcv,sumr,cv,cv0,sumr1,sumce \n";
			cerr << isub << " " << balcv << " " << sumr << " " << cv << " " << cv0 << " " << sumr1 << " " << sumce << '\n';
			balsoil = sumr1 - ((sr - sr0) + sumqv + sumse + sumie + sumsle);
			cerr << "isub,balsoil,sumr1,sr,sr0,sumsle \n";
			cerr << isub << " " << balsoil << " " << sumr1 << " " << sr << " " << sr0 << " " << sumsle << '\n';
			balsat = sumqv - (-(zbar - zbar0)*dth1 + sumqb);
			cerr << "balsat,sumqv,-(zbar-zbar0)*dth1,sumqb \n";
			cerr << balsat << " " << sumqv << " " << -(zbar-zbar0)*dth1 << " " << sumqb << '\n';
			cerr << '\n';
			if ( fabs(bal) > (sumr+sumq)*1.0e-3 ) {
				// can't write to lunco1, so put it in toperror.txt
				cerr << "water balance: bal=" << bal << " sumr=" << sumr << " sumq=" << sumq << " it=" << it << " isub=" << isub << '\n';
			}
		}
	}
	// RPI 16/7/2003 put in this test for the re-arranged version of TOPNET
	// These stmts "push down" the contents of the DR and QINST arrays so that
	// they don't grow beyond the size of the Time-Delay-Histogram.

L4992:	if (ntdh[isub-1] > 1) {
		for (jtdh = 1; jtdh <= ntdh[isub-1]-1; jtdh++) {
			dr(jtdh-1) = dr(jtdh);
			qinst(jtdh-1) = qinst(jtdh);
		}
		dr(ntdh[isub-1]-1)    = 0.0;
		qinst(ntdh[isub-1]-1) = 0.0;
	} else {
			dr(0)    = 0.0;
			qinst(0) = 0.0;
	}
	//-------------------------------------------------------------
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving topmod(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif
	return 0;
}


// **********************************************************************
//   DGT's subroutine for infiltration.
// **********************************************************************
int soil(const double sr, const double r2, double &r3, double &rie, double &rd, double &srn,
	const double zr, const double dth1, const double dth2, const double szf,
	const double ak0fzrdt, const double c, const double soilc,
	const double psif, const bool global, const double ad_cap, double &ad_now)
{
	//  Inputs are
	//    sr  soil storage state variable
	//    r2  surface forcing (precip-ET)
	//    ad_cap   artificial drainage capacity
	//  Outputs are
	//    r3  surface exchange that occurs limited by soil moisture and infiltration capacity
	//    rie infiltration excess runoff
	//    rd  recharge to saturated zone
	//    srn  updated soil storage
	//    ad_now    artificial drainage that occurs limited by availability above field capacity
	//  Other arguments are parameters

	double temp;
	double srd;
	double rdmax, f, zf, sre;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> soil(" << ncalls << ")" << std::endl;
    }
    caller = "soil";
#endif
	temp = sr - zr*dth2;
	//	write(21,*)temp,sr
	if (temp > 0.0) {
		srd = temp/(zr*dth1);
		rd = ak0fzrdt*pow(srd, c); //  + art_drainage     // 12-jan-2005 DGT 6/28/05  Drainage does not go to sat store
		//   Do not drain more than is in drainable capacity + input
		rdmax = temp + max(r2, 0.0);
		rd = min(rdmax, rd);
		ad_now = min(ad_cap, rdmax - rd);  // Drainage is limited by capacity or supply
	} else {
		//	  srd=0.0;   ! not used below
		rd = 0.0;
		ad_now = 0.0;
	}
	if (r2 > 0.0) {
		f = soilc - sr + rd + ad_now;
		zf = sr/(dth1 + dth2);
		if (zf > 0.0)
			f = min(ak0fzrdt*exp(szf*(zr-zf))*(zf + psif)/zf, f);  //  max depth that can be added to root zone
		r3 = min(f, r2);   //  actual depth added
		rie = r2 - r3;    // infiltration excess runoff generated
	} else {
		sre = min(1.0, sr/(zr*dth2));  //  relative plant available saturation
		r3 = sre*r2;
		rie = 0.0;
		if (global && r3 < -sr)r3 = -sr;  //  smallest negative number.
		//      cannot evap more than is there in global case
		//      In local case, i.e. where water table is near to the surface we assume
		//      that due to lateral interflow that sr is not limiting ET.
	}
	srn = sr + r3 - rd - ad_now;
	//   checks to correct for overshooting.  Reduce ad_now first
	if (srn < 0.0) {
		if(-srn > ad_now) {
			ad_now = 0.0;
		} else {
			ad_now += srn;
		}
	}
	srn = sr + r3 - rd - ad_now;
	if (srn < 0) {   // may go negative due to big rd.
		rd += srn;   //  reduce rd accordingly
		srn = 0.0;
	}
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving soil(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif
	return 0;
}


// ***********************************************************************
//     SUBROUTINE INTERCEPT
// ************************************************************************
int intercept(double &cv, const double rit, double &ae, const double delt, const double x1, const double cr, double &r1, double &ud)
{
    double ii;
    double s, dt, p, x2;
	double e, a, b, tt, sn, r;
	//  Subroutine to implement nonlinear interception model

	//  cv is interception state variable [m] (input and returned)
	//  RIT is precipitation amount [m] (input)
	//  AE is reference evaporation amount [m] on input, and upon return is evap
	//     from interception accumulated over time step
	//  dt is time step [hr] (input)
	//  CC = x1 is interception capacity parameter [m] (input)
	//  CR is interception evaporation enhancement factor parameter (input)
	//  r1 net precipitation (throughfall + stemflow)
	//  ud is unsatisfied transpiration demand that goes through to soil.

	//  For the theory see TopNetInterceptionDerivation.docx last saved in
	//   C:\Users\dtarb\Dave\Projects\WRIA1_Water_Budget_Christina\TopNetModel\Documents
	//   DGT 5/28/12

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> intercept(" << ncalls << ")" << std::endl;
    }
    caller = "intercept";
#endif
	if (x1 <= 0.0) {   // case of no interception capacity
		cv = 0.0;
		r1 = rit;
		ud = ae;
		ae = 0.0;
		return 0;
	}
	//  Normalize
	s = cv/x1;
	dt = delt;
	p = rit/(x1*dt);
	// RPI 20/6/2002 put in the "check" on CR to avoid generation of NaNs
	x2 = x1/max(cr, 1.0e-8);
	e = ae/(x2*dt);
	if ((p + e)*dt > 10.0e-9) {
		dt = dt*(p + e);
		if (e > p*10.0e-8) {
			// General case
			// Compute roots
			a = 1.0 + sqrt(e/(p + e));
			b = 2.0 - a;
			//			b = 1 - sqrt(e/(p + e))
			//  The choice of a as larger root >1 is important for
			//  next step to prevent possibility
			//  of divide by zero error when s = a.
			tt = (s - b)/(a - s)*exp( -dt*(a - b));
			sn = (b + a * tt)/(1.0 + tt);
		} else {
			// Special case for e=0
			sn = ((1.0 - s)*dt + s)/((1.0 - s)*dt + 1.0);
		}
		dt = dt/(p + e);
		// Now unnormalize
		sn = sn * x1;
		s = s * x1;
		p = p * x1;
		e = e * x2;
		// Integral of f(S) by mass balance
		ii = (s + p*dt - sn)/(p + x1/x2 * e);
		r = p * ii;
		ae = x1/x2*e*ii;
		ud = e*(dt - ii);
	} else {
		// The case when p and e are both 0 is trivial
		//  no change in s
		sn = s*x1;
		r  = 0.0;
		ae = 0.0;
		ud = 0.0;
	}
	//   Values to return
	cv = sn;
	r1 = r;
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving intercept(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif
	return 0;
}

int irrigation(const double Sr, const ArrayXXd &Sp, const int js, const long int interval, const int maxSlp, double &Dapp)
{
	//REAL*8 Sp(NSP,MAXSLP)
	double dthres, dgoal, zmax, thetathresh;
	//double ieff;
	double zr, dth2, theta, thetagoal;

	//  DGT  8/16/05  Decided to try to reinterpret dthres and dgoal as threshold soil moisture fractions
	//  of the plant available water.  So when soil moisture content is less than dthresh*DTH2 irrigation starts
	//  Irrigation does not apply water that would make the soil wetter than dgoal*dth2.
#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> irrigation(" << ncalls << ")" << std::endl;
    }
    caller = "irrigation";
#endif
    zr      = Sp(5,js-1); //metres
	//	DTH1=Sp(4,js)
    dth2    = Sp(4,js-1);
	//	DTH=dth1+dth2

	dthres = Sp(22,js-1);  // 0.7   min(Sp(23,JS),soilc) metres   DGT 8/16/05 changed to fraction of plant available soil water
	dgoal  = Sp(24,js-1);   // 0.95  min(Sp(25,JS),soilc) metres  DGT 8/16/05 changed to fraction
	//ieff   = Sp[21][js-1]; //dimensionless
	zmax   = Sp(23,js-1); //metres/day
	//		 21 - SprinklerFractionOfIrrigation	(n.b. Fortran index)
	//		 22 - IrrigationEfficiency
	//		 23 - D_Thresh
	//		 24 - Z_Max
	//		 25 - D_Goal
	//		 26-37 - Kc_Jan through to Kc_Dec
	//	if((soilc-sr) > dthres) then
	//		Dapp=min( ((soilc-sr)-dgoal)/ieff,
	//	1				zmax*float(interval)/(3600.*24.)) !m/timestep = m/d * s/timestep / (s/day)
	//      if (sr < dthres*zr*dth2)then   ! DTH2 is the plant available water
	//	    Dapp=min((dgoal*zr*dth2-sr)/ieff,
	//     1                zmax*float(interval)/(3600.*24.))  !m/timestep = m/d * s/timestep / (s/day)
	//   DGT 8/19/05 Gradual approach more commensurate with irrigation of larger areas
	//   Demand 'Dapp' increases linearly from 0 at dgoal to zmax at dt
	theta       = Sr/zr;  // soil moisture content
	thetagoal   = dgoal*dth2;
	thetathresh = dthres*dth2;
	if (theta < thetagoal) {
	    Dapp = min((thetagoal - theta)/(max(thetagoal - thetathresh, 1e-6))*zmax, zmax*(double)(interval)/(3600.0*24.0));
		// the 1e-6 included to avoid divide by 0
		//m/timestep = m/d * s/timestep / (s/day)
	} else {
		Dapp = 0;
	}
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving irrigation(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif
	return 0;
}

