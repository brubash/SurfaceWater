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
#include <sstream>
#include "snow.hh"

using namespace std;
using namespace Eigen;

//  Subroutines and function subprograms for the Utah Energy Balance
//  Snow Accumulation and Melt Model.
//  David G. Tarboton, Utah Water Research Laboratory, Utah State University
//  Version 2.1
//  Last Change 5/4/04 to accommodate glacier surface runoff.

//  This program was written and prepared under agreement with and funding
//  by the U.S. Government , and therefore is in the public domain and not
//  subject to copyright.

// **************READ Bristow Campbell Parameter values **************

int bcparm(vector<double> &dtBar, double &bca, double &bcc, const string bcFile)
{
	istringstream line;
	string inLine;
	int month;
	double val;
	ifstream paramFile(bcFile.c_str());	// unit 88
	if (!paramFile.is_open()) {
		cerr << "Failed to open " << bcFile << '\n';
		exit(EXIT_FAILURE);
	}
	paramFile >> bca >> bcc;
	getline(paramFile, inLine, '\n');			// Read the remainder of the line
	getline(paramFile, inLine, '\n');			// Read the next line.
	while (inLine.length() > 1) {
		line.str(inLine);
		line >> month >> val;
		dtBar[month-1] = val;
		getline(paramFile, inLine, '\n');
	}
	paramFile.close();

	return 0;
}


// ************************** updatetime () ***************************
//                 update time for each time step
int lyear(const int year);

int updatetime(int &year, int &month, int &day, double &hour, const double dt)
{
	int dm;

	const int dmon[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	hour += dt;
	dm = dmon[month-1];
	// check for leap years
	if(month == 2)
		dm = lyear(year);
L10:;
	if (hour >= 24.0) {
		hour -= 24.0;
		day++;
		goto L10;
	}
L20:;
	if(day > dm) {
		day -= dm;
		month++;
		if(month > 12) {
			month = 1;
			year++;
			dm=dmon[month-1];
			if(month == 2)
				dm = lyear(year);
		}
		goto L20;
	}

	return 0;
}

// ************************** lyear () ***************************
//     function to return number of days in February checking for leap years
int lyear(const int year)
{
	if ((year % 4) > 0 || ((year % 100)  == 0 && (year % 400) != 0)) {
	//  Leap years are every 4 years
	//  - except for years that are multiples of centuries (e.g. 1800, 1900)
	//  - except again that when the century is divisible by 4 (e.g. 1600, 2000)
	//    then it is a leap year
        return 28;
	} else {
		return 29;
	}
	return 0;
}

// **************************** atf () ****************************
//     to get the atmospheric transmissivity using the Bristow and Campbell
//     (1984) approach

int atf(double &atff, const double trange, const int month, const vector<double> dtbar, const double a, const double c)
{
	double b;

	b = 0.036*exp(-0.154*dtbar[month-1]);
	atff = a*(1.0 - exp(-b*pow(trange, c)));
	//      write(6,*)trange,month,a,c,dtbar(month),atf
	return 0;
}


// ************************** hourlyRI () **********************
//                To get hourly radiation index
int julian(const int yy, const int mm, const int dd);

int hyri(const int year, const int month, const int day, const double hour, const double dt,
	const double slope, const double azi, const double lat, double &hri, double &coszen)
{
	double lp,lat1, crad, t, delt1, slope1, azi1, fjulian, d, tanprod;
	double td, ddt, tpeqp, tpbeg, tpend, t1, t2;
	//  lp = latitude of equivalent plane in radians
	//  lat1 = latitude in radians
	//  lat = latitude in degrees

	crad = M_PI/180.0;
	//   crad = degree to radian conversion factor
	//     convert times to radians from noon
	t = (hour-12.0)*M_PI/12.0;
	delt1 = dt*M_PI/12.0;
	//     convert angles to radians
	slope1 = slope*crad;
	azi1 = azi*crad;
	lat1 = lat*crad;
	fjulian = (double)(julian(year, month, day));
	d = crad*23.5*sin((fjulian - 82.0)*0.017214206321);
	//  0.017214206321 is 2 pi / 365
	//  d is solar declination
	lp = asin(sin(slope1)*cos(azi1)*cos(lat1) + cos(slope1)*sin(lat1));
	//  lp is latitude of equivalent plane
	//      td=acos(-tan(lat1)*tan(d))  this formula abandoned 1/8/04
	//      to make the code work for polar conditions
	//  td is half day length, i.e. the time from noon to sunset.  sunrise is at -td
	tanprod = tan(lat1)*tan(d);
	if(tanprod > 1.0) {
		td = M_PI;  // this is the condition for perpetual light
	} else if (tanprod < -1.0) {
		td = 0.0;   // the condition for perpetual night
	} else {
		td = acos(-tanprod);  // the condition where there is a sunrise and set
	}
	//   equivalent longitude offset.  modified on 1/8/04
	//   so that it correctly accounts for shift in longitude if equivalent
	//   plane slope goes over a pole.  achieved using atan2.
	//      ddt=atan(sin(azi1)*sin(slope1)/(cos(slope1)*cos(lat1)
	//     *    -cos(azi1)*sin(slope1)*sin(lat1)))
	ddt = atan2(sin(azi1)*sin(slope1), (cos(slope1)*cos(lat1) - cos(azi1)*sin(slope1)*sin(lat1)));

	//   now similar logic as before needs to be repeated for equivalent plane
	//   but with times reflecting
	tpeqp = tan(lp)*tan(d);
	//  keep track of beginning and end of exposure of equiv plane to sunlight
	if (tpeqp > 1.0) {
		tpbeg = -M_PI;   // perpetual light
		tpend = M_PI;
	} else if (tpeqp < -1.0) {
		tpbeg = 0.0;  // perpetual dark
		tpend = 0.0;
	} else {
		tpbeg = -acos(-tpeqp) - ddt;
		tpend = acos(-tpeqp) - ddt;
	}

	//   start and end times for integration of radiation exposure
	//   need to account for both horizon, slope and time step
	t1 = max(t, max(tpbeg, -td));
	t2 = min(t + delt1, min(td, tpend));
	//      write(6,*)t1,t2
	if (t2 <= t1) {
		hri = 0.0;
	} else {
		hri = (sin(d)*sin(lp)*(t2 - t1) + cos(d)*cos(lp)*(sin(t2 + ddt) - sin(t1 + ddt)))/(cos(slope1)*delt1);
		//   in the above the divide by cos slope normalizes illumination to per unit horizontal area
	}
	//   there is a special case if tpbeg is less than -pi that occurs in polar regions
	//   where a poleward facing slope may be illuminated at night more than the day.
	//   add this in
	if (tpbeg < -M_PI) {
		t1 = max(t, max(-tpbeg + 2*M_PI,-td));
		t2 = min(t + delt1, td);
		if (t2 > t1) {
			hri = hri + (sin(d)*sin(lp)*(t2 - t1) + cos(d)*cos(lp)*(sin(t2 + ddt) - sin(t1 + ddt)))/(cos(slope1)*delt1);
		}
	}
	//  for the purposes of calculating albedo we need a cosine of the illumination angle.  this
	//  does not have slope correction so back this out again.  this is an average over the
	//  time step
	coszen = hri*cos(slope1);
	//      write(6,*)hri,coszen

	return 0;
}


// ***************************** JULIAN () ****************************
//              To convert the real date to julian date
//    The Julian are change to a new version to take the Lean Yean into consideration
//     in the old version, there are 365 days each year.
//      FUNCTION JULIAN(MONTH,DAY)
int julian(const int yy, const int mm, const int dd)
{
	int jday, ileap;
	int mmstrt[12] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};
	jday = mmstrt[mm-1] + dd;
	ileap = yy - (int)(yy/4)*4;
	if (ileap == 0 && mm >= 3)
		jday++;

	return jday;
}

// ****************************** QLIF () *******************************
//     Computes the incoming longwave radiation using satterlund Formula
//      Modified 10/13/94 to account for cloudiness.
//      Emissivity of cloud cover fraction is assumed to be 1.
//
int qlif(double &qliff, const double ta, const double rh, const double tk, const double sbc, const double cf)
{
	double ea, tak, ea1;

	ea = svpw(ta)*rh;
	tak = ta + tk;
	ea1 = 1.08*(1.0-exp(-pow(ea/100.0, (tak)/2016.0) ));
	qliff = (cf + (1.0 - cf)*ea1)*sbc*pow(tak, 4);

	return 0;
}

