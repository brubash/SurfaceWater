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
#include <iomanip>

// DATE & TIME ROUTINES FOR TIDEDA

// Now official Tideda date format is yyyymmdd,
// But it still supports old Tideda yyymmdd format

int td8_mowday(const int ky, const int km, const int kd);
int td8_mowymd(int &ky, int &km, int &kd, const int days);

using namespace std;


int td81micdh(int &idate, int &ihour, const long int jsec)
{
	// convert seconds to date & time (midnight is returned as 240000)

	int y, m, d, ihr, imn, jday;
	long int isc, iday;

	iday = jsec/86400;	// days since 1-jan-1940
	if (jsec > 86400*iday) {
		iday++;
	}
	jday = 14609 + iday;	// days since 1-jan-1900
	td8_mowymd(y, m, d, jday);
	idate = y*10000 + m*100 + d;
	isc = jsec - 86400*(long(jday) - 14610);
	ihr = isc/3600;
	isc = isc  -ihr*3600;
	imn = isc/60;
	isc = isc-imn*60;
	ihour = ihr*10000 + imn*100 + isc;

	return 0;
}

//-------------------------------------------------------------------
int td8micsec(const int jdatem, int &jhour, long int &isec)
{
	// convert date & time to seconds

	// time format is "hhmmss"
	// supported date format is "yyyymmdd", e.g. 20010101

	int idatem, ihour, ihr, imn;
	int d, y, m;
	long int isc;
	double difference;

	idatem = jdatem;

	y = idatem/10000;
	d = idatem - y*10000;
	m = d/100;
	d = d - m*100;

	if (jhour == 0) {
		ihour = 0;
	} else {
		ihour = jhour;
	}
	ihr = ihour/10000;
	isc = ihour - ihr*10000;
	imn = isc/100;
	isc = isc - imn*100;
	if (ihr == 24) {
		d += 1;
		ihr = 0;
	}


	struct tm rTime, fileTime;
	time_t refTime, thisTime;
	rTime.tm_year  = 40;	// from 1900
	rTime.tm_mon   = 0;	// 0 - 11
	rTime.tm_mday  = 1;	// 1 - 31
	rTime.tm_hour  = 0;	// 0 - 23
	rTime.tm_min   = 0;	// 0 - 59
	rTime.tm_sec   = 0;	// 0 - 60
	rTime.tm_isdst = 0;

	fileTime.tm_year  = y-1900;
	fileTime.tm_mon   = m-1;
	fileTime.tm_mday  = d;
	fileTime.tm_hour  = ihr;
	fileTime.tm_min   = imn;
	fileTime.tm_sec   = isc;
	fileTime.tm_isdst = 0;
	refTime  = mktime(&rTime);
	thisTime = mktime(&fileTime);

	difference = long(difftime(thisTime, refTime));
	isec = long(difference);

	return 0;
}

	//-----------------------------------------------------------------
	// The guts of the Tideda time algorithm:
	//-----------------------------------------------------------------
int td8_mowday(const int ky, const int km, const int kd)
{
	//	to calculate the number of days since 1st jan 1900
	//	from the year month and day

	int iya, ic;
	int y;
	int m;

	y = ky;
	m = km;

	// adjust months & years to give months in the range 1 to 12
	while (m > 12) {
		if (m > 12) {
			m -= 12;
			y++;
		}
	}

	m -= 3;
	if (m < 0) {
		y--;
		m += 12;
	}
	ic = y/100;
	iya = y - 100*ic;
	return (146097*ic)/4 + (1461*iya)/4 + (153*m + 2)/5 + kd - 693901;
}

	//-----------------------------------------------------------------
int td8_mowymd(int &ky, int &km, int &kd, const int days)
{
	//	returns the year,month,day given the number of days since 1900
	int j, y, m, d;
	j = 4*(693901 + days) - 1;
	y = j/146097;
	j = j - 146097*y;
	d = j/4;
	m = 4*d + 3;
	j = m/1461;
	d = m - 1461*j;
	d = (d + 4)/4;
	y = 100*y + j;
	j = 5*d - 3;
	m = j/153;
	d = j - 153*m;
	d = (d + 5)/5;
	m = m + 3;
	if(m <= 12)
		goto L1002;
	m = m - 12;
	y = y + 1;
	//         calc algorithm stops
L1002:	ky = y;
	km = m;
	kd = d;

	return 0;
}

