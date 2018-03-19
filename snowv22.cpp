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

using namespace std;

//  File  snowx.f
//  Subroutines and function subprograms for the Utah Energy Balance
//  Snow Accumulation and Melt Model.
//  David G. Tarboton, Utah Water Research Laboratory, Utah State University
//  Version 2.1
//  Last Change 5/4/04 to accommodate glacier surface runoff.

//  This program was written and prepared under agreement with and funding
//  by the U.S. Government , and therefore is in the public domain and not
//  subject to copyright.

// **************READ Bristow Campbell Parameter values **************

int bcparm(double dtBar[], double &bca, double &bcc, const string bcFile)
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

// **************READ PARAMETER or state variable VALUES**************
/*
      SUBROUTINE readvals(vals,vfile,n)
      CHARACTER*200 vfile
      REAL vals(1)
      OPEN(8,FILE=vfile,STATUS='OLD')
      READ(8,*)(vals(i),i=1,n)
      CLOSE(8)
      RETURN
      END

C********************* READ THE INITIAL CONDITIONS *********************

C      SUBROUTINE FINITIAL(statev,icfile,nxv,YEAR,MONTH,DAY,HOUR,DT)
C      CHARACTER*20 icfile
C      INTEGER YEAR,DAY
C      REAL statev(1)
C      OPEN(8,FILE=icfile,STATUS='OLD')
C      READ(8,*)(statev(i),i=1,nxv)
C      READ(8,*)YEAR,MONTH,DAY,HOUR,DT
C      CLOSE(8)
C      RETURN
C      END

C********************* READ THE site variables  *********************

      SUBROUTINE readsv(sitev,svfile,nsv,slope,azi,lat)
      CHARACTER*200 svfile
      REAL sitev(1),lat
      OPEN(8,FILE=svfile,STATUS='OLD')
      READ(8,*)(sitev(i),i=1,nsv)
      READ(8,*)slope,azi,lat
      CLOSE(8)
      RETURN
      END
*/
// ************************** updatetime () ***************************
//                 update time for each time step
int lyear(const int year);

int updatetime(int &year, int &month, int &day, double &hour, const double dt)
{
	int dm, i;

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
	int result;
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
/*
      Subroutine atf(atff,trange,month,dtbar,a,c)
      DIMENSION dtbar(12)
      b=0.036*exp(-0.154*dtbar(month))
      atff=a*(1-exp(-b*trange**c))
C      write(6,*)trange,month,a,c,dtbar(month),atf
      RETURN
      END


C************************** hourlyRI () **********************
C                To get hourly radiation index

      SUBROUTINE hyri(YEAR,MONTH,DAY,HOUR,DT,SLOPE,AZI,LAT,
     *                    HRI,COSZEN)
      INTEGER YEAR,DAY
      REAL LP,LAT1,LAT
c  lp= latitude of equivalent plane in radians
c  lat1 = latitude in radians
c  lat = latitude in degrees
c  a number that speaks for itself - every kissable digit
      PI=3.141592653589793238462643383279502884197169399375105820974944592308
      CRAD=PI/180.0
c   crad = degree to radian conversion factor
C     CONVERT TIMES TO RADIANS FROM NOON
      T=(HOUR-12.0)*PI/12.0
      DELT1=DT*PI/12.0
C     CONVERT angles TO RADIANS
      SLOPE1=SLOPE*CRAD
      AZI1=AZI*CRAD
      LAT1=LAT*CRAD
	FJULIAN=FLOAT(JULIAN(year,MONTH,DAY))
      D=CRAD*23.5*SIN((FJULIAN-82.0)*0.017214206321)
c  0.017214206321 is 2 pi / 365
c  D is solar declination
      LP=ASIN(SIN(SLOPE1)*COS(AZI1)*COS(LAT1)
     *   +COS(SLOPE1)*SIN(LAT1))
c  LP is latitude of equivalent plane
c      TD=ACOS(-TAN(LAT1)*TAN(D))  This formula abandoned 1/8/04
c      to make the code work for polar conditions
c  TD is half day length, i.e. the time from noon to sunset.  Sunrise is at -TD
      tanprod=TAN(LAT1)*TAN(D)
	if(tanprod .gt. 1.){
	  td=pi  ! This is the condition for perpetual light
	else if(tanprod .lt. -1.){
	  td=0   ! The condition for perpetual night
	else
	  td=acos(-tanprod)  ! The condition where there is a sunrise and set
	endif
c   Equivalent longitude offset.  Modified on 1/8/04
c   so that it correctly accounts for shift in longitude if equivalent
c   plane slope goes over a pole.  Achieved using atan2.
c      DDT=ATAN(SIN(AZI1)*SIN(SLOPE1)/(COS(SLOPE1)*COS(LAT1)
c     *    -COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))
      ddt=atan2(SIN(AZI1)*SIN(SLOPE1),
     *(COS(SLOPE1)*COS(LAT1)-COS(AZI1)*SIN(SLOPE1)*SIN(LAT1)))

c   Now similar logic as before needs to be repeated for equivalent plane
c   but with times reflecting
      TPeqp=TAN(LP)*TAN(D)
C  Keep track of beginning and end of exposure of equiv plane to sunlight
      IF(tpeqp .gt. 1.0) {
          TPbeg=-pi   ! perpetual light
	    tpend=pi
      ELSEif(tpeqp .lt. -1.){
          TPbeg=0.0  ! perpetual dark
          tpend=0.0
	else
	    tpbeg = -acos(-tpeqp)-ddt
	    tpend = acos(-tpeqp)-ddt
      ENDIF

c   Start and end times for integration of radiation exposure
c   need to account for both horizon, slope and time step
      T1=AMAX1(T,tpbeg,-TD)
      T2=AMIN1(T+DELT1,TD,tpend)
C      write(6,*)t1,t2
      IF(T2.LE.T1) {
        HRI=0.0
      ELSE
        HRI=(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT)
     *       -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
c   In the above the divide by cos slope normalizes illumination to per unit horizontal area
      ENDIF
c   There is a special case if tpbeg is less than -pi that occurs in polar regions
c   where a poleward facing slope may be illuminated at night more than the day.
C   Add this in
      if(tpbeg .lt. -pi){
	  T1=AMAX1(T,-tpbeg+2*pi,-TD)
	  T2=AMIN1(T+DELT1,TD)
	  if(t2 .gt. t1){
	    hri=hri+(SIN(D)*SIN(LP)*(T2-T1)+COS(D)*COS(LP)*(SIN(T2+DDT)
     *       -SIN(T1+DDT)))/(COS(SLOPE1)*DELT1)
	  endif
	endif
c  for the purposes of calculating albedo we need a cosine of the illumination angle.  This
c  does not have slope correction so back this out again.  This is an average over the
c  time step
      COSZEN = HRI*COS(SLOPE1)
C      write(6,*)hri,coszen

      RETURN
      END


C***************************** JULIAN () ****************************
C              To convert the real date to julian date
cYJS  The Julian are change to a new version to take the Lean Yean into consideration
c     in the old version, there are 365 days each year.
c      FUNCTION JULIAN(MONTH,DAY)
      function julian(yy,mm,dd)
      integer yy,dd
      dimension mmstrt(1:12)
      data (mmstrt(i),i=1,12)/0,31,59,90,120,151,181,212,243,273,
     *                        304,334/
      jday = mmstrt(mm) + dd
      ileap = yy - int(yy/4)*4
      if(ileap.eq.0.and.mm.ge.3) jday = jday + 1
      julian = jday
      return
      end

C ****************************** QLIF () *******************************
C     Computes the incoming longwave radiation using satterlund Formula
C      Modified 10/13/94 to account for cloudiness.
C      Emissivity of cloud cover fraction is assumed to be 1.
C
      subroutine qlif(qliff,TA,RH,TK,SBC,cf)
      EA = SVPW(TA)*RH
      TAK = TA + TK
      EA1 = 1.08*(1.0-EXP(-(EA/100.0)**((TAK)/2016.0)))
      QLIFf =(cf+(1.-cf)*EA1)*SBC*TAK**4
      RETURN
      END
*/
