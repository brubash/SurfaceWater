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
#include <iomanip>
#include <sstream>

using namespace std;
using namespace Eigen;

ofstream topErrorFile;

// ***********************************************************************

// This version, V2, has enhanced checking for missing data in the RAMS input and rainfilling
// V3 also has checking/matching of the flow sites with reaches.
int hyData(int &sDate, int &sHour, long &interval, int &m, int &mi, int &mps, int &mpe, int &Ngauge, int &Neq,
	ArrayXXd &bRain, double **flow, int &iret, ArrayXd &temper, double *dewp, double *trange, double **dtBar,
	const int Ns, ArrayXXd &wrg, ArrayXXi &lrg, const double *elevtg, double **bTmax, double **bTmin, double
	**bTdew, double **bdtBar, const ArrayXXd &Sp, const int maxGauge, const int maxInt, const int maxSites,
	const int maxResponse, const int maxTGauge, double *wind2m, ArrayXXd &wrg1,
	int &idebugoutput, int &idebugbasin, int &idebugcase)
{
	int date, hour;
	int  ntri, jj, ij;
	int  js, kg;

	long itemp1, itemp2, itemp3;
	int kk;
	ArrayXd tempr(maxSites), tempr_last(maxSites);
	ArrayXXd tempt(maxSites,3);
	double *tempf, *tempf_last;

	int i, j, nwind, wsite[1];
	//	integer*4 dsite(maxsites) dsite holds the flow site no's
	int nfill[maxGauge];
	int rsite[maxGauge], fillflag=0;
	int tsite[maxGauge];
	string rainfill_file, lineString;
	istringstream line;
	string verno;
	int it_save_old;
	bool i_reset_flag=false, exist;
	struct dirent *dirp;
	DIR *dp;
	string testStr, inLine;

	ArrayXXd rcoeff(maxGauge,maxGauge);
	ArrayXXi fsite(maxGauge,maxGauge);

	//	REAL*4 RAIN_FACTOR           defined in globaly in namespace rain
	//	COMMON /RAIN1/ RAIN_FACTOR

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> hyData(" << ncalls << ")" << std::endl;
    }
	caller = "hyData";
#endif

    topErrorFile.open("results/toperror.txt");

	//------------ topinp.dat -------------------
	std::istringstream iss;
	ifstream topinpFile;
	string topinpFileName;
#ifdef WRIA1
	topinpFileName = "topinpWRIA1.dat";
#else
    topinpFileName = "topinp.dat";
#endif
    topinpFile.open(topinpFileName);
	if (topinpFile.good()) {
        cout << topinpFileName << " opened for reading in hyData()\n";
    } else {
        cout << topinpFileName << " not found\n";
        exit(1);
    }
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> sDate;
    iss.clear();
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> sHour;
    iss.clear();
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> interval;
    iss.clear();
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> m;
    if (m > maxInt) {
        cerr << " The simulation length of " << dec << setw(5) << m;
        cerr << " timesteps exceeds the program limit of " << dec << setw(7) << maxInt;
        cerr << " The program limit will be used.\n";
    }
    if (m <= 0) {
        cerr << " The simulation length < 1, see topinp.dat\n";
        exit(0);
	}
	iss.clear();
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> mps;
    if (mps <= 0) {
        mps = 1;
    }
    if (mps > m) {
        mps = m;
    }

    iss.clear();
    getline(topinpFile, inLine, '\n');
    iss.str(inLine);
    iss >> mpe;
    idebugbasin = 0;
	idebugcase  = 0;
    getline(topinpFile, inLine, '\n');
    iss.clear();
    iss.str(inLine);
	iss >> idebugoutput;
	if(idebugoutput != 0) {
        iss >> idebugbasin;
        iss >> idebugcase;
    }
    topinpFile.close();
	//------------ topinp.dat -------------------

	td8micsec(sDate, sHour, itemp2);

	mi = 0;
	for (i = 1; i <= maxGauge; i++) {
		rsite[i-1] = i;
	}
	string rainFileName;
#ifdef WRIA1
    rainFileName = "rain_allWRIA1.dat";
#else
    rainFileName = "rain.dat";
#endif
	rainFile.open(rainFileName);		// fortran unit lung
	if (!rainFile.is_open()) {
		cerr << "Failed to open " << rainFileName << '\n';
		exit(EXIT_FAILURE);
	} else {
        cerr << rainFileName << " opened for reading.\n";
    }
	getline(rainFile, inLine, '\n');
	getline(rainFile, inLine, '\n');
	rainFile >> verno;

	if (verno == "Ver1" || verno == "ver1" || verno == "Ver2" || verno == "ver2") {
		rainFile >> ntri;
		for (i = 1; i <= ntri; i++) {
			rainFile >> rsite[i-1];
		}
	} else {
		for (i = 1; i <= ntri; i++) {
			rainFile >> rsite[i-1];
		}
		rainFile >> Ngauge;
		ntri = Ngauge;	// need to define ntri for changes below
	}
	rainFile.close();

	for (i = 1; i <= maxGauge; i++) {
		for (j = 1; j <= maxGauge; j++) {
			fsite(i-1,j-1) = 0;
		}
		nfill[i-1] = 0;
	}

	rainfill_file = "rainfill.txt";
	exist = false;
	if ((dp = opendir(".")) == NULL) {
		cerr << "Can't open directory\n";
		exit(EXIT_FAILURE);
	}
	while ((dirp = readdir(dp)) != NULL) {
		testStr = dirp->d_name;
		if (testStr == rainfill_file) {
			exist = true;
			closedir(dp);
			break;
		}
	}
	if (exist == false) {
		cerr << " No " << rainfill_file << " file found, continuing without it.\n";
	} else {
		rainfill_read(rainfill_file, ntri, rsite, rcoeff, fsite, nfill, fillflag, maxGauge);
		if (fillflag < 0) {
            cerr << " Empty " << rainfill_file << " file found, continuing without it.\n";
		}
	}
	date = 0;
	hour = 0;
	rainFile.open(rainFileName);		// fortran unit lung
	if (!rainFile.is_open()) {
		cerr << "Failed to open " << rainFileName << '\n';
		exit(EXIT_FAILURE);
	} else {
        cerr << rainFileName << " opened for reading\n";
	}
	getline(rainFile, inLine, '\n');
	getline(rainFile, inLine, '\n');
	getline(rainFile, inLine, '\n');
	i = 0;
L203: for (jj = 0; jj < ntri; jj++) {
		rainFile >> tempr(jj);
	}
	rainFile >> date >> hour;
    td8micsec(date, hour, itemp1);
	if (i == 0) {
		// data starts after given start - flag as an error
		if ( itemp2 < itemp1-interval) {
			cerr << " Start time of the rainfall, or runoff data " << dec << setw(9) << date;
			cerr << dec << setw(7) << hour;
			cerr << " before start of filed data at " << dec << setw(9) << sDate;
			cerr << dec << setw(7) << sHour << '\n';
			exit(EXIT_FAILURE);	// iret != 0 in the fortran code stops the program
		}
		// present data is before given start - go and read the next value
		if (itemp2 > itemp1)
			goto L203;
	}
	i++;
	// added some code here to allow for breaks in the rain.dat file
	if (i == 1) {  // First time through so need to set up values for itemp3 & TEMPR_LAST
		itemp3 = itemp1 - interval;
	}

	if (itemp1-itemp3 != interval) {
		if (itemp1-itemp3 < interval) {
			cerr << " Check data in " << rainFileName << " near " << date << " " << hour << '\n';
			exit(EXIT_FAILURE);
		}
		cerr << " Data missing in " << rainFileName << " near " << date << " " << hour << '\n';
		cerr << " Missing data replaced with last given value\n";
		itemp3 += interval;
		for (jj = 0; jj < ntri; jj++) {
			tempr(jj) = tempr_last(jj);
		}
	}

		// Fill missing rainfall for this timestep: TEMPr overwritten
		// don't try to fill missing data unless we found a rainfill.txt
		// need to come here even if rainfall ends prematurely to ensure filling.
	// ******************** ignore filling gaps initially RPI 17/3/2004
//L220:
    if (fillflag > 0) {
		rainfill_doit(tempr, i, ntri, rcoeff, fsite, nfill, i_reset_flag, it_save_old, maxGauge);
	}

	// Calculate basin rainfalls - code transferred from CALCTS
	for (js = 0; js < Ns; js++) {				// for each subbasin
		bRain(js,i-1) = 0.0;		 			// "brain" has basin rain time series
		for (kg = 1; kg <= maxGauge; kg++) {	// get brain by adding weighted "rain"
			if (lrg(js,kg-1) > 0) {
				bRain(js,i-1) += tempr(lrg(js,kg-1)-1)*wrg(js,kg-1); // mm to um
			}
		}
	}

	itemp3 = itemp1;	// Finished filling, so update itemp3 AND TEMPR_LAST in case neede at next time step
	for (jj = 0; jj < ntri; jj++) {
		tempr_last(jj) = tempr(jj);
	}
	// reset number of values if missing data at end of record
	if ( i_reset_flag ) {
		m = it_save_old;
		cerr << " +++++++++++ WARNING WARNING ++++++++++++\n";
		cerr << "***** Owing to missing rainfall at all stations to the end of the record\n";
		cerr << " it is not possible to reliably estimate flows beyond interval " << dec << setw(5) << it_save_old << ",\n";
		cerr << " therefore M reset accordingly\n";
	}
	if ( i < m  )
		goto L203;
	rainFile.close();
	// RAIN_FACTOR read in MDDATA
	// from MODELSPC.DAT and the multiplication is done here to because we have
	// to wait until Ross has done any filling required.
	for (i = 0; i < m; i++) {
		for (jj = 0; jj < Ns; jj++) {
			bRain(jj,i) *= rain::rain_factor;
		}
	}

	date = 0;
	hour = 0;

	ifstream windFile("wind.dat");		// fortran unit lung
	if (!windFile.is_open()) {
		cerr << "Failed to open wind.dat\n";
		exit(EXIT_FAILURE);
	}

	getline(windFile, inLine, '\n');
	getline(windFile, inLine, '\n');
	windFile >> verno >> nwind;
	for (i = 1; i <= nwind; i++) {
		windFile >> wsite[i-1];
	}
	getline(windFile, inLine, '\n');	// read the rest of the line
	i = 0;
L3103: if (!(windFile >> tempr(0) >> date >> hour))
		goto L3101;

	td8micsec(date, hour, itemp1);
	if ( i == 0 ) {
		// data starts after given start - flag as an error
		if ( itemp2 < itemp1-interval )
			goto L3100;
		// present data is before given start - go and read the next value
		if ( itemp2 > itemp1 )
			goto L3103;
	}
	i++;
	wind2m[i-1] = tempr(0);
	if (i == 1) {  // First time through so need to set up values for itemp3 & TEMPR_LAST
		itemp3 = itemp1 - interval;
	}
	if (itemp1-itemp3 != interval) {
L1682:	if (itemp1-itemp3 < interval) {
				cerr << " Check data in wind.dat near " << date << " " << hour << '\n';
				exit(EXIT_FAILURE);
		}
		cerr << " Data missing in wind.dat near " << date << " " << hour << '\n';
		cerr << " Missing data replaced with last given value\n";
		itemp3 = itemp3 + interval;
		i++;
		if (itemp1-itemp3 != interval)
			goto L1682;
	}
	itemp3 = itemp1; // Finished filling, so update itemp3 AND TEMPR_LAST in case neede at next time step

	if (i < m)
		goto L3103;
	windFile.close();
	goto L3104;
L3100: cerr << " ***** Warning: wind.dat does not exist or contains an error at ";
	cerr << dec << setw(9) << date << setw(7) << hour << '\n';
	cerr << " Proceeding with constant temperature\n";
L3101: if (m != i) {
		m = i;
        cerr << "Warning: Temperature file Ends prematurely. Number of time steps reduced to " << m << '\n';
	}
	if ( m <= 0 ) {
		cerr << " **** Change to length of temperature data has caused";
		cerr << " a value of data length < 1 - check simulation dates in";
		cerr << " topinp.dat and temper.dat\n";
		exit(EXIT_FAILURE);
	}
L3104: ;
	// put in default temperatures if there's no temperature data
	for (i = 1; i <= m; i++) {
		temper(i-1) = 10.0;
		dewp[i-1]   = 7.0;
		trange[i-1] = 10.0;
	}
	i = m;
	date = 0;
	hour = 0;
	string tmaxtmintdewFileName;
#ifdef WRIA1
    tmaxtmintdewFileName = "tmaxtmintdew_allWRIA1.dat";
#else
    tmaxtmintdewFileName = "tmaxtmintdew.dat";
#endif
	ifstream tmaxtmintdewFile(tmaxtmintdewFileName);		// fortran unit lung
	if (!tmaxtmintdewFile.is_open()) {
		cerr << "Failed to open " << tmaxtmintdewFileName << '\n';
		exit(EXIT_FAILURE);
	}

	getline(tmaxtmintdewFile, inLine, '\n');
	getline(tmaxtmintdewFile, inLine, '\n');
	tmaxtmintdewFile >> verno >> ntri;
	for (i = 1; i <= ntri; i++) {
		tmaxtmintdewFile >> tsite[i-1];
	}
	getline(tmaxtmintdewFile, inLine, '\n');	// read the rest of the line
	i = 0;
L303: for (jj = 1; jj <= ntri; jj++) {
		for (js = 1; js <= 3; js++) {
			if (!(tmaxtmintdewFile >> tempt(jj-1,js-1)))
				goto L301;
		}
	}
	tmaxtmintdewFile >> date >> hour;
	td8micsec(date, hour, itemp1);
	if ( i == 0 ) {
		// data starts after given start - flag as an error
		if ( itemp2 < itemp1-interval )
			goto L300;
		// present data is before given start - go and read the next value
		if ( itemp2 > itemp1 )
			goto L303;
	}
	i++;
	// added some code here to allow for breaks in the input file
	if (i == 1) {  // First time through so need to set up values for itemp3 & TEMPR_LAST
		itemp3 = itemp1 - interval;
	}
	if (itemp1-itemp3 != interval) {
L682: 	if (itemp1-itemp3 < interval) {
			cerr << " Check data in tmaxtmintdew.dat near " << date << hour << '\n';
			exit(EXIT_FAILURE);
		}
		cerr << " Data missing in tmaxtmintdew.dat near " << date << hour << '\n';
		cerr << " Missing data replaced with last given value\n";
		itemp3 = itemp3 + interval;
		//  4/23/07  DGT added code below to implement missing within IF
		for (js = 1; js <= Ns; js++) {			 // for each subbasin
			bTmax[js-1][i-1] = 0.0;		 // "bTmax" has basin tmax time series
			bTmin[js-1][i-1] = 0.0;		 // "bTmin" has basin tmax time series
			bTdew[js-1][i-1] = 0.0;		 // "bTdew" has basin tmax time series
			for (kg = 1; kg <= maxTGauge; kg++) { // get bTmax by adding using lapse to get from gauge_elev to basin_elev, plus spatial weighting on tmax
				if (lrg(js-1,kg-1) > 0) {  // Interpolation weights
            	    bTmax[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1),0) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
            	    bTmin[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1),1) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
            	    bTdew[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1),2) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
				}
			}
		}
		i++;
		if (itemp1-itemp3 != interval)
			goto L682;
	}
	// Calculate basin temps - code adapted from basin rain raw 9-dec-2004
 	for (js = 1; js <= Ns; js++) {			 // for each subbasin
		bTmax[js-1][i-1] = 0.0;		 // "bTmax" has basin tmax time series
		bTmin[js-1][i-1] = 0.0;		 // "bTmin" has basin tmax time series
		bTdew[js-1][i-1] = 0.0;		 // "bTdew" has basin tmax time series
		for (kg = 1; kg <= maxTGauge; kg++) { // get bTmax by adding using lapse to get from gauge_elev to basin_elev, plus spatial weighting on tmax
			if (lrg(js-1,kg-1) > 0) {  //  Interpolation weights
				bTmax[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1)-1,0) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
				bTmin[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1)-1,1) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
				bTdew[js-1][i-1] += wrg1(js-1,kg-1)*(tempt(lrg(js-1,kg-1)-1,2) + Sp(12,js-1)*(elevtg[lrg(js-1,kg-1)-1] - Sp(13,js-1)));
			}
		}
 	}
	itemp3 = itemp1; // Finished filling, so update itemp3 and tempr_last in case neede at next time step
	if ( i < m )
		goto L303;

	for (js = 1; js <= Ns; js++) {
		for (i = 0; i < 12; i++) {
			bdtBar[i][js-1] = 0;
		}
		for (kg = 1; kg <= maxTGauge; kg++) { // get basin average dtbar by averaging dtbar values for gauges
			if (lrg(js-1,kg-1) > 0) {
				for (kk = 1; kk <= 12; kk++) {
					bdtBar[kk-1][js-1] += wrg1(js-1,kg-1)*dtBar[kk-1][lrg(js-1,kg-1)-1];  // WRG1 Interpolation weights
				}
			}
		}
	}

	tmaxtmintdewFile.close();
	goto L304;
L300: cerr << " ***** Warning: temper.dat does not exist or contains an error at ";
	cerr << dec << setw(9) << date << setw(7) << hour << '\n';
	cerr << " Proceeding with constant temperature\n";
L301: if (m != i) {
		m = i;
        cerr << "Warning: Temperature file Ends prematurely\n";
		cerr << "number of time steps set to " << m << '\n';
	}
	if ( m <= 0 ) {
		cerr << " **** Change to length of temperature data has caused";
		cerr << " a value of data length < 1 - check simulation dates in";
		cerr << " topinp.dat and temper.dat\n";
		exit(EXIT_FAILURE);
	}
L304: ;

	// put in zero runoff in case there's no runoff data (e.g. flood forecasting)
	for (jj = 1; jj <= maxResponse; jj++) {
		for (i = 1; i <= m; i++) {
	    	flow[jj-1][i-1] = 0.0;
		}
	}
	i = m;
	date = 0;
	hour = 0;

	ifstream streamflow_calibrationFile("streamflow_calibration.dat");		// fortran unit lung
	if (!streamflow_calibrationFile.is_open()) {
		cerr << "Failed to open streamflow_calibration.dat\n";
		goto L2031;
	}

	getline(streamflow_calibrationFile, inLine, '\n');
	getline(streamflow_calibrationFile, inLine, '\n');
	getline(streamflow_calibrationFile, inLine, '\n');
	streamflow_calibrationFile >> verno;
	// rewind
	streamflow_calibrationFile.clear();
	streamflow_calibrationFile.seekg(0);

	if (verno == "Ver1" || verno == "ver1" || verno == "Ver2" || verno == "ver2") {
		getline(streamflow_calibrationFile, inLine, '\n');
		getline(streamflow_calibrationFile, inLine, '\n');
		getline(streamflow_calibrationFile, inLine, '\n');
		streamflow_calibrationFile >> verno >> ntri;
		getline(streamflow_calibrationFile, inLine, '\n');	// read the rest of the line
		if (ntri != Neq) {
			cerr << " *****ERROR - Number of sites in streamflow_calibration.dat " << ntri;
			cerr << " does not match what modelspc.dat is expecting, " << Neq << "\n";
			cerr << " Setting Neq = ntri\n";
			Neq = ntri;
		}
		tempf = new double[ntri];
		tempf_last = new double[ntri];
		for (jj = 0; jj < ntri; jj++) {
			tempf[jj] = 0.0;
			tempf_last[jj] = 0.0;
		}
	} else {
		getline(streamflow_calibrationFile, inLine, '\n');
	}
	i = 0;
L207: for (jj = 1; jj <= Neq; jj++) {
		if (!(streamflow_calibrationFile >> tempf[jj-1]))
			goto L206;
	}
	if (!(streamflow_calibrationFile >> date >> hour))
		goto L206;

	td8micsec(date, hour, itemp1);
	topErrorFile << dec << setw(4) << i;
	topErrorFile << dec << setw(4) << Neq;
	topErrorFile << dec << setw(6) << m;
	topErrorFile << fixed << setw(9) << setprecision(3) << tempf[0];
	topErrorFile << fixed << setw(9) << setprecision(3) << tempf[1];
	topErrorFile << dec << setw(10) << date;
	topErrorFile << dec << setw(8) << hour;
	topErrorFile << dec << setw(12) << itemp1;
	topErrorFile << dec << setw(12) << itemp2 << endl;

	if ( i == 0 ) {
		// data starts after given start - flag as an error
		if ( itemp2 < itemp1-interval )
			goto L230;
		// present data is before given start - go and read the next value
		if ( itemp2 > itemp1 )
			goto L207;
	}
	i++;
	// added some code here to allow for breaks in the RUNOFF.DAT file
	if (i == 1) {  // First time through so need to set up values for itemp3 & tempf_LAST
			itemp3 = itemp1 - interval;
	}
	if (itemp1-itemp3 != interval) {
L782:	if (itemp1-itemp3 < interval) {
			cerr << " Check data in streamflow_calibration.dat near ";
			cerr << dec << setw(9) << date << setw(7) << hour << '\n';
			exit(EXIT_FAILURE);
		}
		cerr << " Data missing in streamflow_calibration.dat near ";
		cerr << dec << setw(9) << date << setw(7) << hour << '\n';
		cerr << " Missing data replaced with last given value\n";
		itemp3 = itemp3 + interval;
		for (jj = 1; jj <= Neq; jj++) {
			flow[jj-1][i-1] = tempf_last[jj-1]/1000.0;
		}
		i++;
		if (itemp1-itemp3 != interval)
			goto L782;
	}
	itemp3 = itemp1; // Finished filling, so update itemp3 AND tempf_LAST in case neede at next time step
	for (jj = 1; jj <= Neq; jj++) {
		tempf_last[jj-1] = tempf[jj-1];
	}

	for  (jj = 1; jj <= Neq; jj++) {
		flow[jj-1][i-1] = tempf[jj-1]/1000.0;
	}
	if ( i < m  ) {
		goto L207;
	}
	streamflow_calibrationFile.close();
	//	close(20) not closed so that top1 can read from this
L2031: if (m != i) {
		cerr << "***** Warning: streamflow_calibration.dat does not exist or contains an error at ";
		cerr << dec << setw(9) << date << setw(7) << hour << '\n';
		cerr << " Proceeding with measured runoff = 0\n";
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving hyData @ 1(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif
	return 0;
//L250:  cerr << "Error or End of file reading topinp.dat\n";
//	iret = 4;
//	return 0;
//L200:  cerr << "***** Error in rain.dat or streamflow_calibration.dat at ";
//	cerr << dec << setw(9) << date << setw(7) << hour << '\n';
//	iret = 4;
//	return 0;
//L201: cerr << " ***** Rainfall file is empty\n";
//	iret = 4;
//	return 0;
//L202: m = i;
//	if ( m <= 0 ) {
//		cerr << " **** Premature end to data has caused a value of data";
//		cerr << " length < 1 - check simulation dates in topinp.dat and/or";
//		cerr << " the dates in rain.dat and streamflow_calibration.dat";
//		exit(EXIT_FAILURE);
//	}
//	cerr << " ***** Rainfall ends prematurely - number of values reduced to " << dec << setw(6) << m << '\n';
//	goto L220;
//L204: cerr << " ***** Error in file streamflow_calibration.dat at ";
//	cerr << dec << setw(9) << date << setw(7) << hour << '\n';
//	iret = 4;
//	for (i = 0; i < maxGauge; i++) {
//		delete [] rcoeff[i];
//		delete [] fsite[i];
//	}
//	delete rcoeff;
//	delete fsite;
//	return 0;
L206: for (ij = i+1; ij <= m; ij++) {
 		for (jj = 1; jj <= Neq; jj++) {
 			flow[jj-1][ij-1] = -1.0;
 		}
	}
 	cerr << " ***** RUNOFF DATA ENDS PREMATURELY - Flows set to -1 - Ignore this message if forecasting\n";
	if ( m <= 0 ) {
		cerr << " **** Premature end to data has caused a value of data";
		cerr << " length < 1 - check simulation dates in topinp.dat and/or";
		cerr << " the dates in rain.dat and streamflow_calibration.dat";
		exit(EXIT_FAILURE);
	}
	// need to close this because we will be opening it again later
	streamflow_calibrationFile.close();
L230:
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving hyData @ 2(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}
