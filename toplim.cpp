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

//  This version, V4, automatically estimates zbar0 from initial flows

#include "topnet.hh"
#include <dirent.h>
#include <iomanip>

// Subroutines to link topmodel to NLFIT
//  The purpose of this routine is to read input and parameter data needed to
//  run TOPMODEL within nlfit, as well as outputs to be used in calibration.
//  The first 7 variables are required outputs returned to NLFIT for calibration.
//  nrx and iex are array dimensions. Common is used to pass input and parameter
//  data to the subroutine MODEL.

//  The main output is qact, a two dimensional array of data for use in calibration.
//  this has dimension (nrx,iex) with each column, j, representing a response vector
//  with values from init(j) to ined(j).  Times of rows in qact are given by actime.
//  Npar gives the number  of model parameters that NLFIT may adjust in calibration
//  modelid gives a model name.

using namespace std;
using namespace Eigen;

// Fortran common (to inputT() and Model() )
ofstream lunmodFile;
ofstream luntopFile;
ofstream lundatFile;
ofstream lunpFile;

int maxInt, maxGauge, maxSlp, maxTGauge;
int maxResponse, maxA, maxC, maxChn;
int maxSites, maxRchAreas;
int max_lakes, max_lheads;
int max_of_shifts;

int nRchSav;

// The above variables are declared extern to ensure these are saved outside of file scope.

int inputT(int initT[], int iend[], int &Neq, ArrayXXd &Qact, double acTime[],
	string &modelId, int &Npar, const int Nrx, const int iex)
{
	int n, topinp;
	int chn;
	int  maxA2, maxA1, maxC2, maxC1, Neq0, Nout0, lineCount, Ntri;
	int  Rpairs1, Rpairs2=0;

	string text;
	string verno;

	bool exist;
	int *ll, *Ntr, *llOut, iret;
	int *Nts;
	double **Si, **Rp, **atb, **pka, *tl, **pd, **cl;
	double units;
	double **Sid, **Rpd;	// The arrays Sid, Spd and Rpd record the original values.
	int *Nka, *Nd, *iBout;
	int maxFgauge;
	int relFlag, neq_temp;
	int *Qmap;
	double *rel;	// for forecasting e
	long int stim8;

	int *ishift0;
	double **bTmax;
	double **bTmin;
	double  *wind2m;
	double **bTdew;
	double  *bXlat;
	double  *bXlon;
	double **flow;
	double dummyq;

	//   ET variables  DGT
	double  *Xlat, *Xlong;
	double stdlon;
	double  *elevtg, **dtBar, **bdtBar;
	int *temper_id;
	double *dewp, *trange;
	int ns_temper, idebugoutput, idebugbasin, idebugcase;
	int i, j, k;
	struct dirent *dirp;
	DIR *dp;
	string testStr, inLine;

	// ***********************************************************************
#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> inputT(" << ncalls << ")" << std::endl;
    }
    caller = "inputT";
#endif
	// Read values from data files
	// iDoIt tells model what operations it is to do only once
	model_master::iDoIt = 0;

	//------------ topinp.dat -------------------
	ifstream topinpFile;
	topinpFile.open("topinp.dat");
	if (topinpFile.good()) {
        cout << "topinp.dat opened for reading in inputT()\n";
    } else {
        cout << "topinp.dat not found\n";
        exit(1);
    }
    getline(topinpFile, inLine, '\n');
    getline(topinpFile, inLine, '\n');
    getline(topinpFile, inLine, '\n');
    topinpFile >> topinp;
    topinpFile.close();
    maxInt = topinp + 200;

	//------------ modelspc.dat -------------------
	modelspcFile.open("modelspc.dat");  // Fortran unit 2
	if (modelspcFile.good()) {
        cout << "modelspc.dat opened for reading in inputT()\n";
    } else {
        cout << "modelspc.dat not found\n";
        exit(1);
    }
    getline(modelspcFile, inLine, '\n');
    getline(modelspcFile, inLine, '\n');
    modelspcFile >> maxGauge >> maxSlp >> chn;
    maxChn   = maxSlp + chn;

    maxA2 = 0;
	maxC2 = 0;

	while (!modelspcFile.eof()) {
        getline(modelspcFile, testStr);
        if (testStr.compare("Number of points in a/tan b distribution") == 0) {
            modelspcFile >> maxA1;
            if (maxA2 < maxA1) {
                maxA2 = maxA1;
            }
        } else if  (testStr.compare("Following identifies the reaches to be used as response time series and output time series") == 0) {
            modelspcFile >> Neq0 >> Nout0;
        }
        if (testStr.compare("The number of points in the overland flow distance distribution") == 0) {
            modelspcFile >> maxC1;
            if (maxC2 < maxC1) {
                maxC2 = maxC1;
            }
        }
	}
	maxA = maxA2;
	maxC = maxC2;

	maxResponse = max(Neq0, Nout0);
	modelspcFile.close();

	//------------ streamflow_calibration.dat -------------------
	// How many flow gauges?
	ifstream streamflow_calibrationFile("streamflow_calibration.dat");	// Fortran unit 2 (again)
		if (!streamflow_calibrationFile.is_open()) {
			cerr << "Failed to open streamflow_calibration.dat\n";
			exit(EXIT_FAILURE);
		}
	for (i = 1; i <= 3; i++) {	// Discard first three lines.
		getline(streamflow_calibrationFile, inLine, '\n');
	}
	streamflow_calibrationFile >> verno >> maxFgauge;
	streamflow_calibrationFile.close();
	maxResponse = max(maxResponse, maxFgauge);

	// How many temperature gauges?
	ifstream cliparFile("clipar.dat");	// Fortran unit 2 (again)
		if (!cliparFile.is_open()) {
			cerr << "Failed to open clipar.dat\n";
			exit(EXIT_FAILURE);
		}
	for (i = 1; i <= 2; i++) {	// Discard first two lines.
		getline(cliparFile, inLine, '\n');
	}
	cliparFile >> maxTGauge;
	cliparFile.close();

	//------------ rchareas.txt -------------------
	// Added to determine maxRchAreas
	// which is the no of lines in rchareas.txt
	exist = false;
	if ((dp = opendir(".")) == NULL) {
		cerr << "Can't open directory\n";
		exit(EXIT_FAILURE);
	}
	while ((dirp = readdir(dp)) != NULL) {
		testStr = dirp->d_name;
		if (testStr == "rchareas.txt") {
			exist = true;
			closedir(dp);
			break;
		}
	}
	if (exist == false) {
		cerr << " No 'rchareas.txt' file found\n";
		exit(EXIT_FAILURE);
	}
	if (!rchareasFile.is_open()) {
		rchareasFile.open("rchareas.txt");	// Fortran unit 9
		if (!rchareasFile.is_open()) {
			cerr << "Failed to open rchareas.txt\n";
			exit(EXIT_FAILURE);
		}
	} else {	// rewind
		rchareasFile.clear();
		rchareasFile.seekg(0);
	}
	lineCount = 0;
	while (getline(rchareasFile, testStr, '\n')) {
		lineCount++;
	}
	maxRchAreas = lineCount;
	rchareasFile.close();
	//-------------------------------------------

	//------------ rain.dat -------------------
	// Added the code relating to rain.dat to cope with when
	// there are more sites in rain.dat than are required by the model.
	string rainFileName;
#ifdef WRIA1
    rainFileName = "rain_allWRIA1.dat";
#else
    rainFileName = "rain.dat";
#endif
	exist = false;
	if ((dp = opendir(".")) == NULL) {
		cerr << "Can't open directory\n";
		exit(EXIT_FAILURE);
	}
	while ((dirp = readdir(dp)) != NULL) {
		testStr = dirp->d_name;
		if (testStr == rainFileName) {
			exist = true;
			closedir(dp);
			break;
		}
	}
	if (exist == false) {
		cerr << " No " << rainFileName << " file found\n";
		exit(EXIT_FAILURE);
	}
	if (!rainFile.is_open()) {
		rainFile.open(rainFileName);	// Fortran unit 9 (again)
		if (!rainFile.is_open()) {
			cerr << "Failed to open " << rainFileName << '\n';
			exit(EXIT_FAILURE);
		}
	} else {	// rewind
		rainFile.clear();
		rainFile.seekg(0);
	}
	for (i = 1; i <= 2; i++) {	// Discard first two lines.
		getline(rainFile, inLine, '\n');
	}
	rainFile >> verno;
	// rewind
	rainFile.clear();
	rainFile.seekg(0);
	if (verno == "Ver1" || verno == "ver1" || verno == "Ver2" || verno == "ver2") {
		for (i = 1; i <= 2; i++) {	// Discard first two lines.
			getline(rainFile, inLine, '\n');
		}
		rainFile >> verno >> Ntri;
	}
	maxGauge = max(maxGauge, Ntri);
	rainFile.close();
	//-------------------------------------------

	// Error check on iex and nrx.
	if (Nrx < maxInt) {
		cerr << "Number of observations in longest response record has to be at least " << maxInt << endl;
		exit(EXIT_FAILURE);
	}
	if (iex < maxResponse) {
		cerr << "Maximum number of responses used in fitting of model has to be at least" << maxResponse << endl;
		exit(EXIT_FAILURE);
	}

	//------------ lakes.dat -------------------
	// Read values from data files
	exist = false;
	if ((dp = opendir(".")) == NULL) {
		cerr << "Can't open directory\n";
		exit(EXIT_FAILURE);
	}
	while ((dirp = readdir(dp)) != NULL) {
		testStr = dirp->d_name;
		if (testStr == "lakes.dat") {
			exist = true;
			closedir(dp);
			break;
		}
	}
	if (exist) {
		ifstream lakesFile("lakes.dat");	// Fortran unit 5
		if (!lakesFile.is_open()) {
			cerr << "Failed to open lakes.dat\n";
			exit(EXIT_FAILURE);
		}
		getline(lakesFile, inLine, '\n');
		lakesFile >> max_lakes;
		// rewind
		lakesFile.clear();
		lakesFile.seekg(0);
		Rpairs1 = 0;
		while (getline(lakesFile, testStr, '\n')) {
			if (testStr.find("Give the number of rating \"pairs\"") < testStr.npos) {
				lakesFile >> Rpairs1;
				if (Rpairs2 < Rpairs1) {
					Rpairs2 = Rpairs1;
				}
			}
		}
		lakesFile.close();
		max_lheads = Rpairs2;
	} else {
		max_lakes = 1;
		max_lheads = 1;
	}

	maxSites = max(maxGauge, max(maxResponse, 3));
	// ************************************************************************************************

	// Allocate dynamic arrays
    ArrayXXi ClinkR(3,maxChn);
	ArrayXXi  linkR(4,maxChn);
	rel  = new double[maxRchAreas];
	for (n = 0; n < maxRchAreas; n++) {
		rel[n] = 0.0;
	}
	ishift0 = new int[maxRchAreas];
	bTmax = new double*[maxSlp];
	bTmin = new double*[maxSlp];
	bTdew = new double*[maxSlp];
	ArrayXXd   wrg(maxSlp,maxGauge);
	ArrayXXd  wrg1(maxSlp,maxGauge);
	ArrayXXi   lrg(maxSlp,maxGauge);
	ArrayXXd bRain(maxSlp,maxInt);
	for (i = 0; i < maxSlp; i++) {
		bTmax[i] = new double[maxInt];
		bTmin[i] = new double[maxInt];
		bTdew[i] = new double[maxInt];
	}
	wind2m = new double[maxInt];
	bXlat  = new double[maxSlp];
	bXlon  = new double[maxSlp];
	flow = new double*[maxResponse];
	for (i = 0; i < maxResponse; i++) {
		flow[i] = new double[maxInt];
	}
	ArrayXd temper(maxInt);
	dewp   = new double[maxInt];
	trange = new double[maxInt];
	Nts = new int[maxSlp];
	ArrayXXi linkS(2,maxSlp);
	Si  = new double*[Nsi];
	Sid = new double*[Nsi];
	for (i = 0; i < Nsi; i++) {
		Si[i]  = new double[maxSlp];
		Sid[i] = new double[maxSlp];
	}
	ArrayXXd Sp(Nsp,maxSlp);
	ArrayXXd Spd(Nsp,maxSlp);
	ArrayXXd bp(num_basinpars,maxSlp);
	Nka   = new int[maxSlp];
	Nd    = new int[maxSlp];
	iBout = new int[maxSlp];
	atb = new double*[maxA];
	pka = new double*[maxA];
	for (i = 0; i < maxA; i++) {
		atb[i] = new double[maxSlp];
		pka[i] = new double[maxSlp];
	}
	tl = new double[maxSlp];
	pd = new double*[maxC];
	cl = new double*[maxC];
	for (i = 0; i < maxC; i++) {
		pd[i] = new double[maxSlp];
		cl[i] = new double[maxSlp];
	}

	ll    = new int[maxChn];
	Ntr   = new int[maxChn];
	llOut = new int[maxChn];
	Rp  = new double*[Nrp];
	Rpd = new double*[Nrp];
	for (i = 0; i < Nrp; i++) {
		Rp[i]  = new double[maxChn];
		Rpd[i] = new double[maxChn];
	}
	ArrayXi kllout(maxResponse);
	Qmap   = new int[maxResponse];
	Xlat   = new double[maxTGauge];
	Xlong  = new double[maxTGauge];
	elevtg = new double[maxTGauge];
	temper_id = new int[maxTGauge];
	dtBar  = new double*[12];
	bdtBar = new double*[12];
	for (i = 0; i < 12; i++) {
		dtBar[i]  = new double[maxTGauge];
		bdtBar[i] = new double[maxSlp];
	}

	// ***********************************************************************

	//Define model identification string and number of parameters
	modelId = "TOP model";
	cliParam(Xlat, Xlong, stdlon, elevtg, dtBar, ns_temper, temper_id);
	iret = 0;
	mdData(model1::nGauge, model1::Ns, model1::Nrch, Nka, tl, atb, pka, Nd, cl, pd, units, ll, Ntr, Nts, linkS,
		linkR, Sid, Spd, Rpd, iret, model4::pmap, Npar, lrg, wrg, iex, llOut, Neq, model4::Nout, model4::nBout, iBout,
		nRchSav, Qmap, rel, relFlag, constr::minSp, constr::maxSp, constr::minSi, constr::maxSi,
		constr::minRp, constr::maxRp, constr::limitC, ClinkR, kllout, ishift0, maxInt,
		maxGauge, maxSlp, maxResponse, maxA, maxC, maxChn, maxRchAreas, maxSites, bp, bXlat, bXlon, wrg1);
	// Added the next bit for Neq = 0;
	neq_temp = Neq;
	hyData(model2::sDate, model2::sHour, model2::interval, model2::m, model2::mi, model2::mps, model2::mpe, model1::nGauge, Neq,
		bRain, flow, iret, temper, dewp, trange, dtBar, model1::Ns, wrg, lrg, elevtg, bTmax, bTmin, bTdew, bdtBar, Spd,
		maxGauge, maxInt, maxSites, maxResponse, maxTGauge, wind2m, wrg1, idebugoutput, idebugbasin, idebugcase);
	if ( iret != 0 ) {
		cerr << " **** an error has occurred on data input (hyData())\n";
		exit(EXIT_FAILURE);
	}
	// Rearrange the FLOW array so that the "columns" match the reach order
	// read in from MODELSPC.DAT using Qmap.
	// Qmap(i) is the column in runoff.dat corresponding to reach j
	for (i = 1; i <= neq_temp; i++) {
		for (j = 1; j <= neq_temp; j++) {
			if (Qmap[j-1] != j) {
				if (Qmap[j-1] == i) { // need to change order of ith and jth  columns
					// Interchange columns i and Qmap(i) and flag which columns are affected using imap
					for (k = 1; k <= model2::m; k++) {
						dummyq = flow[i-1][k-1];
						flow[i-1][k-1] = flow[Qmap[i-1]-1][k-1];
						flow[Qmap[i-1]-1][k-1] = dummyq;
					}
					Qmap[j-1] = Qmap[i-1];
					Qmap[i-1] = i;
						goto L150;
				}
			}
		}
L150: ;
	}

	// At this point we have all the data needed to sort out the how to update
	// the ZBAR0s for each measured sub-basin
	if (neq_temp > maxResponse) {
		cerr << " ****ERROR, maxResponse needs to be set to: " << neq_temp << '\n';
		exit(EXIT_FAILURE);
	}

	td8micsec(model2::sDate, model2::sHour, stim8);
	model2::stim = stim8;
	//       Neq = 1
	//       initT(Neq) = 1
	//       iend(Neq) = model2::m

	// Read flow data and precip data

	for (i = 1; i <= max(Neq, 1); i++) {
		initT[i-1] = 1;
		iend[i-1] = model2::m;
		if ( iend[i-1] <= Nrx ) {
			for (j = initT[i-1]; j <= iend[i-1]; j++) {
				Qact(j-1,i-1) = flow[i-1][j-1];
				//  p[j-1] = PREC(1,j)
			}
		} else {
            cerr << " **** You are trying to read too much data\n";
			exit(EXIT_FAILURE);
		}
	}

	// Assign actual time to each observation

	//        for (i = initt(1), iend(1)
	for (i = 1; i <= model2::m; i++) {
           acTime[i-1] = double(i);
	}

	// All these variables have been deallocated before they enter MODEL
	// Deallocate dynamic arrays
	for (i = 0; i < maxSlp; i++) {
		delete [] bTmax[i];
		delete [] bTmin[i];
		delete [] bTdew[i];
	}
	delete [] bTmax;
	delete [] bTmin;
	delete [] bTdew;

	delete [] wind2m;
	for (i = 0; i < maxResponse; i++) {
		delete [] flow[i];
	}
	delete [] flow;
	for (i = 0; i < Nsi; i++) {
		delete [] Si[i];
		delete [] Sid[i];
	}
	delete [] Si;
	delete [] Sid;
	for (i = 0; i < Nrp; i++) {
		delete [] Rp[i];
		delete [] Rpd[i];
	}
	delete [] Rp;
	delete [] Rpd;
	for (i = 0; i < maxA; i++) {
		delete [] atb[i];
		delete [] pka[i];
	}
	delete [] atb;
	delete [] pka;
	for (i = 0; i < maxC; i++) {
		delete [] pd[i];
		delete [] cl[i];
	}
	delete [] pd;
	delete [] cl;

	delete [] tl;
	delete [] Nka;
	delete [] Nd;
	delete [] iBout;
	delete [] Qmap;
	delete [] Nts;
    //delete [] temper;
    delete [] dewp;
    delete [] trange;
	delete [] rel;
	delete [] ishift0;
	delete [] ll;
	delete [] Ntr;
	delete [] llOut;

	// Following are variable arrays left undeleted in the fortran code
	// which are then reallocated in subroutine model().
	delete [] bXlat;
	delete [] bXlon;
	delete [] Xlat;
	delete [] Xlong;
	delete [] temper_id;
	delete [] elevtg;
	for (i = 0; i < 12; i++) {
		delete [] dtBar[i];
		delete [] bdtBar[i];
	}
	delete [] dtBar;
	delete [] bdtBar;
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving inputT(" << ncalls << ")" << "\n\n";
	}
	ncalls++;
#endif

	return 0;
}


// ************************************************************************72
//
//     Subroutine model designed to be called to compute response given a
//     time set and set of parameter values.  This is called once for each time step
//      Arguments
//      it - integer time step number
//      iflag - initialization flag set to 1 when model called first in an interation
//       0 otherwise.  This logic is circumvented here.  The model is run
//       only once each time step.  When iflag = 0 all this routine does is
//       copy computed results into the output array qfit
//      iopt - derivative flag set to 1 when nlfit requires response derivatives.
//       This is not possible here so an error is written
//      prt, a print indicatior 'y' or 'n' to control I/O
//      Neq - number of response variables
//      Npar - number of model parameters
//      npx - dimensioned size of number of parameters
//      iex - dimensioned size of number of response variables
//      nfor - fortran unit for I/O  - not used
//      qfit - Vector output flows - length iex, giving each response
//      par - vector of model parameters
//      mfit - flag specifying fit status indicating when a response
//        should be given - not used here - we always return flow
//        and never return derivative.
//      ifit - vector specifying fit status of each parameter.
//         Could apply parameter mapping only to these parameters - taking
//         others from modelspc.dat file
//      ibeale - flag intended to control parameter constraint checking
//         No param constraint checking yet implemented.
// **************************************************************************

int Model(const int it, const int iflag, const int iopt, const char prt, int &Neq, int &Npar, const int npx, const int iEx,
	const int nfor, double qfit[], const double par[], const int mfit[], const int ifit, const int ibeale)
{
    int n;
	// Note that Nchn is the same as Nrch by way of fortran common blocks within toplim_v7.f but not including calv46sn_v7.f:9:
	// Passed variables
    ArrayXXi  linkR(4,maxChn);
    ArrayXXi ClinkR(3,maxChn);
	int *ll, *Ntr, *llOut;
	int *Nts;
	int *Qmap;
	int iret;
	double **Si, **Rp;
    ArrayXXd Sp(Nsp,maxSlp);
    ArrayXXd Spd(Nsp,maxSlp);
    ArrayXXd bp(num_basinpars,maxSlp);
    ArrayXXd bRain(maxSlp,maxInt);

	//     ***********
	double **Sid;
	double **Rpd;
	//       real*8 barea(maxchn)  !  dgt 5/27/12  seems to be a problem in that maxchn is dynamic and barea has not been passed in
	double *bArea;
	//     the arrays sid, spd and rpd record the original values  dgt
	int relflag;
	int ngut, jg1, jg, irch, jr, jr1, ipos;
	double *rel;
	int *Nka, *Nd, *iBout;
	int *lOut;
	int *iReach, *knt;
	int neq_temp;

	double **atb, **pka, *tl;
	double **pd, **cl;
	double units;

	double dummyq;

	int *ishift0;
	int ipsub, ipatb;

	//     time series info
	//      integer m, mi
    bool reinit, modwrt;
    bool limitS; // limits is local, limitc global for constraints
	//     outputs
	double **q4;  // made allocatable as was not passed in
	double qmod, mse, mae;
	//       real*8 q5(maxint+1,maxchn)  ! rpi 31/7/2002 for Neq=0 + maxchn
	bool ok;
	int *ntdh;

	int ndump;
	double **bTmax;
	double **bTmin;
	double *wind2m;
	double **bTdew;
	double *bXlat;
	double *bXlon;
	double **flow;
	//     etvariables  dgt
	double *Xlat, *Xlong;
	double stdlon;
	double *elevtg, **dtBar, **bdtBar;
	int *temper_id;
	ArrayXd temper(maxInt);
	double *dewp, *trange;
	double *dFlow, *suma;

	string fname;
	//     declarations for lake varibales
	int ns_temper;
	int *lzero, *lk_line;
	int *lake_beach_slps, *num_rat_vals;
	int **lheads, **loflows;
	double *lake_areas, *ini_levels;
	int j, i, k;
	int idebugoutput, idebugbasin, idebugcase;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> Model(" << ncalls << ")" << std::endl;
    }
    caller = "Model";
#endif
	//  Allocate barea
	//  write(6,*)maxchn,maxint,maxchn,maxslp
	bArea = new double[maxChn];
	q4 = new double*[maxInt];
	for (i = 0; i < maxInt; i++) {
		q4[i] = new double[maxChn];
	}
	ntdh = new int[maxSlp];
	//     ********************************************
	//     To solve the problem created by common block
	iret = 0;

    if (model_master::iDoIt == 0) {
		model_master::iDoIt = 1;

		//       All these variables have been dynamically allocated to solve the problem
		//       ALLOCATE DYNAMIC ARRAYS
		rel  = new double[maxRchAreas];
		ishift0 = new int[maxRchAreas];
		bTmax = new double*[maxSlp];
		bTmin = new double*[maxSlp];
		bTdew = new double*[maxSlp];
		ArrayXXd  wrg(maxSlp,maxGauge);
		ArrayXXd wrg1(maxSlp,maxGauge);	//  Interpolation weights
		ArrayXXi     lrg(maxSlp,maxGauge);
		Nts = new int[maxSlp];
		for (i = 0; i < maxSlp; i++) {
			bTmax[i] = new double[maxInt];
			bTmin[i] = new double[maxInt];
			bTdew[i] = new double[maxInt];
		}
		wind2m = new double[maxInt];
		bXlat  = new double[maxSlp];
		bXlon  = new double[maxSlp];
		flow = new double*[maxResponse];
		for (i = 0; i < maxResponse; i++) {
			flow[i] = new double[maxInt];
		}
		elevtg = new double[maxTGauge];
		temper_id = new int[maxTGauge];
		dtBar  = new double*[12];
		bdtBar = new double*[12];
		for (i = 0; i < 12; i++) {
			dtBar[i] = new double[maxTGauge];
			bdtBar[i] = new double[maxSlp];
		}
		dewp   = new double[maxInt];
		trange = new double[maxInt];
		ArrayXXi linkS(2,maxSlp);
		Si  = new double*[Nsi];
		Sid = new double*[Nsi];
		for (i = 0; i < Nsi; i++) {
			Si[i]  = new double[maxSlp];
			Sid[i] = new double[maxSlp];
		}
		Nka   = new int[maxSlp];
		Nd    = new int[maxSlp];
		iBout = new int[maxSlp];
		atb = new double*[maxA];
		pka = new double*[maxA];
		for (i = 0; i < maxA; i++) {
			atb[i] = new double[maxSlp];
			pka[i] = new double[maxSlp];
		}
		tl = new double[maxSlp];
		pd = new double*[maxC];
		cl = new double*[maxC];
		for (i = 0; i < maxC; i++) {
			pd[i] = new double[maxSlp];
			cl[i] = new double[maxSlp];
		}
		ll    = new int[maxChn];
		Ntr   = new int[maxChn];
		llOut = new int[maxChn];
		Rp  = new double*[Nrp];
		Rpd = new double*[Nrp];
		for (i = 0; i < Nrp; i++) {
			Rp[i]  = new double[maxChn];
			Rpd[i] = new double[maxChn];
		}
		ArrayXi kllout(maxResponse);
		Qmap   = new int[maxResponse];
		lOut   = new int[maxResponse];
		iReach = new int[maxResponse];
		knt    = new int[maxResponse];
		dFlow  = new double[maxResponse];
		suma   = new double[maxResponse];
		Xlat   = new double[maxTGauge];
		Xlong  = new double[maxTGauge];

		for (n = 0; n < maxResponse; ++n) {
            Qmap[n]   = 0;
            lOut[n]   = 0;
            iReach[n] = 0;
            knt[n]    = 0;
            dFlow[n]  = 0.0;
            suma[n]   = 0.0;
        }
        for (n = 0; n < maxTGauge; ++n) {
            Xlat[n]   = 0.0;
            Xlong[n]  = 0.0;
        }
		cliParam(Xlat, Xlong, stdlon, elevtg, dtBar, ns_temper, temper_id);
		mdData(model1::nGauge, model1::Ns, model1::Nrch, Nka, tl, atb, pka, Nd, cl, pd, units, ll, Ntr,
			Nts, linkS, linkR, Sid, Spd, Rpd, iret, model4::pmap, Npar, lrg, wrg, iex, llOut,
			Neq, model4::Nout, model4::nBout, iBout, nRchSav, Qmap, rel, relflag, constr::minSp, constr::maxSp,
			constr::minSi, constr::maxSi, constr::minRp, constr::maxRp, limitS, ClinkR, kllout, ishift0, maxInt,
			maxGauge, maxSlp, maxResponse, maxA, maxC, maxChn, maxRchAreas,
			maxSites, bp, bXlat, bXlon, wrg1);
		max_of_shifts = 0;
		for (i = 1; i <= maxRchAreas; i++) {
			max_of_shifts = max(max_of_shifts, ishift0[i-1]);
		}
        neq_temp = Neq;
		hyData(model2::sDate, model2::sHour, model2::interval, model2::m, model2::mi, model2::mps, model2::mpe, model1::nGauge, Neq,
			bRain, flow, iret, temper, dewp, trange, dtBar, model1::Ns, wrg, lrg, elevtg, bTmax, bTmin, bTdew, bdtBar, Spd,
			maxGauge, maxInt, maxSites, maxResponse, maxTGauge, wind2m, wrg1, idebugoutput, idebugbasin,  idebugcase);

		if (iret != 0) {
			cerr << " **** An error has occurred on data input\n";
			exit(EXIT_FAILURE);
		}
		//       Rearrange the FLOW array so that the "columns" match the reach order
		//       read in from MODELSPC.DAT using Qmap.
		//       Qmap(i) is the column in runoff.dat corresponding to reach j
        for (i = 1; i <= neq_temp; i++) {
			for (j = 1; j <= neq_temp; j++) {
				if (Qmap[j-1] != j) {
					if (Qmap[j-1] == i) { // need to change order of ith and jth  columns
						// Interchange columns i and Qmap(i) and flag which columns are affected using imap
						for (k = 1; k <= model2::m; k++) {
							dummyq = flow[i-1][k-1];
							flow[i-1][k-1] = flow[Qmap[i-1]-1][k-1];
							flow[Qmap[i-1]-1][k-1] = dummyq;
						}
						Qmap[j-1] = Qmap[i-1];
						Qmap[i-1] = i;
						goto L150;
					}
				}
			}
L150: ;
        }

        if (nooksack == 0) {
			//        Calculate sum of catchment areas above each reach
			ngut = model1::Ns;
			for (jg1 = 1; jg1 <= ngut; jg1++) {
				jg = ll[jg1-1];
				//  inserted the following line to ensure that area(jg) defined
				bArea[jg-1] = 0.0;
				for (j = 1; j <= linkR(3,jg-1); j++) {
					bArea[jg-1] += fabs(Spd(0,linkR(j,jg-1)-1));
				}
			}
			irch = model1::Nrch - ngut;
			for (jr1 = 1; jr1 <= irch; jr1++) {
				jr = ll[jr1+ngut-1];
				bArea[jr-1] = 0.0;
				for (j = 1; j <= linkR(3,jr-1); j++) {
					ipos = iposn(maxChn, linkR, linkR(j,jr-1));
					bArea[jr-1] += bArea[ipos-1];
				}
			}
			// At this point we have all the data needed to sort out the how to update
			// the ZBAR0s for each measured sub-basin


			set_cor_data(ClinkR, Neq, model4::Nout, nRchSav, model1::Ns, kllout, flow, Spd,
				bArea, ll, llOut, dFlow, suma, knt, iReach, lOut, maxInt, maxSlp, maxChn,maxResponse);
			//  ************************************************
        }
		// ALLOCATE DYNAMIC ARRAYS
		Array<int,Dynamic,1> lake_reach(max_lakes);
		lzero           = new int[max_lakes];
		lk_line         = new int[max_lakes];
		lake_beach_slps = new int[max_lakes];
		num_rat_vals    = new int[max_lakes];
		lheads          = new int*[max_lakes];
		loflows         = new int*[max_lakes];
		for (i = 0; i < max_lakes; i++) {
			lheads[i]   = new int[max_lheads];
			loflows[i]  = new int[max_lheads];
		}
		lake_areas = new double[max_lakes];
		ini_levels = new double[max_lakes];
		// ***********************************************************************

		// See if there is a Lakes.dat file to be processed
		// removed to make UNIX version work
		lakes1::nlakes = 0; // Initialise in case there is no LAKES.DAT
		// Read in lake data
		read_lakes(lake_reach, lzero, lake_areas, lake_beach_slps, lk_line, num_rat_vals,
			lheads, loflows, max_lakes, max_lheads);
		// RPI 17/05/2002 added in the initial lake levels stuff
		// See if there is a initial lake levels file Lakes_levels.dat file to be processed
		// Initialise initial lake levels to -1 in case none are read in
		if (lakes1::nlakes > 0) {
			for (i = 1; i <= lakes1::nlakes; i++) {
              ini_levels[i-1] = -1.0;
			}
		}
		//       Read in initial lake level data
		read_lakes_levels(lake_reach, ini_levels, max_lakes);
	} //once

	//      call SET_COR_DATA(clinkr,Neq,nout,nrchsav,model1::Ns,kllout,flow,spd,
	//     *  barea,ll,llout,maxInt,MAXSLP,MAXCHN,MAXRESPONSE)
	// ************************************************************************


	// Labelled common passed from subroutine inputt

	//CC RAW 7/4/97 changed p(9000) to p(maxInt)
	//      double precision p(maxInt)
	//      common /abc_model/ p
	//      save
	//-----
	// p(it) is rainfall input for current time step it
	// dq is a vector of runoff derivatives wrt a,b,c and so
	// ds is a vector of storage derivatives wrt a,b,c and s0
	//
	// If this is the first call in an iteration initialize groundwater storage
	// and its derivatives and give parameters meaningful names to make code
	// more readable
	//
	// so is antecedent storage

	// Note use of the save statement to ensure the assignments are not lost

	//-----
	// If Beale nonlinearity measure is being computed check for
	// infeasible parameters
	//
    if (iflag == 1) {    //  Only do work first time through
		//   copy parameters from defaults  ! Inserted by DGT at to work with interpretation
		//     of parameters varied in calibration as multiplication factors.
		// Keep track of calling parameters in toperror.txt

		for (i = 1; i <= model1::Ns; i++) {
			for (j = 1; j <= Nsi; j++) {
				Si[j-1][i-1] = Sid[j-1][i-1];
			}
			for (j = 1; j <= Nsp; j++) {
				Sp(j-1,i-1) = Spd(j-1,i-1);
			}
		}
		for (i = 1; i <= model1::Nrch; i++) {
			for (j = 1; j <= Nrp; j++) {
				Rp[j-1][i-1] = Rpd[j-1][i-1];
			}
		}

		//  Now treat calibration parameters as multiplication factors without worrying
		//  about losing information
		// Unscramble the parameter map
		for (i = 1; i <= Npar; i++) {
			if (model4::pmap(2,i-1) <= 0) {
				k = model1::Ns;
				if (model4::pmap(0,i-1) == 2)
					k = model1::Nrch;
				for (j = 1; j <= k; j++) {
					// dealing with a sub-basin's properties
					if (model4::pmap(0,i-1) == 1)
						Sp(model4::pmap(1,i-1)-1,j-1) = par[model4::pmap(3,i-1)-1]*Sp(model4::pmap(1,i-1)-1,j-1);
					// dealing with a river reach
					if (model4::pmap(0,i-1) == 2)
						Rp[model4::pmap(1,i-1)-1][j-1] = par[model4::pmap(3,i-1)-1]*Rp[model4::pmap(1,i-1)-1][j-1];
					// dealing with an initial condition
					if (model4::pmap(0,i-1) == 3)
						Si[model4::pmap(1,i-1)-1][j-1] = par[model4::pmap(3,i-1)-1]*Si[model4::pmap(1,i-1)-1][j-1];
				}
			} else {
				// dealing with a sub-basin's properties
				if (model4::pmap(0,i-1) == 1)
					Sp(model4::pmap(1,i-1)-1,model4::pmap(2,i-1)-1) = par[model4::pmap(3,i-1)-1]*Sp(model4::pmap(1,i-1)-1,model4::pmap(2,i-1)-1);
				// dealing with a river reach
				if (model4::pmap(0,i-1) == 2)
					Rp[model4::pmap(1,i-1)-1][model4::pmap(2,i-1)-1] = par[model4::pmap(3,i-1)-1]*Rp[model4::pmap(1,i-1)-1][model4::pmap(2,i-1)-1];
				// dealing with an initial condition
				if (model4::pmap(0,i-1) == 3)
					Si[model4::pmap(1,i-1)-1][model4::pmap(2,i-1)-1] = par[model4::pmap(3,i-1)-1]*Si[model4::pmap(1,i-1)-1][model4::pmap(2,i-1)-1];
			}
		}
		// check if constraints apply and if so check them
		if (constr::limitC) {
			for (i = 1; i <= model1::Ns; i++) {
				for (j = 1; j <= Nsi; j++) {
					Si[j-1][i-1] = max(constr::minSi[j-1], min(Si[j-1][i-1], constr::maxSi[j-1]));
				}
				for (j = 1; j <= Nsp; j++) {
					Sp(j-1,i-1) = max(constr::minSp[j-1], min(Sp(j-1,i-1), constr::maxSp[j-1]));
				}
			}
			for (i = 1; i <= model1::Nrch; i++) {
				for (j = 1; j <= Nrp; j++) {
					Rp[j-1][i-1] = max(constr::minRp[j-1], min(Rp[j-1][i-1], constr::maxRp[j-1]));
				}
			}
		}

		// Calculate the initial value of ZBAR0 for each sub-basin. A uniform
		// value is assigned to all sub-basins in each measured sub-basin

		ndump = 0;
		setup_zbar0s(model1::Ns, Si, Sp, tl, model2::interval, lOut, iReach, dFlow, suma, knt, maxResponse, maxSlp);

		//      if (ibeale == 1) then
		//               write(*,*) ' **** Constraint checking still to be done'
		//       for (i = 1,11
		//           if ( par[i-1] .le. 0.) then
		//           ibeale = 2
		//           return
		//            end if
		//       end do
		//      end if

		//-----
		// Print option which is enabled in output mode in nlfit and in
		// predict.  Note that this is optional

 	if (prt != 'Y') {
			reinit = true;
			ipsub = model1::Ns/2 + 1;
			ipatb = Nka[model1::Ns-1]/2 + 1;
			modwrt = false;
			calcts( Si,            Sp,           Rp,               linkR,         ll,
					model1::Ns,                  Nka,              tl,            atb,
					pka,           Nd,           cl,               pd,            units,
					ipsub,         ipatb,        reinit,           modwrt,        model2::stim,
					bRain,        model2::interval, model2::m,     model2::mi,
					model2::mps,   model2::mpe,  ok,               Xlat,          Xlong,
					stdlon,        elevtg,       bdtBar,           model2::sDate, model2::sHour,
					temper,        dewp,         trange,           Neq,           model4::Nout,
					model4::nBout, iBout,        wind2m,           bTmin,         bTmax,
					bTdew,         bXlat,        bXlon,            ntdh,          ndump,
					maxInt,        maxSlp,       maxA,             maxC,          maxChn,
					idebugoutput,  idebugbasin,  idebugcase);
		} else {	// l1
			if (idebugoutput >= 1) {
				lunmodFile.open("results/topsbd_v8.cpp.txt");
				luntopFile.open("results/topreachd_v7.cpp.txt");
				lundatFile.open("results/detail_v7.cpp.txt");
				lunpFile.open("results/topinfo_v7.cpp.txt");
				lunpFile << "Topmodel output\n";
				lunpFile << " Simulation starts at" << dec << setw(9) << model2::sDate;
				lunpFile << dec << setw(8) << model2::sHour << " uses a time interval of ";
				lunpFile << dec << setw(6) << model2::interval << " seconds\n";
				lunpFile << " and contains " << dec << setw(8) << model2::m << " steps.";
				lunpFile << " Detailed printout begins at interval " << dec << setw(8) << model2::mps;
				lunpFile << " and ends at " << dec << setw(8) << model2::mpe << '\n';
			}
			if (constr::limitC) {   // tidying up
                cerr << " Limits are being applied: Min Value Max\n";
                for (i = 1; i <= Nsp; i++) {
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::minSp[i-1];
					cerr << " " << fixed << setw(12) << setprecision(4) << Sp(i-1,0);
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::maxSp[i-1] << '\n';
                }
                for (i = 1; i <= Nsi; i++) {
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::minSi[i-1];
					cerr << " " << fixed << setw(12) << setprecision(4) << Si[i-1][0];
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::maxSi[i-1] << '\n';
                }
                for (i = 1; i <= Nrp; i++) {
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::minRp[i-1];
					cerr << " " << fixed << setw(12) << setprecision(4) << Rp[i-1][0];
					cerr << " " << fixed << setw(12) << setprecision(4) << constr::maxRp[i-1] << '\n';
                }
			}
			reinit = true;
			ipsub = model1::Ns/2 + 1;
			ipatb = Nka[model1::Ns-1]/2 + 1;
			modwrt = true;
			calcts( Si,            Sp,           Rp,               linkR,         ll,
					model1::Ns,                  Nka,              tl,            atb,
					pka,           Nd,           cl,               pd,            units,
					ipsub,         ipatb,        reinit,           modwrt,        model2::stim,
					bRain,        model2::interval, model2::m,     model2::mi,
					model2::mps,   model2::mpe,  ok,               Xlat,          Xlong,
					stdlon,        elevtg,       bdtBar,           model2::sDate, model2::sHour,
					temper,        dewp,         trange,           Neq,           model4::Nout,
					model4::nBout, iBout,        wind2m,           bTmin,         bTmax,
					bTdew,         bXlat,        bXlon,            ntdh,          ndump,
					maxInt,        maxSlp,       maxA,             maxC,          maxChn,
					idebugoutput,  idebugbasin,  idebugcase);
			// added two close statements. Without them, files get very big
			if (idebugoutput >= 1) {
				lundatFile.close();
				lunmodFile.close();
				luntopFile.close();
				lunpFile.close();
			}
		}	// l1

		//    Keep track of some stats in file toperror.txt
		for (j = 1; j <= max(1, Neq); j++) {
			mse = 0.0;
			mae = 0.0;
			for (i = 1; i <= model2::m; i++) {
				qmod = 0; // q4(i,ll(iabs(llout(j))))
				mse += (qmod - flow[j-1][i-1])*(qmod - flow[j-1][i-1]);
				mae += fabs(qmod - flow[j-1][i-1]);
			}
			mse = mse/double(model2::m);
			mae = mae/double(model2::m);
			cerr << " Mean Sq. Error and Mean Abs. Error " << mse << " " << mae << '\n';
		}


	} // flag
	//-----
	// If required compute model derivatives algebraically
	// Algebraic derivatives are OPTIONAL; you can use NLFIT's finite
	// difference option. The advantage is speed and better accuracy,
	// but coding is very error prone.

	if (iopt == 1)
		cerr << " ***** Warning - analytical derivatives unavailable\n";

	// Put observed responses into array qfit for return

	for (j = 1; j <= max(1, Neq); j++) {
		qfit[j-1] = 0; // q4(it,ll(iabs(llout(j)))) pass back variables to calibrate on
	}

	// On first call to model write heading for table
	//       if (iflag == 1) write(nfor,1)
	// Write rainfall with predicted and observed runoff and storage

	//        if (mfit(1) == -1.or.mfit(2) == -1) then
	//           write(nfor,2) it, p(it), (qfit(j),qa(j),j=1,2), 'Yes'
	//        else
	//           write(nfor,2) it, p(it), (qfit(j),qa(j),j=1,2), ' No'
	//        end if
	//      end if

	//...
	//1     format(' Time  Rainfall  Pred runoff  Obs runoff  Pred storage',
	//     &       '  Obs storage Censored'/
	//     &       ' -----------------------------------------------------',
	//     &       '----------------------')
	//2     format(' ',i4,f10.0,f13.0,f11.0,f14.0,f13.0,6x,a)
	//-----
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving Model(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

// ************************************************************************

int set_cor_data(const ArrayXXi &linkR, const int Neq, const int Nout, const int Nrch, const int Ns,
	const Array<int,Dynamic,1> &kllout, double **flow, const ArrayXXd &Spd, const double *bArea, const int *ll,
	const int *llOut, double *dFlow, double *suma, int *knt, int *iReach, int *lOut, const int maxInt,
	const int maxSlp, const int maxChn, const int maxResponse)
{
    double fArea[maxResponse], nflow_mean,  pArea[maxResponse];
    int idFlow[maxResponse];

	double area_max;
	double dflow_max, dflow_mean;
	int nReach, i;
	int i1, i2, i3, i6;
	int ii, ib, i4, i5, neq2;
	int i_area=-1, iprev, ix;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> set_cor_data()" << std::endl;
    }
	caller = "set_cor_data";
#endif
	// Calculate the number of measured sub-basins removing possible duplicates
	// caused by having both flows and levels - remember only have data for a
	// maximum of Neq basins

	model5::Neq1 = 0;
	for (i = 1; i <= Neq; i++) {
		if (kllout(i-1) > 0) {
			model5::Neq1++;
			// lout contains the locations in LL of the reaches that are measured
              lOut[model5::Neq1-1]  = kllout(i-1);
              fArea[model5::Neq1-1] = bArea[ll[llOut[i-1]-1]-1];
              pArea[model5::Neq1-1] = fArea[model5::Neq1-1];  // RPI 3/7/2002 - make a copy for later use
		}
	}
	for (i = 1; i <= model1::Ns; i++) {
		iReach[i-1] = 0;
	}

	for (i2 = model1::Nrch; i2 >= 1; i2--) { 	// Search backwards since most d/s reach is the first in
													// in the list and this reach does not appear elsewhere in the list
		for (i3 = 2; i3 <= 3; i3++) {
			ii = i2;
			if ((linkR(i3-1,i2-1) > model1::Ns) || (linkR(i3-1,i2-1) <= 0))
				goto L3; // Not a sub-basin
			if (iReach[linkR(i3-1,i2-1)-1] > 0)
				goto L3;  // This basin already assigned to a reach
			// have a basin not yet assigned to a reach - ib is the next reach d/s
L10: 		ib = linkR(0,ii-1);
			for (i4 = model5::Neq1; i4 >= 1; i4--) {	// Note reaches in modelspc.dat might need to have
														// u/s reaches to the right of d/s reaches
				if (ib == lOut[i4-1]) {
					iReach[linkR(i3-1,i2-1)-1] = lOut[i4-1];
					goto L3;
				}
			}
			//  else
			// find next connecting reach d/s
			for (i5 = model1::Nrch; i5 >= 1; i5--) { // Search backwards since most d/s reach is the first in
				// in the list and this reach does not appear elsewhere in the list
				for (i6 = 2; i6 <= 3; i6++) {
					if (linkR(i6-1,i5-1) == ib) {  // found next d/s reach
						ii = i5;
						goto L10;
					}
				}
			}
			// If we get here we have a problem
			cerr << " ***** Could not associate a measured reach with basin " << dec << setw(4) << linkR(i3-1,i2-1) << '\n';
			cerr << " Hint - this might occur if the most d/s model reach is below the most d/s measured reach\n";
			iReach[linkR(i3-1,i2-1)-1] = 1;
L3: ;
		}
	}

	for (i1 = 1; i1 <= model5::Neq1; i1++) {
		suma[i1-1] = 0;
		knt[i1-1] = 0;
		for (i4 = 1; i4 <= model1::Ns; i4++) {
			if (iReach[i4-1] == lOut[i1-1]) {
				suma[i1-1] = suma[i1-1] + Spd(0,i4-1); // area in mm*2
				knt[i1-1] = knt[i1-1] + 1;
			}
		}
	}
	neq2 = 0;
	// Take note of the flow from the biggest area,
	// so we can use it for ungauged subcats
	area_max = 0.0;
	for (i = 1; i <= Neq; i++) {
		if (kllout(i-1) > 0) {
			neq2++;
			// setup initial flows for the ZBAR0 calcs
			idFlow[neq2-1] = 1;
			dFlow[neq2-1] = flow[neq2-1][0]*fArea[neq2-1];  // flow in mm**3/interval
			if (fArea[neq2-1] > area_max) {
				area_max = fArea[neq2-1];
				i_area = neq2;
			}
		}
	}
	dflow_max = flow[i_area-1][0]*fArea[i_area-1];
	// now subtract flows from u/s from d/s flows, but do it only once!
	for (i1 = model1::Nrch; i1 >= 1; i1--) {
		if ((linkR(1,i1-1) <= model1::Ns) || (linkR(2,i1-1) <= model1::Ns)) { // have a basin
			iprev = 0;
			nReach = linkR(0,i1-1);  // This is the reach the basin flows into
			ix = i1 - 1;
L500:		if (ix == 0)
				goto L51;  //  Am at d/s node, go to end of loop
			for (i2 = ix; i2 >= 1; i2--) {
				if ((linkR(1,i2-1) == nReach) || (linkR(2,i2-1) == nReach))  {
					// found next d/s reach
					ix = i2 - 1;
					nReach = linkR(0,i2-1);
					for (i3 = 1; i3 <= model5::Neq1; i3++) {
						if (nReach == lOut[i3-1]) { // at a measured reach
							// is there an u/s reach feeding into this one?
							if (iprev > 0) {
								if (idFlow[iprev-1] > 0) {
									dFlow[i3-1] = dFlow[i3-1] - dFlow[iprev-1];
									// Parea is the partial area between 2 sites - farea = full area to a site
									pArea[i3-1] -= pArea[iprev-1];
									idFlow[iprev-1] = 0;
								}
								//  lout(i3,2)=nrchprev
							}
							iprev = i3;
							goto L500;
							// calculate tributary inflow her
						}
					}
				}
			}
		}
L51: ;
	}
	dflow_mean = 0.0;
	nflow_mean = 0;
	for (i = 1; i <= model5::Neq1; i++) {
		// RPI 3/7/2002 - weight flow with area contributing this flow
		dflow_mean += dFlow[i-1]*pArea[i-1];
		nflow_mean += pArea[i-1];
	}
	dflow_mean = dflow_mean/nflow_mean;
	if (dflow_mean <= 0)
		dflow_mean = dflow_max;
	dflow_max = dflow_mean;
	for (i = 1; i <= model5::Neq1; i++) {
		//  if (dflow(i).le.0.) dflow(i)=1.d-4  ! RPI 12/4/01 to avoid -ve tributary flows
		if (dFlow[i-1] <= 0.0) {
			// C RAW 28/6/01 1e-4 gave awful results - RPI tried dflow_max and found that
			// it was not much better and substituted the mean of the dflows instead on - 1/7/2002
			dFlow[i-1] = dflow_max*pArea[i-1]/nflow_mean;
			cerr << " Test - for response " << dec << setw(3) << i << " dflow set to " << fixed << setw(12) << setprecision(4) << dflow_max << '\n';
		}
		// dflow(i) = dabs(dflow(i))
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving set_cor_data(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}


// *****************************************************************************
//     SUBROUTINE SETUP_ZBAR0S
// *****************************************************************************
int setup_zbar0s(const int Ns, double **Si, const ArrayXXd &Sp, const double *tl, const long int dt,
	const int *lOut, const int *iReach, const double *dFlow, const double *sumA, const int *knt, const int maxResponse, const int maxSlp)
{
	double sumF[maxResponse], sumC[maxResponse];
	double area[maxSlp], k[maxSlp], f[maxSlp], Lambda[maxSlp];
	double **Si_copy;
	double sum_zbar0, av_zbar0;

	int temp_lout[maxResponse];
	int i, i1, i2, i3, i6, i7, ikeep_i6, inb, inz_zbar0;

	Si_copy = new double*[Nsi];
	for (i = 0; i < Nsi; i++) {
		Si_copy[i] = new double[maxSlp];
	}

	for (i1 = 1; i1 <= model5::Neq1; i1++) {
		sumF[i1-1] = 0;
		sumC[i1-1] = 0;
	}
	for (i2 = 1; i2 <= Ns; i2++) {
		Si_copy[1][i2-1] = Si[1][i2-1];
		Si[1][i2-1] = -100000.0;
		for (i3 = 1; i3 <= model5::Neq1; i3++) {
			if (iReach[i2-1] == lOut[i3-1]) {
				sumF[i3-1] += Sp(1,i2-1);
				sumC[i3-1] += (Sp(0,i2-1)/1.0e6)*Sp(2,i2-1)*exp(-tl[i2-1])/Sp(1,i2-1);
				if (sumC[i3-1] <= 0) {
					cerr << " **** Error - sumc(i) is negative - this probably means that f, or k needs to be constrained\n";
					exit(EXIT_FAILURE);
				}
				goto L5;
			}
		}
L5: ;
	}
	inz_zbar0 = 0;
	sum_zbar0 = 0.0;
	// The code below solves the third eqn from the bottom of the ZBAR0.DOC report
	// It does this for each of the neq1 basins for which there are measurements
	// The DFLOW variables are the incremental flow between 2 sites in VOLUME/time
	// inb is the number of sub-basins contributing flow to the i7'th site

	// save the whole of lout for later reference RPI 1/11/2001
	for (i = 0; i < maxResponse; i++) {
		temp_lout[i] = lOut[i];
	}
	for (i7 = 1; i7 <= model5::Neq1; i7++) {
		inb = 0;
		if (dFlow[i7-1] > 0) { // we can't cope with influent reaches
			// Pick out which basins contribute to the outlet reach defined by lout(7)
			for (i6 = 1; i6 <= Ns; i6++) {
				if (iReach[i6-1] == lOut[i7-1]) {
					// the basin in ireach(i6) contributes directly to outlet reach lout(i7)
					// WARNING the units of dflow need sorting out should be mm*3/dt
					Si[1][i6-1] = (-knt[i7-1]/sumF[i7-1])*log(dFlow[i7-1]/(((dt/3600.0)*1.0e+9)*sumC[i7-1]));
					inb++;
					area[inb-1] = Sp(0,i6-1);
					f[inb-1]    = Sp(1,i6-1);
					k[inb-1]    = Sp(2,i6-1);
					Lambda[inb-1] = tl[i6-1];
				}
			}
			// At this point we have all the relevant data from all the upstream basins
			// for solving for ZBAR0 for the reach defined by lout(i7)
			// Since we assume the same ZBAR0 applies to all the upstream basins we
			// only need to solve once for each set and then assign the solution to
			// each basin - RPI 1/11/2001
			ikeep_i6 = 0;
			for (i6 = 1; i6 <= Ns; i6++) {
				if (iReach[i6-1] == lOut[i7-1]) {
					if (temp_lout[i7-1] > 0) {
						solve_for_zbar0(Si[1][i6-1], area, f, k, Lambda, dFlow[i7-1], dt, inb);
						temp_lout[i7-1]= -temp_lout[i7-1];
						ikeep_i6 = i6;
					} else {
						Si[1][i6-1] = Si[1][ikeep_i6-1];
					}
					inz_zbar0++;
					sum_zbar0 += Si[1][i6-1];
				}
			}
		}
	}
	// guard against division by zero
	if (inz_zbar0 > 0) {
		av_zbar0 = sum_zbar0/inz_zbar0;
	} else {
		av_zbar0 = 1.0;
	}
	for (i6 = 1; i6 <= Ns; i6++) {
		if (Si[1][i6-1] < -99999.0) {
			// following line moved to inside new if statement below,
			// to cope with case when no flow data is available and so
			// ALL the zbar0 values would be undefined (so av_zbar0 undefined)
			//                                   Si(2,i6)=av_zbar0
			if (inz_zbar0 > 0) {
				Si[1][i6-1] = av_zbar0;
				cerr << " Test - zbar0 for basin " << dec << setw(4) << i6;
				cerr << " set to average of " << fixed << setw(12) << setprecision(4) << av_zbar0 << '\n';
			} else {
				Si[1][i6-1] = Si_copy[1][i6-1];
			}
		}
	}
    for (i = 0; i < Nsi; i++) {
    	delete [] Si_copy[i];
    }
    delete [] Si_copy;

	return 0;
}

