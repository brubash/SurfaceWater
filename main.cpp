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

//Main program to drive topmodel in standalone mode
#include "topnet.hh"

using namespace Eigen;

namespace model_master {
	int iDoIt;
}
//-----------------------
namespace mddata {
	int maxRno;
}
namespace lakes1 {
	int nlakes;
}
namespace ts_save {
	double Ts_old;
	double Tave_old;
	double Ts_Ave;
	double Tave_ave;
}
namespace rain {
	double rain_factor;
}
namespace ccompr2 {
	double t0;
}
namespace model1 {	// shared by inputT() and model()
	int nGauge;
	int Ns;
	int Nrch;
}
namespace model2 {	// shared by inputT() and model()
	int sDate;
	int sHour;
	long int interval;	// int
	int m;
	int mi;
	int mps;
	int mpe;
	int stim;
}
namespace model4 {	// shared by inputT() and model()
	ArrayXXi pmap(4,dpm);
	int Nout;
	int nBout;
}
namespace model5 {
	int Neq1;
}
namespace constr {	// shared by inputT() and model()
	double minSp[Nsp];
	double maxSp[Nsp];
	double minSi[Nsi];
	double maxSi[Nsi];
	double minRp[Nrp];
	double maxRp[Nrp];
	bool limitC;
}
//COMMON /LAKES1/ NLAKES
//COMMON /LAKES2/ LAKE
//COMMON /CCOMPR1/ RANGE,V0,V1,V2
//COMMON /CCOMPR2/ T0,T1,T2,SCOUNT
//COMMON /MODEL_MASTER/ IDOIT
//COMMON /MODEL1/ NGAUGE,NS,NRCH
//COMMON /MODEL2/ SDATE,SHOUR,INT,M,MI,MPS,MPE,STIM,AREA    (
//COMMON /MODEL4/ pmap,nout,nbout
//COMMON /CONSTR/ minsp,maxsp,minsi,maxsi,minrp,maxrp,limitc
//COMMON /MODEL5/ NEQ1
//COMMON /RAIN1/ RAIN_FACTOR
//COMMON /MDDATA1/ MAXRNO
//common /ts_save/ ts_old, tave_old, Ts_Ave, Tave_ave

using namespace std;

ifstream topinpFile, modelspcFile, rchareasFile, rainFile;
#if TRACE
    ofstream traceFile;
    string caller;
#endif
#ifdef ZBAR_OUT
    ofstream zbarFile;
#endif
//ofstream debugFile;

int main()
{
	int initT[iex], iend[iex], Neq=0, i, Npar=0;
    int mfit[iex];
	int it, iflag, iopt, npx, nfor;
	int ibeale, ifit=0;
	ArrayXXd qact(Nrx, iex);
	double actime[Nrx];
	double qfit[iex], par[dpm];

	char prt;
	string modelid;
#if TRACE
    traceFile.open("results/trace.dat");
    caller = "main()";
    traceFile << caller << endl;
#endif
#ifdef ZBAR_OUT
    zbarFile.open("results/zbar.dat");
#endif
    //debugFile.open("results/debug.dat");
    inputT(initT, iend, Neq, qact, actime, modelid, Npar, Nrx, iex);
	if (Neq > iex)
		cout << "Too many responses - increase iex dimension\n";
    if (Npar > dpm)
		cout << "Too many parameters - increase dpm dimension\n";
	it    = 1;
	iflag = 1;
	iopt  = 0;
	prt   = 'Y';
	npx   = dpm;
	nfor  = 6;

	//------------ topinp.dat -------------------
	std::string inLine;
	std::istringstream iss;
	ifstream topinpFile;
	topinpFile.open("topinp.dat");
	if (topinpFile.good()) {
        cout << "topinp.dat opened for reading in main()\n";
    } else {
        cout << "topinp.dat not found\n";
        exit(1);
    }
    for (i = 0; i < 7; ++i) {
        getline(topinpFile, inLine, '\n');
    }
	for (i = 0; i < Npar; i++) {
	    getline(topinpFile, inLine, '\n');
        iss.str(inLine);
        iss >> par[i];
        iss.clear();
	}
	topinpFile.close();
    //------------ topinp.dat -------------------

	for (i = 0; i < max(1,Neq); i++) { // When Neq = 0, i <= 1.
		mfit[i] = 1;	// mfit is an array of flags
	}
	ibeale = 0;
    cout << "Starting model\n";  // For an unkonwn reason the release version crashes without this line
	Model(it, iflag, iopt, prt, Neq, Npar, npx, iex, nfor, qfit, par, mfit, ifit, ibeale);
#if TRACE
    traceFile.close();
#endif
#ifdef ZBAR_OUT
    zbarFile.close();
#endif
    //debugFile.close();
    return 0;
}
