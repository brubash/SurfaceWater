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
#include <dirent.h>

// This version, V8, has had extra checking of rainfall sites added plus,
// input for lakes.
// V9 also has checking/matching of the flow sites with reaches

using namespace std;
using namespace Eigen;

int mdData(int &Ngauge, int &Ns, int &Nrch, int *Nka, double *tl, double **atb, double **pka, int *Nd, double **cl2, double **pd2,
    double &units, int *ll, int *Ntr, int *Nts, ArrayXXi &linkS, ArrayXXi &linkR, double **si,
    ArrayXXd &Sp, double **Rp,  const int iret, ArrayXXi &pMap, int &Npar, ArrayXXi &lrg,
    ArrayXXd &wrg, const int iex, int *llout, int &Neq, int &Nout, int &nBout, int *iBout, int &nRchSav, int *Qmap,
    double *rel, int &relFlag, double *minSp, double *maxSp, double *minSi, double *maxSi, double *minRp, double *maxRp, bool &limitC,
    ArrayXXi &ClinkR, ArrayXi &kllout, int *ishift0, const int maxInt, const int maxGuage, const int maxSlp,
    const int maxResponse, const int maxA, const int maxC, const int maxChn, const int maxRchAreas, const int maxSites,
    ArrayXXd &bp, double *bXlat, double *bXlon, ArrayXXd &wrg1)
{
    // *********************************************************************
    //  This subroutine represents an amalgamation of subroutines READPS and
    //  SET_LINKR . It reads the disk file containing the hydrological data
    //  and returns the data in several arrays and variables.

    //      Variable list :-
    //  a)  Variables received ( via the include file 'tdims_v7.INC' )
    //      ------------------
    //  MAXSLP -- allowable number of subbasins
    //  MAXCHN --   "      "    "  reaches.

    //  b)  Variables and arrays returned
    //      -----------------------------
    //  NGAUGE -- number of raingauges used.
    //  NS -- actual number of subbasins.
    //  NRCH --  "     "    "  reaches.
    //  LL(MAXCHN) -- pointer giving the order for processing all reaches.
    //  NTR(MAXCHN) -- no. of times this reach is used as input to another
    //                 reach.
    //  NTS(MAXSLP) -- no. of times this basin is used
    //  LINKS(2,MAXSLP) -- array containing basin-raingauge link data ,
    //    i.e. LINKS[0][i-1] = arbitrary numbering (must be between 1 and NS).
    //         LINKS[1][i-1] = no. of the raingauge to be used for this basin.
    //                       (must be between 1 and NGAUGE)

    //  LINKR(4,MAXCHN) -- array containing reach network information ,
    //    i.e. linkR(0,i-1) = arbitrary numbering (must be greater than NS).
    //         linkR(1,i-1) = basin or upstream reach feeding into this reach.
    //         LINKR(3,I) = second basin or upstream reach feeding into this
    //                      reach.
    //         LINKR(4,I) = number of basins or upstream reaches feeding into
    //                      this reach.

    //  SI(NSI,MAXSLP) -- basin initial conditions ,
    //  minsi(NSI) -- lower limit of basin initial condition
    //  maxsi(NSI) -- upper limit of basin initial condition

    //  SP(NSP,MAXSLP) -- basin properties , Only ISP actually used.
    //  minsp(NSP) -- lower limit of basin properties
    //  maxsp(NSP) -- upper limit of basin properties

    //  RP(NSP,MAXCHN) -- reach properties ,
    //    i.e. RP[0][i-1] = basin.
    //         RP[1][i-1] = Manning's n.
    //         RP(3,I) = width .
    //         RP(4,I) = length
    //  minrp(NRP) -- lower limit of reach properties
    //  maxrp(NRP) -- upper limit of reach properties

    //  PAR(4,dpm) -- Parameter map.  Only Npar actually used.
    //   This controls the mapping of subbasin properties that may be
    //   controlled by an optimizer.
    //
    //  LRG(MAXSLP,MAXGAUGE), WRG(MAXSLP,MAXGAUGE)
    //             -- Extended basin-raingauge link data
    //    i.e. LRG(0,i-1) = raingauge numbers to be used for sub-basin I.
    //         WRG(0,i-1) = raingauge weight to be used for
    //                      raingauge LRG[0][i-1] on sub-basin I
    // limitc = true if a file of constraint values is to be used, false otherwise

    //   llresponse(maxchn) - vector of reaches whose flow will be used as response time series
    //   llout(maxchn) - vector of reaches whose flow is to be output

    //  IRET -- This is for if user presses ESCAPE
    // *********************************************************************

    //  Passed variables
    //  ****************
    int n0rch;
    int iMatch;
    double temp;
    double Rtemp1, Rtemp2, Ptemp1, Ptemp2, dummy;
    int rchNo[maxRchAreas];
    int Qsite1[maxRchAreas], Qsite2[maxResponse];
    int Qeast[maxRchAreas], Qnorth[maxRchAreas], Qsite[maxRchAreas];
    double relTemp[maxRchAreas];
    int iShift0temp[maxRchAreas];
    double flat,clat,flong,clong;

    //  Local variables
    //  ***************
    int i, j, ii, jj, Line, maxG, n1rch, iCheck;
    Array<int,Dynamic,1> lln(maxChn);
    int Ntsg[maxSlp][maxChn], iTbout[maxSlp];

    string verno, ver_msp;
    string fName;
    bool rchAreas;
    bool allzer, us1, us2, high;
    int sites[maxSites];

    int isp, isum, Ntri, k;
    int Nrchno, krch, neq_temp, neq_chk;
    int NgBase;
    char cr, comma; // carriage return, comma separator
    bool exist;
    struct dirent *dirp;
    DIR *dp;
    string testStr, inLine;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> mdData(" << ncalls << ")" << std::endl;
    }
    caller = "mdData";
#endif
    // COMMON /RAIN1/ RAIN_FACTOR
    // added this because we now need the value of MAXRNO in CALCTS
    // COMMON /MDDATA1/ MAXRNO

    for (i = 0; i < maxRchAreas; i++) {
        rchNo[i] = 0;
    }
    // ************************************************
    ifstream modelspcFile("modelspc.dat");	// unit lunrd
    if (!modelspcFile.is_open()) {
        cerr << "Failed to open modelspc.dat\n";
        exit(EXIT_FAILURE);
    }
    ifstream interpweightFile("interpweight.dat");	// Interpolation weights
    if (!interpweightFile.is_open()) {
        cerr << "Failed to open interpweight.dat\n";
        exit(EXIT_FAILURE);
    }
    getline(interpweightFile, inLine, '\n');           // Comment line

    // We need a separate grid of interpolation weights for temperature etc because precip weights
    // do not work for temperature

    modelspcFile >> ver_msp;
    getline(modelspcFile, inLine, '\n');           // Read the remainder of the heading line
    getline(modelspcFile, inLine, '\n');           // and the next comment line.
    Line = 1;

    // Read no. of raingauges , basins and reaches used
    //  RAIN_FACTOR is the multiplier to be applied to meso-scale rainfall
    rain::rain_factor = 1.0;   // Set default
    if (ver_msp == "Ver2" || ver_msp == "ver2") {
        modelspcFile >> Ngauge >> Ns >> Nrch;
        modelspcFile >> rain::rain_factor;
        getline(modelspcFile, inLine, '\n');    // Read the rest of the line, or at least the line end.
        // check RAIN_FACTOR for reasonableness
        if (rain::rain_factor < 0.5 || rain::rain_factor > 4.0) {
            cerr << " ***** Your rainfall multiplier, 4th number on ";
            cerr << "the 3rd line of MODELSPC.DAT, is unreasonable. ";
            cerr << "Its value is" << fixed << setw(10) << setprecision(2) << rain::rain_factor << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        modelspcFile >> Ngauge >> Ns >> Nrch;
        getline(modelspcFile, inLine, '\n');    // Read the rest of the line, or at least the line end.
    }

    if (Ngauge > maxGauge) {
        cerr << " Number of raingauges exceeds maximum of " << dec << setw(3) << maxGauge << '\n';
        cerr << " Error in data in the model description file \n";
        exit(EXIT_FAILURE);
    }
    if (Nrch > maxChn) {
        cerr << " Number of reaches exceeds maximum of " << dec << setw(3) << maxChn << '\n';
        cerr << " Error in data in the model description file \n";
        exit(EXIT_FAILURE);
    }
    if (Ns > maxSlp) {
        cerr << " Number of basins exceeds maximum of " << dec << setw(3) << maxSlp << '\n';
        cerr << " Error in data in the model description file \n";
        exit(EXIT_FAILURE);
    }
    // Read number of basin properties and parameters
    getline(modelspcFile, inLine, '\n');    //  Comment line
    modelspcFile >> isp >> Npar;	//  isp is number of slope properties
    getline(modelspcFile, inLine, '\n');    // Read the rest of the line, or at least the line end.
    if (isp > Nsp) {
        cerr << " Number of basin properties exceeds maximum of " << dec << setw(3) << Nsp << '\n';
        cerr << " Error in data in the model description file \n";
        exit(EXIT_FAILURE);
    }
    if (Npar > dpm) {
        cerr << " Number of parameters exceeds maximum of " << dec << setw(3) << dpm << '\n';
        cerr << " Error in data in the model description file \n";
        exit(EXIT_FAILURE);
    }
    // Moved two lines below from later to be logically before
    // all element properties repeat
    getline(modelspcFile, inLine, '\n');	//  Comment line
    modelspcFile >> units;
    getline(modelspcFile, inLine, '\n');    // Read the rest of the line, or at least the line end.
    // Read parameter map here before all element properties repeat
    // modified below to include variable length parameter map
    for (i = 1; i <= 6; i++) {   			//  6 comment lines
        getline(modelspcFile, inLine, '\n');
    }
    for (i = 1; i <= Npar; i++) {
        for (j = 1; j <= 4; j++) {
            modelspcFile >> pMap(j-1,i-1);
        }
        getline(modelspcFile, inLine, '\n');	// Read the remainder of the line.
        //      checks to ensure valid parameter map values
        if (pMap(0,i-1) == 1) {
            //	    if(pmap[1][i-1] .le. 0 .or. pmap[1][i-1] .gt. isp)go to 296
        } else if (pMap(0,i-1) == 2) {
            if (pMap(1,i-1) <= 0 || pMap(1,i-1) > Nrp) {
                cerr << " Error in the parameter map data in the model description file\n";
                exit(EXIT_FAILURE);
            }
        } else if (pMap(0,i-1) == 3) {
            if(pMap(1,i-1) <= 0 || pMap(1,i-1) > Nsi) {
                cerr << " Error in the parameter map data in the model description file\n";
                exit(EXIT_FAILURE);
            }
        } else {
            cerr << " Error in the parameter map data in the model description file\n";
            exit(EXIT_FAILURE);
        }
    }

    // check that sufficient and correct values given for pmap
    isum = 0;
    for (i = 0; i < Npar; i++) {
        isum += pMap(3,i);
    }

    if (isum != (Npar*(Npar+1)/2)) {
        cerr << " Error in the parameter map data in the model description file\n";
        exit(EXIT_FAILURE);
    }
    //  modified above for general number of parameters

    //  Read and check basin-raingauge data
    //  ***********************************
    getline(modelspcFile, inLine, '\n');	//  Comment line

    // USE LLN TEMPORARILY TO CHECK THAT ALL BASINS HAVE GOT A RAINGAUGE
    for (i = 0; i < Ns; i++) {
        lln(i) = 0;
    }
    maxG = 0;
    // get raingauge numbers from rain.dat and set up mapping to column #s
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
        rainFile.open(rainFileName);	// Fortran unit lunOpt
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
    if (verno == "Ver1" || verno == "ver1") {
        if (ver_msp != "Ver1" && ver_msp != "ver1" &&  ver_msp != "Ver2" && ver_msp != "ver2") {
            cerr << " Incompatible modelspc.dat and " << rainFileName << '\n';
            exit(EXIT_FAILURE);
        }
    }

    if((verno == "Ver1" || verno == "ver1" || verno == "Ver2" || verno == "ver2") &&
            (ver_msp == "Ver1" || ver_msp == "ver1" || ver_msp == "Ver2" || ver_msp == "ver2")) {

        getline(rainFile, inLine, '\n');
        getline(rainFile, inLine, '\n');
        rainFile >> verno >> Ntri;
        for (i = 1; i <= Ntri; i++) {
            rainFile >> sites[i-1];
        }
    } else {
        getline(rainFile, inLine, '\n');
        getline(rainFile, inLine, '\n');
        rainFile >> verno;
    }
    rainFile.close();
//===================================================================================

    for (i = 1; i <= Ns; i++) {
        Line = 2;
        // Allow >1 raingauge per sub-basin 4-May-98
        // if LINKS[1][i-1]<0, it means there is a list of -LINKS[1][i-1] pairs
        // of (raingauge #, weight) which are to be used for that sub-basin
        // So if the line in the file contains '1 2', then we have usual case:
        // LINKS(1:2,I) is [1 2], and we use gauge 2 for basin 1 as usual.

        // But if the line in the file contains '1 -2   3 0.3   4 0.5', then
        // the '-2' indicates two pairs of (gauge#,weight), which are taken
        // to mean 'subbasin 1 rainfall = gauge#3*0.3 + gauge#4*0.5
        for (k = 1; k <= maxGauge; k++) {
            lrg(i-1,k-1) = 0;
            wrg(i-1,k-1) = 0.0;
            wrg1(i-1,k-1) = 0.0;  // Interpolation weights
        }
        interpweightFile >> linkS(0,i-1) >> linkS(1,i-1);
        for (k = 1; k <= -linkS(1,i-1); k++) {
            interpweightFile >> lrg(i-1,k-1) >> wrg1(i-1,k-1);
        }
        interpweightFile.get(cr);	// read line end
        if (int(cr) == 13) {  // If MSDOS file end ( carriage return )
            interpweightFile.get(cr);	// read line feed
        }

        // Read the same data in again to the same arrays from a different file
        modelspcFile >> linkS(0,i-1) >> linkS(1,i-1);
        for (k = 1; k <= -linkS(1,i-1); k++) {
            modelspcFile >> lrg(i-1,k-1) >> wrg(i-1,k-1);
        }
        getline(modelspcFile, inLine, '\n');    // Read the rest of the line, or at least the line end.

        //if(verno.eq.'Ver1'.or.verno.eq.'ver1') then    ! we have new format rain.dat
        if(verno == "Ver1" || verno == "ver1" || verno == "Ver2" || verno == "ver2") {
            // In a new modelspc.dat the sub-basin raingauge relationship is based
            // on raingauge numbers, not the column numbers of data. We need to find
            // column positions in rain.dat corresponding to each raingauge number
            // to be compatible with early code.
            for (k = 1; k <= -linkS(1,i-1); k++) {
                // Set IMATCH to 0 to show that not yet found a matching site
                iMatch = 0;
                for (j = 1; j <= Ntri; j++) {
                    // Introduced the following code to allow for meso-scale
                    // rainfall grid variations

                    // Introduced the if(imatch .eq. 0) below to make sure that the chk_sites is only used
                    // to find a match once, because upon finding a match the value lrg(j,k) is changed and an
                    // error can occur if after it being changed the changed value matches to a different site.
                    if (iMatch == 0) {
                        chk_sites(lrg(i-1,k-1), sites[j-1], j, iMatch);
                    }
                }
                if (iMatch == 0) {
                    cerr << Ntri << " " << -linkS(1,i-1) << " " << lrg(i-1,k-1) << " " << sites[j-1] << " " << j << endl;
                    cerr << " Error in data in the model description file:\n";
                    cerr << " The number of raingauges used for the basins is inconsistent with (more than) \n";
                    cerr << " that specified in the first data line\n";
                    cerr << " This may arise because of a mismatch in the grid, check ";
                    cerr << " around site number " << dec << setw(9) << maxG << '\n';
                    exit(EXIT_FAILURE);
                }
            }
            // Translate lrgs which are site #s to column position in rain.dat
            if (linkS(1,i-1) > 0) {
                lrg(i-1,0)  = linkS(1,i-1);
                wrg(i-1,0)  = 1.0;
                wrg1(i-1,0) = 1.0;
            }
        } else {
            // RPI 5/6/2003 added this else clause for compatibility with old rain.dat files
            // Old rain.dat files assume column 1 is for the first mentioned site etc.
            // Therefore need to set up  a correct set of lrg values, & discard site #s
            for (k = 1; k <= -linkS(1,i-1); k++) {
                lrg(i-1,k-1) = k;
            }
        }
        // translate lrgs which are site #s to column position in rain.dat RPI & RAW 25/8/00
        if (linkS(1,i-1) > 0) {
            lrg(i-1,0)  = linkS(1,i-1);
            wrg(i-1,0)  = 1.0;
            wrg1(i-1,0) = 1.0;
        }

        // USE LLN TEMPORARILY TO CHECK THAT ALL BASINS HAVE A RAINGAUGE
        lln(linkS(0,i-1)-1)++;
        // Altered following test to allow negative raingauge#, if done properly
        if (linkS(1,i-1) <= 0 && lrg(i-1,0) == 0) {
            cerr << " Error reading the model description file in mdData() ~line 420\n";
            exit(EXIT_FAILURE);
        }
    }

    // Altered to allow use of fewer raingauges than are in the data file

    // check that all basins have a raingauge
    // we should have seen each  basin exactly once
    for (i = 0; i < Ns; i++) {
        if ( lln(i) != 1 ) {
            cerr << " no raingauge specified for basin, lln(" << dec << setw(3) << i << ") " << lln(i) << '\n';
            exit(EXIT_FAILURE);
        }
    }

    //  read reach link data
    //  ********************

    getline(modelspcFile, inLine, '\n');               // comment line
    mddata::maxRno = 0;
    for (i = 1; i <= Nrch; i++) {
        Line = 3;
        for (j = 1; j <= 3; j++) {
            modelspcFile >> linkR(j-1,i-1);
        }
        getline(modelspcFile, inLine, '\n');	// read the rest of the line
        for (j = 1; j <= 3; j++) {
            ClinkR(j-1,i-1) = linkR(j-1,i-1);
        }
        mddata::maxRno = max(mddata::maxRno, linkR(0,i-1));  // keep track of largest reach number
        if (linkR(0,i-1) <= 0 || linkR(1,i-1) <= 0) {
            cerr << " Error reading the model description file in mdData() ~line 468\n";
            exit(EXIT_FAILURE);
        }
    }

    //  sequence to read response and output reach numbers 5.26.98
    getline(modelspcFile, inLine, '\n');	// comment line
    modelspcFile >> Neq >> Nout;            // nFlowRecorders, nFlowOutLocations
    getline(modelspcFile, inLine, '\n');	// Read the rest of the line.
    if (Nout > maxChn) {
        cerr << "ERROR - cannot output more than " << maxChn << " channels\n";
    }

    if (max(Neq, Nout) > maxResponse) {
        cerr << " Error in data in the model description file.\n";
        cerr << " neq or nout more than maxresponse " << Neq << " " << Nout << " " << maxResponse << '\n';
        exit(EXIT_FAILURE);
    }
    if (Neq > iex) {
        cerr << " Error in data in the model description file.\n";
        cerr << "neq more than dimension allows" << Neq << '\n';
        exit(EXIT_FAILURE);
    }
    //  the first neq of the following list will be used as response time series
    //  Any remaining are output.
    kllout = 0;	// array assignment (response to valgrind complaint)
    for (i = 1; i <= max(Neq, Nout); i++) {
        modelspcFile >> kllout(i-1);
    }
    getline(modelspcFile, inLine, '\n');	// Read the rest of the line.

    // If nudging a forecast the some of the rel(i)  will != zero. Set
    // relflag to one in this case, zero otherwise. Initial relflag here
    // in case there is no rchareas.txt
    relFlag = 0;

    // This is where we check the reach numbers against those
    // in runoff.dat using rchareas.txt. Initially rchareas.txt contained
    // reach no., area upstream of reach, and site no. These have been supplemented
    // with easting and northing. This addition should not affect the Matlab
    // Topplot procedure since it reads only the first 3 values.
    rchAreas = false;
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
        cerr << " ****Warning, there is no rchareas.txt file\n";
        cerr << " Data in streamflow_calibration.dat assumed in correct order\n";
    } else {
        if (!rchareasFile.is_open()) {
            rchareasFile.open("rchareas.txt");	// Fortran unit 82
            if (!rchareasFile.is_open()) {
                cerr << "Failed to open rchareas.txt\n";
                exit(EXIT_FAILURE);
            }
        } else {	// rewind
            rchareasFile.clear();
            rchareasFile.seekg(0);
        }

        rchAreas = true;
        //	Read in rcharea.txt - may need to revise this if more or less than
        // maxresponse values in rchareas.txt
        for (i = 1; i <= maxRchAreas; i++) {
            rchareasFile >> rchNo[i-1] >> dummy >> Qsite[i-1] >> Qeast[i-1] >> Qnorth[i-1] >> rel[i-1] >> ishift0[i-1];
        }
        Nrchno = i;
        rchareasFile.close();
    }

    // get flow site numbers from runoff.dat and set up mapping to reaches RPI 18/9/01
    exist = false;
    if ((dp = opendir(".")) == NULL) {
        cerr << "Can't open directory\n";
        exit(EXIT_FAILURE);
    }
    while ((dirp = readdir(dp)) != NULL) {
        testStr = dirp->d_name;
        if (testStr == "streamflow_calibration.dat") {
            exist = true;
            closedir(dp);
            break;
        }
    }
    if (exist == false) {
        cerr << "  No streamflow_calibration.dat available in MDDATA\n";
        exit(EXIT_FAILURE);
    }
    ifstream streamflow_calibrationFile("streamflow_calibration.dat");	// Fortran unit 82
    if (!streamflow_calibrationFile.is_open()) {
        cerr << "Failed to open streamflow_calibration.dat\n";
        exit(EXIT_FAILURE);
    }

    getline(streamflow_calibrationFile, inLine, '\n');
    getline(streamflow_calibrationFile, inLine, '\n');
    getline(streamflow_calibrationFile, inLine, '\n');
    streamflow_calibrationFile >> verno;
    // rewind
    streamflow_calibrationFile.clear();
    streamflow_calibrationFile.seekg(0);

    // Expanded the commented out if stmt below for ver2 modelspc.dats
    if (verno == "Ver1" || verno == "ver1") {
        if (ver_msp != "Ver1" && ver_msp != "ver1" && ver_msp != "Ver2" && ver_msp != "ver2") {
            cerr << " Incompatible MODELSPC.DAT and streamflow_calibration.dat\n";
            exit(EXIT_FAILURE);
        }
    } else {
        //           if(ver_msp == "Ver1" || ver_msp == "ver1" ||
        //     +		ver_msp == "Ver2" || ver_msp == "ver2") then
        if (ver_msp == "Ver1" || ver_msp == "ver1") {
            cerr << "Warning - Incompatible MODELSPC.DAT and streamflow_calibration.dat\n";
            cerr << "    It is the users responsibility to get order of\n";
            cerr << "    sites and data columns in streamflow_calibration.dat in correct order\n";

        }
    }
    //	if( ((verno == "Ver1" || verno == "ver1") &&
    //     +      (ver_msp != "Ver1" && ver_msp != "ver1")) .or.
    //     +	((verno != "Ver1" && verno != "ver1") &&
    //     +      (ver_msp == "Ver1" || ver_msp == "ver1"))) then
    //			write(21,*) " Incompatible MODELSPC.DAT and RUNOFF.DAT"
    //			stop
    //	endif
    // Initialise qmap so that OK if an old file
    for (i = 1; i <= Neq; i++) {
        Qmap[i-1] = i;
    }

    // Expanded the commented out if stmt below for ver2 modelspc.dats
    // Changed the "if" statement to allow for old style "runoff.dats"
    if (verno == "Ver1" || verno == "ver1") {
        if (ver_msp == "Ver1" || ver_msp == "ver1" || ver_msp == "Ver2" || ver_msp == "ver2") {
            //      if((verno == "Ver1" || verno == "ver1")  &&
            //     +   (ver_msp == "Ver1" || ver_msp == "ver1" ||
            //     +    ver_msp == "Ver2" || ver_msp == "ver2")) then
            // We have a new format Runoff.dat and a new format Modelspc.dat - we must also
            // have a Rchareas.txt - if not stop
            if (!rchAreas) {
                cerr << " No RCHAREAS.TXT file - execution stops\n";
                exit(EXIT_FAILURE);
            }	// l8
            //	if((verno == "Ver1" || verno == "ver1")  &&
            //     +   (ver_msp == "Ver1" || ver_msp == "ver1")) then
            streamflow_calibrationFile >> verno >> Ntri;
            for (i = 1; i <= Ntri; i++) {
                streamflow_calibrationFile >> Qsite2[i-1];
            }
            // There must be neq columns of data in runoff.dat

            if(Ntri != Neq) {
                cerr << " ****ERROR The number of columns in streamflow_calibration.dat\n";
                cerr << " is different from the number," << Neq << " of reaches\n";
                cerr << " in modelspc.dat\n";
                exit(EXIT_FAILURE);
            } // l3
            // Compare reaches in modelspc.dat, in kllout, with those just read in
            // and select the subset of sites from rchareas.txt. Can't do this until
            // we have qsite2 because some sites may fall in the same reach RPI 22/5/2002
            krch = 0;
            // put in the neq=0 test etc.
            neq_temp = Neq;
            //		if(neq == 0) neq_temp=ntri
            for (i = 1; i <= neq_temp; i++) {
                for (j = 1; j <= Nrchno; j++) {
                    if ( kllout(i-1) ==  rchNo[j-1]) {
                        // Make sure you have the correct site as well
                        for (k = 1; k <= neq_temp; k++) {
                            if (Qsite[j-1] == Qsite2[k-1])
                                goto L3300; // We have a match
                        }
                        goto L3302;
L3300:
                        krch++;
                        Qsite1[krch-1] = Qsite[j-1];
                        // We need the rel(i) in the same order as the qsite1 otherwise
                        // in CALCTS the wrong rel values are applied to weighting modelled and corrected
                        // modelled flows in the "nudging" of forecasts.
                        relTemp[krch-1] = rel[j-1];			// j>=krch so this should not be a problem
                        iShift0temp[krch-1] = ishift0[j-1];	// j>=krch so this should not be a problem
                        if (relTemp[krch-1] > 0)
                            relFlag = 1;
                        goto L3301;
                    }
L3302:
                    ;
                }
L3301:
                ;
            }
            // Now put reltemps into rel
            for (k = 1; k <= krch; k++) {
                rel[k-1] = relTemp[k-1];
                ishift0[k-1] = iShift0temp[k-1];
            }
            // We now have the site numbers of the data in Runoff.dat (the qsite2s)
            // and the mapping of these to reach numbers from rchareas.txt. Now we
            // need to work out which column in runoff.dat corresponds to each reach
            // number in modelspc.dat. Initially we will ignore the easting and northing
            // data and focus on the site number to reach number relationship.

            // We can't proceed any further here because the flow data have not yet
            // been read in, and since we are going to
            // take the cautious approach, i.e. rearrange the columns in
            // array flow to match the order of the sites in kllout, this will have
            // to be done in inputt after the flows have been read in by hyData.
            // This will effectively
            // set flow up the way it would have been read prior to this development.
            // The alternative high risk approach is to leave the columns in flow the
            // way they were read in, and alter every occurrence of the use of flow to
            // reflect the new mapping.

            // Set up the mapping as follows: use an array called qmap. Define
            // qmap(i) to be the column in runoff.dat corresponding to i'th reach
            // in the list read in from modelspc.dat, e.g. through the reach-site no.
            // linkage data for the 2nd reach in the list might be in the 4th column,
            // in which case qmap(2)=4
            // read in from modelspc.dat.
            // We can only do this if there is a rchareas.txt
            if (rchAreas) {
                neq_chk = 0;
                // Introduced the max and changed k to krch because there
                // are legimate cases where more than one site number can be attached to
                // the same reach and in this case qsite1 has more than neq entries in it.
                for (i = 1; i <= max(Neq, krch); i++) {
                    for (j = 1; j <= neq_temp; j++) {
                        if(Qsite1[i-1] == Qsite2[j-1]) {
                            // 	qmap(i)=j
                            neq_chk++;
                            Qmap[neq_chk-1] = j;
                        } // l6
                    }
                }
                if(neq_chk != neq_temp) {
                    cerr << " ****ERROR - Not all the reaches in MODELSPC.DAT";
                    cerr << " could be found in RCHAREAS.TXT\n";
                    exit(EXIT_FAILURE);
                } // l7

                // If both verno and ver_msp are != Ver1 then assume old files
                // and it is the users responsibility to get things right!!!
                // RPI 6/6/2003 - next 2 lines are irrelevant here!
                //		else
                //			read(82,*) verno
            } // l5
        } // l2
        // RPI 6/6/2003 added the "else" clause for backwards compatibility with old
        // style format files. Come here if old style "modelspc.dat" file
    } else {
        // Decide if there is a RCHAREAS.TXT. If there isn't then there can be no nudging
        // of forecasts and we have to assume that the order of data in RUNOFF.DAT is the
        // same as the order of the reaches in MODELSPC.DAT. If there is a RCHAREAS.TXT then
        // we can do nudging but still have to assume the reaches in RUNOFF.DAT are the
        // same order as in MODELSPC.DAT
        if (rchAreas) {
            for (i = 1; i <= Neq; i++) {
                for (j = 1; j <= Nrchno; j++) {
                    if ( kllout(i-1) == rchNo[j-1] ) {
                        Qsite1[i-1] = Qsite[j-1];
                        relTemp[i-1] = rel[j-1];
                        iShift0temp[i-1] = ishift0[j-1];
                        if (rel[j] > 0.0)
                            relFlag = 1;
                        goto L4301;
                    }
                }
L4301:
                ;
            }
            // Now put reltemps into rel
            //for (k = 1; k <= Neq; k++) {

            //		rel(k)=reltemp(k)
            //		ishift0(k)=ishift0temp(k)
            //}
        } else {
            // No RCHAREAS.TXT so no nudging possible
            relFlag = 0;
        }
        //  end of additions
    } // l1
    //	endif
    streamflow_calibrationFile.close();

    // End of changes needed to match flows in runoff.dat to reaches in
    // modelspc.dat using site numbers
    // -------------------------------------------------------------------------

    // At this point we have all the data needed to sort out the how to update
    // the ZBAR0s for each measured sub-basin
    //     call SET_COR_DATA(linkR,neq,nout,nrch,ns,kllout)
    //                       ----- NOTE ------
    //  From this point onwards until label 220 the reach link data will be
    //  checked and processed , and only after that more data will be read
    //  from the disk file.

    // TOPMODEL: Create artificial gutters between subcatchments and reaches
    // The first NRCHSAV lines are all reaches fed by subcatchments or
    // reaches. Whenever we find a subcatchment, convert it to a gutter, and
    // and add in a new line which says our new gutter is fed by that
    // subcatchment. The values of N0RCH and LINKR(4,.) will be done later on

    nRchSav = Nrch;
    //  DGT - Changes below to generalize the numbering of gutters
    //   in case there is a reach number bigger than 100.

    //      ngbase = max(100, maxrno)
    if (nooksack == 1) {
        for (i = 1; i <= Nrch; i++) {
            ll[i-1] = i;
        } // not the usual routing
        n0rch = 0;
    } else {
        NgBase = max(100, mddata::maxRno);
        for (i = 1; i <= nRchSav; i++) {
            for (j = 1; j <= 2; j++) {
                ii = linkR(j,i-1);
                if (ii <= Ns && ii > 0) {
                    // Have found a subcatchment, so convert it to a gutter
                    // LINKR(J+1,I) = 100+II
                    linkR(j,i-1) = NgBase + ii;
                    // Add a new line for the gutter
                    Nrch++;
                    // The gutter (100+ii) is fed by subcatchment ii
                    // LINKR(1,NRCH) = 100+II
                    linkR(0,Nrch-1) = NgBase + ii;
                    linkR(1,Nrch-1) = ii;
                    linkR(2,Nrch-1) = 0;
                }
            }
        }
        //  Find 0th order reaches and fill in 4th column of LINKR
        //  ******************************************************
        n0rch = 0;
        for (i = 1; i <= Nrch; i++) {
            linkR(3,i-1) = 0;
            if (linkR(1,i-1) > 0 && linkR(2,i-1) == 0)
                linkR(3,i-1) = 1;
            if (linkR(1,i-1) > 0 && linkR(2,i-1) > 0)
                linkR(3,i-1) = 2;
            if (linkR(3,i-1) == 0) {
                //       ELSE
                //         IF(LINKR(3,I).EQ.0) LINKR(4,I) = 1
                //         IF(LINKR(3,I).GT.NS) LINKR(4,I) = 2
                //         IF(LINKR(3,I).LE.NS.AND.LINKR(3,I).NE.0) THEN
                cerr << " Error reading the model description file, i = " << dec << setw(3) << i << '\n';
                cerr << " error in the reach link data , reach number " << dec << setw(3) << linkR(0,i-1) << '\n';
                exit(EXIT_FAILURE);
            }
            if (linkR(1,i-1) <= Ns && linkR(2,i-1) <= Ns) {
                n0rch++;
                lln(n0rch-1) = linkR(0,i-1);
                ll[n0rch-1]  = i;
            }
        }
        //  Fill in 4th column of LINKR for 0th order reaches
        //  *************************************************
        //      DO 65 I=1,N0RCH
        //       II = LL(I)
        //       IF(LINKR(2,II).LE.0) THEN
        //         WRITE(21,55)LINKR(1,II)
        //         GOTO 290
        //       ELSE
        //         IF(LINKR(3,II).GT.0.AND.LINKR(3,II).LE.NS) LINKR(4,II) = 2
        //         IF(LINKR(3,II).LE.0) LINKR(4,II) = 1
        //         IF(LINKR(3,II).GT.NS) THEN
        //         WRITE(21,55)LINKR(1,II)
        //           GOTO 290
        //         ENDIF
        //       ENDIF
        //  65 CONTINUE
        //  Check to see if all the reaches are not zero order reaches.
        //  Check also if N0RCH > MAXGUT or N0RCH = 0 , both fatal errors.
        //  **************************************************************
        allzer = false;
        //     IF(N0RCH.EQ.NRCH) THEN
        //       WRITE(21,67)
        //  67   FORMAT(" WARNING! All the reaches in the river network are ",
        //    1         "zero order reaches")
        //       ALLZER = .TRUE.
        //       GOTO 180
        //     ENDIF
        //  Test each "non-zero" order segment or link to see if it gets all its
        //  inputs from segments above* it in the stream network. (* 'Above' is
        //  only figurative and it really refers to those reaches which already
        //  have their relative positions allocated in the array LL. This is a
        //  recursive process and it carries on until all the segments have been
        //  allocated a position in the array LL. )
        //  ********************************************************************
        n1rch = n0rch;
        iCheck = n1rch;
L70:
        ;
        for (i = 1; i <= Nrch; i++) {
            //  Goto the end of the loop if reach I has already been assigned a position.
            for (j = 1; j <= n1rch; j++) {
                if (ll[j-1] == i)
                    goto L100;
            }
            if(linkR(3,i-1) == 1) {
                //  Enter here if only one upstream link.
                for (j = 1; j <= n1rch; j++) {
                    if (linkR(1,i-1) == lln(j-1)) {
                        n1rch++;
                        lln(n1rch-1) = linkR(0,i-1);
                        ll[n1rch-1] = i;
                        goto L100;
                    }
                }
            } else if (linkR(3,i-1) == 2) {
                //  Enter here if two upstream links.
                us1 = false;
                us2 = false;
                for (i = 1; i <= n1rch; i++) {
                    if (linkR(1,i-1) == lln(j-1) || linkR(1,i-1) <= Ns)
                        us1 = true;
                    if (linkR(2,i-1) == lln(j-1) || linkR(2,i-1) <= Ns)
                        us2 = true;
                    if (us1 && us2) {
                        n1rch++;
                        lln(n1rch-1) = linkR(0,i-1);
                        ll[n1rch-1] = i;
                        goto L100;
                    }
                }
            }
            //  Enter here to check if entire array LL has been filled out.
L100:
            if (n1rch == Nrch)
                goto L120;
        }
        //  Enter here after looping and check if any new position has been
        //  out in array LL . If not then the link data is corrupt.
        if (n1rch == iCheck) {
            cerr << " Link data corrupt after link " << dec << setw(3) << lln(n1rch-1) << '\n';
            cerr << " This is the order in which the river links will be processed so far :\n";
            for (i = 1; i <= n1rch; i++) {
                cerr << dec << setw(5) << lln(i-1);
            }
            exit(EXIT_FAILURE);
        }
        //  Start looping from beginning again.
        iCheck = n1rch;
        goto L70;
L120:
        ;
        //  Fill in the vector NTR
        //  **********************
        for (i = 1; i <= Nrch; i++) {
            Ntr[i-1] = 0;
        }
        for (i = n0rch+1; i <= Nrch; i++) {
            ii = ll[i-1];
            for (j = 1; j <= Nrch; j++) {
                jj = linkR(0,j-1);
                if (linkR(1,ii-1) == jj)
                    Ntr[j-1]++;
                if (linkR(2,ii-1) == jj)
                    Ntr[j-1]++;
            }
        }
        //  Check vector NTR i.e. find out if all 0th order reaches are used (as
        //  inputs), or if a higher order reach is used more than once , or if
        //  more than one higher order reach is not used as an input, or if there
        //  is no outlet to the system i.e. NTR(I).NE.0 for all I.
        //  ******************************************************
        for (i = 1; i <= n0rch; i++) {
            ii = ll[i-1];
            if (Ntr[ii-1] == 0) {
                cerr << "2 WARNING ! Reach number " << dec << setw(3) << linkR(0,ii-1);
                cerr << " (a 0th order reach) is not used as an input in the river network.\n";
            }
        }
        high = false;
        for (i = n0rch+1; i <= Nrch; i++) {
            ii = ll[i-1];
            if (Ntr[ii-1] == 0) {
                // RAW 16/7/2002 wanted to simulate more than one catchment at once so removed
                // the test on HIGH. We will only re-instate if we have problems - RPI
                //	  IF(HIGH) THEN  ! RAW commented this out 16/7/2002
                //	    WRITE(21,555) LINKR(1,II)   ! RAW commented this out 16/7/2002
                //  555 FORMAT(" Reach link data corrupt! Reach no. ",I3," is the second r
                //     1each encountered"/" that does not have a downstream link.")
                //	    GOTO 290  ! RAW commented this out 16/7/2002
                //	  ENDIF   ! RAW commented this out 16/7/2002
                high = true;
            }
            if (Ntr[ii-1] > 1) {
                cerr << "2 WARNING ! Reach number " << dec << setw(3) << linkR(0,ii-1);
                cerr << " (a higher order reach) is used more than\n";
                cerr << " once as an input. This will cause errors in the water balance calculations.\n";
            }
        }
        if (!high) {
            cerr << " Reach link data corrupt! There is no final outlet in theriver network.\n";
            cerr << " The last reach processed was reach no. " << dec << setw(3) << linkR(0,ii-1);
            exit(EXIT_FAILURE);
        }
        // 180 CONTINUE (This was active here in the fortran code, but no active line referrs to it.)
        //	write(6,5090)NGAUGE,NS,NRCH,NKA,ND,LL,NTR,NTS,LINKS,LINKR
        //  Fill in the vector NTS by first filling in the matrix NTSG which
        //  indicates how many times each slope is used for each gutter.
        //  ***********************************************************
        for (i = 1; i <= Ns; i++) {
            Nts[i-1] = 0;
            for (j = 1; j <= n0rch; j++) {
                jj = ll[j-1];
                Ntsg(i,j) = 0;
                if (linkR(1,jj-1) == i)
                    Ntsg(i,j)++;
                if (linkR(2,jj-1) == i)
                    Ntsg(i,j)++;
            }
        }
        for (i = 1; i <= Ns; i++) {
            for (j = 1; j <= n0rch; j++) {
                jj = ll[j-1];
                if (allzer)
                    Ntr[jj-1] = 1;
                Nts[i-1] += Ntsg(i,j)*Ntr[jj-1];
            }
        }
        //  Sort out llout array so that it is an index to ll array
        // Array LL defines the processing order of reaches. Array LLOUT contains
        // the position within LL of reaches with measured values or for which
        // output is required.
        for (i = 1; i <= max(Neq, Nout); i++) {
            ii = kllout(i-1);
            llout[i-1] = 0;
            for (j = 1; j <= Nrch; j++) {
                jj = ll[j-1];
                if (linkR(0,jj-1) == abs(ii)) {
                    if (llout[i-1]  ==  0) {
                        llout[i-1] =  ii < 0 ? -j : j;
                    } else {
                        cerr << "Logic error or duplicate output reach\n";
                    }
                }
            }
        }
    } // nooksack == 0 (the 'else' part of if nooksack == 1
    //  Sequence to read subbasin output numbers 6.1.98
    getline(modelspcFile, inLine, '\n');                // Comment line
    modelspcFile >> nBout;
    getline(modelspcFile, inLine, '\n');	// Read the rest of the line.
    for (i = 0; i < nBout; i++) {
        modelspcFile >> iTbout[i];
    }
    getline(modelspcFile, inLine, '\n');    // Read the rest of the line.
    for (i = 0; i < Ns; i++) {
        iBout[i] = 0;
    }
    for (i = 1; i <= nBout; i++) {
        iBout[iTbout[i-1]-1] = 1;  //  Set flags for basins to output
    }

    //  Resume reading main data disk file
    //  Read basin properties
    //  *********************

    ifstream basinparsFile("basinpars.txt");	// Fortran unit 88
    if (!basinparsFile.is_open()) {
        cerr << "Failed to open basinpars.txt\n";
        exit(EXIT_FAILURE);
    }

    getline(basinparsFile, inLine, '\n');                // Comment line

    for (i = 0; i < Ns; i++) {
        getline(basinparsFile, inLine, '\n');
        char ch = inLine.back();    // Is there a comma at the end?
        // In C++03, std::string::back is not available due to an oversight,
        // but you can get around this by dereferencing the reverse_iterator you get back from rbegin:
        // char ch = *myStr.rbegin();
        if (ch == ',') {
            //cout << "A comma was found at the end of the line." << endl;;
            inLine.pop_back();
        }
        std::istringstream iss(inLine);
        for (j = 0; j < num_basinpars-1; j++) {
            iss >> bp(j,i);
            iss >> comma;
        }
        j = num_basinpars-1;
        iss >> bp(j,i); // no comma.
    }
    basinparsFile.close();
    i = Ns;

    //  Insertion to get latitude and longitude for each basin

    //   Read latlongfromxy.dat
    ifstream latlongfromxyFile("latlongfromxy.txt");	// Fortran unit 88 (again)
    if (!latlongfromxyFile.is_open()) {
        cerr << "Failed to open latlongfromxy.txt\n";
        exit(EXIT_FAILURE);
    }
    getline(latlongfromxyFile, inLine, '\n');			// Comment line
    getline(latlongfromxyFile, inLine, '\n');			// Comment line
    latlongfromxyFile >> flat;
    getline(latlongfromxyFile, inLine, '\n');			// Read the rest of the line.
    latlongfromxyFile >> clat;
    getline(latlongfromxyFile, inLine, '\n');			// Read the rest of the line.
    latlongfromxyFile >> flong;
    getline(latlongfromxyFile, inLine, '\n');			// Read the rest of the line.
    latlongfromxyFile >> clong;
    latlongfromxyFile.close();

    //   calculate basin lat long.
    //   Assumptions
    //   1.  Linear tranformation between local coords and lat long is sufficient
    //   2.  Outlet local coordinates is representative of subbasin
    for (j = 0; j < Ns; j++) {
        bXlat[j] = flat*bp(6,j) + clat;
        bXlon[j] = flong*bp(5,j) + clong;
    }
    //   END  Lat and long modifications

    //   DGT modified below to streamline subbasin property input  3/97
    // See if there is a constraints file

    exist = false;
    ifstream modelconFile;
    if ((dp = opendir(".")) == NULL) {
        cerr << "Can't open directory\n";
        exit(EXIT_FAILURE);
    }
    while ((dirp = readdir(dp)) != NULL) {
        testStr = dirp->d_name;
        if (testStr == "modelcon.dat") {
            exist = true;
            closedir(dp);
            break;
        }
    }

    // RPI 12/8/2002 introduced limitc here and rplaced exist in the argument list
    // because exist is used when opening many files but we need to keep the value
    // associated with modelcon for use elsewhere. It leaves here as limitc, undergoes
    // a name change to exist in input, i.e. we change the name in the calls so as
    // not to have to alter code elsewhere
    limitC = exist;   // exist is local, limitc is glpobal for constraints
    // Open constraints file if it exists
    if (exist == false) {
        cerr << " No 'modelcon.dat' file found\n";
    } else {
        ifstream modelconFile("modelcon.dat");	// unit LUNopt
        if (!modelconFile.is_open()) {
            cerr << "Failed to open modelcon.dat\n";
            exit(EXIT_FAILURE);
        }
    }
    getline(modelspcFile, inLine, '\n');			// Comment line
    getline(modelspcFile, inLine, '\n');			// Comment line
    for (i = 1; i <= Ns; i++) {
        getline(modelspcFile, inLine, '\n');			    // Comment line (Basin #)
        if (exist && i == 1) {                              // Only if there is a modelcon.dat file
            getline(modelconFile, inLine, '\n');			// Comment line in modelcon.dat file
        }
        Line = 4;
        for (j = 1; j <= isp; j++) {
            if (exist && i == 1) {                          // Only if there is a modelcon.dat file
                modelconFile >> minSp[j-1] >> maxSp[j-1];
            }
            modelspcFile >> Sp(j-1,i-1);	//  One property per line
            getline(modelspcFile, inLine, '\n');            // Read the rest of the line.
        }
        for (j = 1; j <= Nsp; j++) {		// overwrite using basinpars
            Sp(j-1,i-1) = bp(j+6,i-1);		// bp(1:7,:) has other stuff we also need
        }

        // READ A/TANB DISTBN
        getline(modelspcFile, inLine, '\n');	// Comment line
        modelspcFile >> Nka[i-1];
        getline(modelspcFile, inLine, '\n');            // Read the rest of the line.
        //	write(6,5090)ngauge,ns,nrch,nka,nd,ll,ntr,nts,links,linkr
        if ( Nka[i-1] > maxA ) {
            cerr << " too many points in ln(a/tanb) distbn, maxA = " << dec << setw(3) << maxA << endl;
            exit(EXIT_FAILURE);
        }
        getline(modelspcFile, inLine, '\n');			// Comment line
        for (j = 1; j <= Nka[i-1]; j++) {
            modelspcFile >> atb[j-1][i-1] >> pka[j-1][i-1];
            getline(modelspcFile, inLine, '\n');        // Carriage return or comment.
        }

        temp = 0.0;
        //	do 4000 j=1,nka(i)
        for (j = 2; j <= Nka[i-1]; j++) { // first one is not used
            //	if ( atb(j,i).lt.0d0.or.pka(j,i).lt.0d0.or.pka(j,i).gt.1d0 )
            if ( pka[j-1][i-1] < 0.0 || pka[j-1][i-1] > 1.0 ) {  // atb < 0 is a valid number
                cerr << " Error in line " << dec << setw(2) << j << " of LOG(A/TANB) for subcat ";
                cerr << dec << setw(3) << i << " PKA(J,I) < 0 OR PKA(J,I) > 1\n";
                exit(EXIT_FAILURE);
            }
            if (j > 1 && atb[j-1][i-1] <= atb[j-2][i-1] ) {
                cerr << " Error in line " << dec << setw(2) << j << " of LOG(A/TANB) for subcat ";
                cerr << dec << setw(3) << i << " : ATB(J,I) <= ATB(J-1,I)\n";
                exit(EXIT_FAILURE);
            }
            if (nooksack == 1) {
                pka[Nka[i-1]-1][i-1] = 0;// nooksack doesn't accept the flat areas
                temp += pka[j-1][i-1];
            }
        }
        if ( nooksack != 1 && fabs(temp - 1.0) > 0.01 ) {
            cerr << " Sum of area proportions for a/tanb <> 1: ";
            cerr << fixed << setw(16) << setprecision(8) << temp << '\n';
            exit(EXIT_FAILURE);
        }

        //  Here temp is the sum of proportions and is close to 1.
        //  Normalize all proportions so that the proportions are numerically
        //  equal to 1 to avoid mass balance errors.
        for (j = 2; j <= Nka[i-1]; j++) {
            pka[j-1][i-1] = pka[j-1][i-1]/temp;
        }
        // areal integral of a/tanb (=lambda)
        // use values at centre of class ranges to calculate mean
        tl[i-1] = atb[0][i-1]*pka[0][i-1];
        for (j = 2; j <= Nka[i-1]; j++) {
            tl[i-1] += pka[j-1][i-1]*(atb[j-1][i-1] + atb[j-2][i-1])/2.0;
        }

        // read overland flow routing distbn
        getline(modelspcFile, inLine, '\n');			// Comment line
        modelspcFile >> Nd[i-1];
        //	write(6,5090)ngauge,ns,nrch,nka,nd,ll,ntr,nts,links,linkr
        if ( Nd[i-1] > maxC ) {
            cerr << " Too many points in channel distbn, MAXC = " << dec << setw(3) << maxC << '\n';
            exit(EXIT_FAILURE);
        }
        getline(modelspcFile, inLine, '\n');            // Read the rest of the Nd line.
        getline(modelspcFile, inLine, '\n');			// Comment line
        for (j = 1; j <= Nd[i-1]; j++) {
            modelspcFile >> cl2[j-1][i-1] >> pd2[j-1][i-1];
        }
        for (j = 1; j <= Nd[i-1]; j++) {
            Rtemp1 = cl2[j-1][i-1];
            Ptemp1 = pd2[j-1][i-1];
            if (j > 1) {
                Rtemp2 = cl2[j-2][i-1];
                Ptemp2 = pd2[j-2][i-1];
            }
            if (Rtemp1 < 0 || Ptemp1 < 0 || Ptemp1 > 1) {
                cerr << " Error in line " << dec << setw(2) << j;
                cerr << " of Flow Dist. for subcat " << dec << setw(3) << i;
                cerr << " : cl(j,i).lt.0 || pd(j,i).lt.0 || pd(j,i).gt.1\n";
                exit(EXIT_FAILURE);
            }
            if (j > 1) {
                if(Rtemp1 < Rtemp2 || Ptemp1 < Ptemp2) {
                    cerr << " Error in line " << dec << setw(2) << j;
                    cerr << " of Flow Dist. for subcat " << dec << setw(3) << i;
                    cerr << " : cl(j,i).lt.cl(j-1,i) || pd(j,i).lt.pd(j-1,i)\n";
                    exit(EXIT_FAILURE);
                }
            }
        }
        Rtemp1 = cl2[0][i-1];
        Ptemp1 = pd2[0][i-1];
        Ptemp2 = pd2[Nd[i-1]-1][i-1];
        if (Rtemp1 != 0 || Ptemp1 != 0 || Ptemp2 != 1) {
            cerr << " Error in Flow Dist. for subcat " << dec << setw(3) << i;
            cerr << " : cl[0][i-1] != 0 || pd[0][i-1] != 0 || pd(last,i) != 1\n";
            exit(EXIT_FAILURE);
        }
        //  Read basin initial properties
        //  **********************************
        //	write(6,7712) LUNRD
        // 7712	format(" LUNRD",i6)
        getline(modelspcFile, inLine, '\n');            // Read the rest of the line.
        getline(modelspcFile, inLine, '\n');			// Comment line
        if (exist && i == 1) {
            getline(modelconFile, inLine, '\n');			// Comment line
            for (j = 1; j <= Nsi; j++) {
                modelconFile >> minSi[j-1] >> maxSi[j-1];
            }
        }
        //      DO 230 I=1,NS
        for (j = 1; j <= Nsi; j++) {
            modelspcFile >> si[j-1][i-1];	// One property per line.
        }
        getline(modelspcFile, inLine, '\n');            // Read the rest of the line.
    }	// i in Ns loop
    modelspcFile.unget();

    // Read in the parameter map. The first number tells if the following data
    // relates to a sub-basin (1), a reach (2), or an initial condition (3)
    // The second value tells which hillslope, reach or initial condition is
    // to be altered. The third indicates which sub-basin the parameter relates to.
    // If the third parameter is .le. 0 then the parameter value is to be assigned
    // to all sub-basins. The fourth value indicates where the parameter lies in the
    // NLFIT parameter box list.

    //   DGT Parameter map read was here - moved.

    // Check the use of negative area values to indicate  basin copies

    for (i = 1; i <= Ns; i++) {
        if ( Sp(0,i-1) <= 0.0 ) {
            ii = lround(Sp(1,i-1));	// fortran nint(SP(2,I))
            if ( ii < 1 || ii > Ns ) {
                cerr << " basin" << dec << setw(3) << i;
                cerr << " area is <= 0 : must have 1 <= length <= ";
                cerr << dec << setw(3) << Ns << '\n';
                exit(EXIT_FAILURE);
            }
            if ( Sp(0,ii-1) < 0.0 ) {
                cerr << " basin" << dec << setw(3) << i;
                cerr << "area is <=0:length must point to real basin\n";
                exit(EXIT_FAILURE);
            }
        }
    }

    //  Read reach properties
    //  *********************
    // Note that because we added artificial gutters, we don't expect to
    // read any data for them, so replace NRCH by NRCH-N0RCH
    modelspcFile.get(cr);	// read line end after last cl2, pd2 read
    getline(modelspcFile, inLine, '\n');			// Comment line

    for (i = 1; i <= Nrch-n0rch; i++) {
        Line = 6;
        if (exist && i == 1) {
            getline(modelconFile, inLine, '\n');			// Comment line
            for (j = 1; j <= Nrp; j++) {
                modelconFile >> minRp[j-1] >> maxRp[j-1];
            }
        }
        for (j = 1; j <= Nrp; j++) {
            modelspcFile >> Rp[j-1][i-1];	// One property per line
        }
    }
    modelspcFile.close();
    if (exist) {
        modelconFile.close();
    }
    //  *************************
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving mdData(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

    return 0;
}


//-----------------------------------------------------------------------------
//         SUBROUTINE READ_LAKES
//-----------------------------------------------------------------------------
int read_lakes(Array<int,Dynamic,1> &lake_reach, int *lzero, double *lake_areas, int *lake_beach_slps, int *lk_line, int *num_rat_vals,
               int **lheads, int **loflows, const int max_lakes, const int max_lheads)
{
    // lake variables
    int i, j;
    bool exist;
    struct dirent *dirp;
    DIR *dp;
    string testStr, inLine;
    char cr;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> read_lakes(" << ncalls << ")" << std::endl;
    }
    caller = "read_lakes";
#endif
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
        cerr << " Found lakes.dat\n";
        ifstream lakesFile("lakes.dat");	// Fortran unit 5
        if (!lakesFile.is_open()) {
            cerr << "Failed to open 'lakes.dat'\n";
            exit(EXIT_FAILURE);
        }

        getline(lakesFile, inLine, '\n');
        lakesFile >> lakes1::nlakes;
        lakesFile.get(cr);	// read line end after last read
        getline(lakesFile, inLine, '\n');
        for (i = 0; i < lakes1::nlakes; i++) {
            lakesFile >> lake_reach(i);
        }
        lakesFile.get(cr);	// read line end after last read
        getline(lakesFile, inLine, '\n');
        for (i = 1; i <= lakes1::nlakes; i++) {
            getline(lakesFile, inLine, '\n');
            lakesFile >> lzero[i-1] >> lake_areas[i-1] >> lake_beach_slps[i-1] >> lk_line[i-1];
            lakesFile.get(cr);	// read line end after last read
            getline(lakesFile, inLine, '\n');
            lakesFile >> num_rat_vals[i-1];
            lakesFile.get(cr);	// read line end after last read
            getline(lakesFile, inLine, '\n');
            for (j = 1; j <= num_rat_vals[i-1]; j++) {
                lakesFile >> lheads[i-1][j-1] >> loflows[i-1][j-1];
            }
        }
        lakesFile.close();
    } else {
        cerr << " No 'lakes.dat' file found\n";
    }
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving read_lakes(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

    return 0;
}

// *************************************************************************
int chk_sites(int &lrg, const int sites, const int j, int &iMatch)
{
    int itol, lpe, lpn, spe, spn, de, dn;
    // itol=1 means we will tolerate an error of 1 part in 1000 in the easting or northing
    // component of the site number for a mesoscale grid point
    itol = 1;
    // we want an exact test for equality if we are using raingauges
    if (lrg <= 9999999)
        itol = 0;
    lpe = lrg/10000;
    lpn = lrg - lpe*10000;
    spe = sites/10000;
    spn = sites - spe*10000;
    de  = abs(lpe - spe);
    dn  = abs(lpn - spn);
    if (de+dn == 0) {
        lrg = j;
        iMatch = 1;
        return 0;
    }
    if (de <= itol && dn <= itol) {
        lrg = j;
        iMatch = 1;
        return 0;
    }

    return 0;
}

// ******************************************************************************
// Removed subroutine variables lzero,lake_areas, lake_beach_slps,lk_line,num_rat_vals,
// lheads,loflows, and max_lheads for C++ version
int read_lakes_levels(const Array<int,Dynamic,1> &lake_reach, double *ini_levels, const int max_lakes)
{
    // lake variables
    bool exist;
    int lake_reach_ini[max_lakes], nlakes_ini;
    int i, j, itemp;
    double temp;
    struct dirent *dirp;
    DIR *dp;
    string testStr, inLine;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> read_lakes_levels(" << ncalls << ")" << std::endl;
    }
    caller = "read_lakes_levels";
#endif
    exist = false;
    if ((dp = opendir(".")) == NULL) {
        cerr << "Can't open directory\n";
        exit(EXIT_FAILURE);
    }
    while ((dirp = readdir(dp)) != NULL) {
        testStr = dirp->d_name;
        if (testStr == "lakes_levels.dat") {
            exist = true;
            closedir(dp);
            break;
        }
    }
    if (exist) {
        cerr << " Found lakes_levels.dat\n";
        ifstream lakes_levelsFile("lakes_levels.dat");	// Fortran unit 5
        if (!lakes_levelsFile.is_open()) {
            cerr << "Failed to open 'lakes_levels.dat'\n";
            exit(EXIT_FAILURE);
        }
        for (i = 1; i <= lakes1::nlakes; i++) {
            ini_levels[i-1] = -1.0;
        }

        getline(lakes_levelsFile, inLine, '\n');
        lakes_levelsFile >> nlakes_ini;
        for (i = 1; i <= nlakes_ini; i++) {
            lakes_levelsFile >> lake_reach_ini[i-1];
        }
        // negative values mean program is to calculate the initial level
        for (i = 1; i <= nlakes_ini; i++) {
            lakes_levelsFile >> ini_levels[i-1];
        }

        // sort the initial levels into the same order as the rest of the lake data
        for (i = 1; i <= lakes1::nlakes; i++) {
            for (j = 1; j <= nlakes_ini; j++) {
                if(lake_reach(i-1) == lake_reach_ini[j-1]) {
                    temp                = ini_levels[i-1];
                    itemp               = lake_reach_ini[i-1];
                    ini_levels[i-1]     = ini_levels[j-1];
                    lake_reach_ini[i-1] = lake_reach_ini[j-1];
                    ini_levels[j-1]     = temp;
                    lake_reach_ini[j-1] = itemp;
                }
            }
        }
        lakes_levelsFile.close();
    } else {
        cerr << " No 'lakes_levels.dat' file found\n";
    }
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving read_lakes_levels(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

    return 0;
}
