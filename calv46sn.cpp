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
#include <iomanip>
#include "snow.hh"

using namespace std;
using namespace Eigen;

ofstream oFile[15];
ofstream snowcontrol3_File;	// unit 10

// This version,V3, has lakes and snow modelling added to it

// ******************************************************
// *  subroutine  calcts
// ******************************************************
int calcts( double **Si,            const ArrayXXd &Sp,           double **Rp,           ArrayXXi &linkR,        const int *ll,
            const int Nsub,         const int *Nka,        const double *tl,   double **atb,
            double **pka,           const int *nd,         double **cl,           double **pd,        const double units,
            const int ipsub,        const int ipatb,       const bool reinit,     const bool modwrt,  const int stim,
            ArrayXXd &bRain, const long int interval, const int m,        const int mi,
            const int mps,          const int mpe,         bool &ok,              const double *xlat, const double *xlong,
            const double stdlon,    const double *elevtg,  double **bdtBar,       const int sDate,    int &sHour,
            const ArrayXd &temper,   const double *dewp,    const double *tRange,  const int Neq,      const int Nout,
            const int nBout,        const int *iBout,      const double *wind2m,  double **bTmin,     double **bTmax,
            double **bTdew,         const double *bXlat,   const double *bXlon,   int *ntdh,    const int ndump,
            const int maxInt,       const int maxSlp,      const int maxA,        const int maxC,     const int maxChn,
            const int idebugoutput, const int idebugbasin, const int idebugcase)
{
    istringstream line;

    // The Model
    int ievap_method, n_irrig, n_drainage, kk, i_irrig, i_drainage, kp;
    double wt0, wt1, wt2, wt12, depth_irrig_dem, rate_irrig, art_drainage;
    double dep, qinst_out_0, dr_out_d, dr_out_0;
    double zbm_d, acsem_d, s1_d, sr_d, cv_d, bal_d, qinst_out_d;
    double zbar_d;
    double sumr_d, sumq_d, sumae_d, s0_d, sumpe_d, dth1;
    double art_drainage_out;   // artificial drainage
    double scalefactor;

    // DGT 5/27/12 allocating arrays that seemed to be allocated dynamically before
    double *quc;
    double *groundwater_to_take; 					// groundwater_to_take(maxSlp)
    double evap_mm, qlat_mm;
    double *evap_for_watermgmt, *precip_for_watermgmt;
    const int nreg = 6;								// how many regions can we model within a sub-basin?
    // use these regions for irrigation and drainage

    //flood comment out next line
    static Array<int,Dynamic,1> irr_l(Nip1);
    static ArrayXd rirr(Nip1);
    static int istep, j, n;

    //   locals
    int js;
    //   ET variables
    static double albedo, elevsb, rlapse;
    static double PETsngl, elev;
    static double temper_temp, dewp_temp, trange_temp;
    // snow
    static double ddf;
    static int irad;

    //         REACH ROUTING ( ROUTE )
    static int jr, jr1;

    // locals:
    static double prec, pet=0.0;
    static double *snowst;

    int jsub, i;
    static int jj;

    static double smin;

    // This is special for CALCTS
    static int  iyear, month = 0, iday, ihr, imm, isec, ihh, iss;
    static double hour1;
    static double **sr, **cv;
    static double *q0, **s0;
    static double **tdh;
    static double **aciem;
    static double **acsem;
    static double **zbm;
    static double **sumr;
    static double **sumae;
    static double **sumpe;
    static double **sumq;
    static double **sumie;
    static double **sumse;
    static double **sumqb;
    static double **sumce;


    static double sumad;
    static double **sumsle;
    static double **sumr1;
    static double **sumqv;

    static double *qb;
    static double *zr;
    static double *ak0fzrdt;
    static double *logoqm;
    static double *qvmin;
    static double *dth;

    static int ngut, natball;

    // SNOWUEB
    const int Nsv = 9, Npar = 25, Nxv = 11, Niv = 7;
    static int ndepletionpoints, ipflag, mQgOption;
    static int nintstep, isurftmpoption, nstepday, isnow_method;
    static double bca, bcc, a, b, timestep;
    double surfacewaterinput;
    static double snowevaporation, areafractionsnow;
    static double ta, p, ws, qsiobs, qnetob;
    static string inFile, outFile, pFile, svFile, bcFile;
    static string aFile, dfcFile, outFileDet;
    static string inLine;
    static ArrayXd snowsitev(Nsv);
    //double **snowstatev;

    static ArrayXd snowforcing(Niv);
    static ArrayXd snowparam(Npar);
    static Array<int,Dynamic,1> snowcontrol(7);
    static ArrayXd dtbar(12);
    static double *cump;
    static double *cume;
    static double *cummr;
    static double *errmbal;
    static double *w1;

    // arrays to keep records of surface temperature and snowpack average
    // temperature. This is for the fourth model (Modified force restore approach)
    static double **dfc;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> calcts(" << ncalls << ")" << std::endl;
    }
    caller = "calcts";
#endif

    // Allocate and initialize all the now dynamic arrays
    ArrayXd         baseflow(maxSlp);
    ArrayXd      totalrunoff(maxSlp);
    ArrayXd    potentialevap(maxSlp);
    ArrayXd          tempave(maxSlp);
    ArrayXd           surfro(maxSlp);
    ArrayXd         canstore(maxSlp);
    ArrayXd        soilstore(maxSlp);
    ArrayXd            tiled(maxSlp);
    ArrayXd           ditchd(maxSlp);
    ArrayXd      ArtDrainage(maxSlp);
    ArrayXd volume_irrig_sup(maxSlp);
    ArrayXd vol_irrig_demand(maxSlp);
    baseflow         = 0.0;
    totalrunoff      = 0.0;
    potentialevap    = 0.0;
    tempave          = 0.0;
    surfro           = 0.0;
    canstore         = 0.0;
    soilstore        = 0.0;
    tiled            = 0.0;
    ditchd           = 0.0;
    ArtDrainage      = 0.0;
    volume_irrig_sup = 0.0;
    vol_irrig_demand = 0.0;
    groundwater_to_take  = new double[maxSlp];
    evap_for_watermgmt   = new double[maxSlp];
    precip_for_watermgmt = new double[maxSlp];
    for (i = 0; i < maxSlp; i++) {
        groundwater_to_take[i]  = 0.0;
        evap_for_watermgmt[i]   = 0.0;
        precip_for_watermgmt[i] = 0.0;
    }
    quc = new double[maxInt];
    for (i = 0; i < maxInt; i++) {
        quc[i] = 0.0;
    }
    ArrayXXd bdtBarR(12,maxSlp);
    ArrayXXd zbar(maxSlp,nreg);
    snowst   = new double[maxSlp];
    q0       = new double[maxSlp];
    qb       = new double[maxSlp];
    qvmin    = new double[maxSlp];
    zr       = new double[maxSlp];
    ak0fzrdt = new double[maxSlp];
    logoqm   = new double[maxSlp];
    dth      = new double[maxSlp];
    cump     = new double[maxSlp];
    cume     = new double[maxSlp];
    cummr    = new double[maxSlp];
    errmbal  = new double[maxSlp];
    w1       = new double[maxSlp];
    sr     = new double*[maxSlp];
    cv     = new double*[maxSlp];
    s0     = new double*[maxSlp];
    aciem  = new double*[maxSlp];
    acsem  = new double*[maxSlp];
    zbm    = new double*[maxSlp];
    sumr   = new double*[maxSlp];
    sumae  = new double*[maxSlp];
    sumpe  = new double*[maxSlp];
    sumq   = new double*[maxSlp];
    sumie  = new double*[maxSlp];
    sumse  = new double*[maxSlp];
    sumqb  = new double*[maxSlp];
    sumce  = new double*[maxSlp];
    sumsle = new double*[maxSlp];
    sumr1  = new double*[maxSlp];
    sumqv  = new double*[maxSlp];
    for (j = 0; j < maxSlp; j++) {
        sr[j]     = new double[nreg];
        cv[j]     = new double[nreg];
        s0[j]     = new double[nreg];
        aciem[j]  = new double[nreg];
        acsem[j]  = new double[nreg];
        zbm[j]    = new double[nreg];
        sumr[j]   = new double[nreg];
        sumae[j]  = new double[nreg];
        sumpe[j]  = new double[nreg];
        sumq[j]   = new double[nreg];
        sumie[j]  = new double[nreg];
        sumse[j]  = new double[nreg];
        sumqb[j]  = new double[nreg];
        sumce[j]  = new double[nreg];
        sumsle[j] = new double[nreg];
        sumr1[j]  = new double[nreg];
        sumqv[j]  = new double[nreg];
    }

    double ***dr = new double**[maxSlp];
    double ***qinst = new double**[maxSlp];
    for (js = 0; js < maxSlp; ++js) {
        dr[js] = new double*[nreg];
        qinst[js] = new double*[nreg];
        for (kk = 0; kk < nreg; ++kk) {
            dr[js][kk] = new double[MAX_NTDH];
            qinst[js][kk] = new double[MAX_NTDH];
        }
    }
    ArrayXd dr1(MAX_NTDH);
    ArrayXd qinst1(MAX_NTDH);

    double***zbar_in = new double**[m+1];
    for (int i = 0; i <= m; ++i) {
        zbar_in[i] = new double*[Nsub];
        for (int j = 0; j < Nsub; ++j) {
            zbar_in[i][j] = new double[nreg];
        }
    }

    //Matrix<double,Dynamic,3,3>    dr(MAX_NTDH, nreg, maxSlp);
    //Matrix<double,Dynamic,3,3> qinst(MAX_NTDH, nreg, maxSlp);

    tdh = new double*[MAX_NTDH];
    for (j = 0; j <= MAX_NTDH; j++) {
        tdh[j] = new double[maxSlp];
    }
    ArrayXXd snowstatev(Nxv, maxSlp);

    // *******************************************************************************
    // THIS SECTION IS RELEVANT TO THE WHOLE MODEL

    // ASSUME EVERYTHING WILL TURN OUT OK
    ok = true;

    //  Calculate sum of catchment areas above each reach
    ngut = Nsub;
    ievap_method = 2; //0=P-T, 1 and 2 are P-M
    // initialise snowueb model
    aFile = "snow.in";
    ifstream snowFile(aFile.c_str());	// unit 1
    if (!snowFile.is_open()) {
        cerr << "Failed to open " << aFile << '\n';
        exit(EXIT_FAILURE);
    }
    snowFile >> inFile >> outFile >> outFileDet >> pFile >> svFile >> bcFile >> dfcFile;
    outFile = "results/" + outFile;
    snowFile >> irad;
    getline(snowFile, inLine, '\n');           // Read the remainder of the line
    //   Flag to control radiation
    //    0 is no measurements - radiation estimated from diurnal temperature range
    //    1 is incoming shortwave radiation read from file (measured), incoming longwave estimated
    //    2 is incoming shortwave and longwave radiation read from file (measured)
    //    3 is net radiation read from file (measured)
    snowFile >> ipflag;  						// Flag to control printing (1=print)
    getline(snowFile, inLine, '\n');			// Read the remainder of the line
    snowFile >> nintstep;  						// Number of internal time steps to use
    getline(snowFile, inLine, '\n');			// Read the remainder of the line
    snowFile >> isurftmpoption; 				// Surface temperature algorithm option
    getline(snowFile, inLine, '\n');			// Read the remainder of the line
    snowFile >> mQgOption;  					// ground heat input
    getline(snowFile, inLine, '\n');			// Read the remainder of the line
    snowFile.close();

    ifstream snowparamFile(pFile.c_str());	// unit 88
    if (!snowparamFile.is_open()) {
        cerr << "Failed to open " << pFile << '\n';
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < Npar; i++) {
        snowparamFile >> snowparam(i);
        getline(snowparamFile, inLine, '\n');			// Read the remainder of the line
    }
    snowparamFile.close();

    bcparm(dtbar, bca, bcc, bcFile);	// warning we wont use the dtbar values we read here but use this read to get bca and bcc
    // dtbar values used are read in hyData and averaged for nearby gages
    ifstream dcFile(dfcFile.c_str());	// unit 88
    if (!dcFile.is_open()) {
        cerr << "Failed to open " << dfcFile << '\n';
        exit(EXIT_FAILURE);
    }
    i = 0;
    getline(dcFile, inLine, '\n');			// Read the first line.
    while (inLine.length() > 1) {
        line.str(inLine);
        line >> a >> b;
        i++;
        getline(dcFile, inLine, '\n');
    }
    ndepletionpoints = i;
    dfc = new double*[ndepletionpoints];
    for (j = 0; j < ndepletionpoints; j++) {
        dfc[j] = new double[2];
    }
    // rewind
    dcFile.clear();
    dcFile.seekg(0);

    for (i = 0; i < ndepletionpoints; i++) {
        dcFile >> dfc[i][0] >> dfc[i][1];
    }
    dcFile.close();

    for (i = 0; i < Nxv; i++) {
        for (j = 0; j < maxSlp; j++) {
            snowstatev(i,j) = 0.0;
        }
    }
    //	month,day,year,hour,
    timestep = 24;

    snowcontrol(0) = irad;
    snowcontrol(1) = ipflag;
    snowcontrol(2) = 10;
    snowcontrol(3) = nintstep;

    snowcontrol(6) = isurftmpoption;

    nstepday = 24/timestep*nintstep;
    ArrayXXd snowsurfacetemp(nstepday, maxSlp);
    ArrayXXd snowaveragetemp(nstepday, maxSlp);
    snowsurfacetemp = -9999.0;
    snowaveragetemp = -9999.0;

    for (i = 0; i < maxSlp; i++) {
        snowsurfacetemp(nstepday-1,i) = 0.0;
        snowaveragetemp(nstepday-1,i) = 0.0;
    }

    if (snowcontrol(1) >= 1) {
        snowcontrol3_File.open(outFile.c_str());
        if (!snowcontrol3_File.is_open()) {
            cerr << "Failed to open " << outFile << '\n';
            exit(EXIT_FAILURE);
        }
    }

    //   initialize variables for mass balance
    for (i = 0; i < maxSlp; i++) {
        //w1[i] = snowstatev[1][i];
        w1[i] = snowstatev(1,i);
        cump[i]    = 0.0;
        cume[i]    = 0.0;
        cummr[i]   = 0.0;
        errmbal[i] = 0.0;
    }

    // ***********************************************************************

    // SUBCATCHMENT SECTION

    // Landcare
    natball = 0;
    // implementing multiple raingauges for each sub-basin
    // (see also top.f, mddatav4.f)
    for (js = 0; js < Nsub; js++) {		// for each subbasin
        // Landcare
        natball += Nka[js];
    }		// Now pass 'bRain' instead of 'rain' and use row js of bRain, rather than row links(2,js) of rain
    //	do 3551 it=1,m
    // 3551	write(21,*)(bRain(js,it),js=1,Nsub)

    // write subbasin output headers
    if ( modwrt ) {
        if (idebugoutput >= 1) {
            lunmodFile << nBout << " " << m << '\n';
            lunmodFile << " Basin";
            lunmodFile << " TimeStep IrrDrainCat Afrac SWInput_mm Qlat_mm Qtot_mm Qb_mm";
            lunmodFile << " Recharge_mm SatEx_mm InfEx_mm SurfRo_mm SatAfrac InfAfrac";
            lunmodFile << " IntStore_mm WTDepth_mm SoilStore_mm Pet_mm Aet_mm";
            lunmodFile << " Irrig_mm GWTake_mm IrrDem_mm Prec_mm SWE_mm Sublim_mm";
            lunmodFile << " Tave_C Tdew_C Trange_C ErrClosure_mm\n";

            // Landcare
            if (mpe == -1) {
                lundatFile << Nsub << '\n';
            }
            luntopFile << "Measured and modeled flows\n";
            luntopFile << "Timestep (column 1), Measured reaches (columns 2 to " << dec << setw(3) << Neq+1;
            luntopFile << ") Modeled reaches (columns " << dec << setw(4) << Neq+2 << " to " << dec << setw(4) << Neq+Nout+1 << ")\n";
            luntopFile << "units ARE um/interval normalized by each basins own area";
            luntopFile << " unless they have a -ve site No. in which case they are lake";
            luntopFile << " levels in metres\n";
            luntopFile << "The next two rows give the number of flows to be used for fitting";
            luntopFile << " and printing, followed by the reach numbers to which they relate\n";
        }
    }
    // get any snow input parameters
    ddf = -1.0;
    ifstream snowinpFile("snowinp.txt");
    if (!snowinpFile.is_open()) {
        cerr << "Failed to open 'snowinp.txt'\n";
    }
    else {
        snowinpFile >> ddf;
        snowinpFile.close();
    }

    ofstream snowoutFile("snowout.txt");	// unit 78

    // ***********************************************************************

    // INITIALISE FOR REACH ROUTING

    // The routing cannot handle 0 slopes so fix those.
    // The relationship used is
    //  S = C A^theta = .133351 (A[m^2])^(-.425) = 47.31513 (A[mm^2])^(-.425)
    //  This was fit as a "lower bound" to the scatter in a slope vs area
    //  plot for the Grey river, New Zealand.
    ofstream testnFile("testn_v7.ed");	// unit 78 for testing only
    for (jr1 = 1; jr1 <= model1::Nrch; jr1++) {
        jr = ll[jr1+ngut-1];
        // smin=47.31513 * area(jr)**(-0.425)    Area not valid due to reachlogic above inoperable
        Rp[0][jr-1] = max(smin, Rp[0][jr-1]);
    }
    for (i = 0; i < maxSlp; i++) {
        snowst[i] = 0.0;
    }
#ifdef ZBAR_IN
    double t0 = (double)clock()/(double)CLOCKS_PER_SEC;
    char cr;
    ifstream zbarInFile("results/zbar_in.dat");
    if (!zbarInFile) {
        cerr << "Failed to open zbar_in.dat\n";
        exit(1);
    }
    int in_istep, in_basin;
    for (istep = 0; istep <= m; istep++) {
        for (jsub = 0; jsub < Nsub; jsub++) {
            zbarInFile >> in_istep >> in_basin; // discard the first two integers.
            for (n = 0; n < nreg; n++) {
                zbarInFile >> zbar_in[istep][jsub][n];
            }
            zbarInFile.get(cr);	// read line end
        }
    }
    zbarInFile.close();
    /*for (istep = 0; istep <= m; istep++) {
        for (jsub = 1; jsub <= Nsub; jsub++) {
            cout << dec << setw(4) << istep << dec << setw(3) << jsub << " ";
            for (n = 0; n < nreg; n++) {
                cout << fixed << setw(15) << setprecision(9) << zbar_in[istep][jsub-1][n];
            }
            cout << '\n';
        }
    }
    exit(0); */
    double t1 = (double)clock()/(double)CLOCKS_PER_SEC;
    cout << t1 - t0 << " seconds to read depth to water file \n";
#endif
    // start of time loop +++++++++++++++++++++++++++++++++++++++++
    cout << "Starting time loop\n";
    for (istep = 0; istep <= m; istep++) {
        if (istep%100 == 0) {
            cout << "timestep " << istep << '\n';
        }
#if TRACE
        static int nsteps = 0;
        if (nsteps < MAX_TRACE) {
            traceFile << setw(10) << "timestep " << istep << "\n\n";
        }
        nsteps++;
#endif
        //   LOOP OVER SUBBASINS
        // Work through basins in "stream order" order, i.e.,
        // do basins that feed a first order channels before second order and so on.
        // The order is already defined by array LL, and linkR(2,:) gives the basin number.
        // The first ngut rows in linkR refer to basins (basin # are transformed by
        // having ngut added to them). The next nrch rows of linkR tell which reach
        // is fed by which basin and/or upstream reach. We want the basin numbers
        // from column 2, rows ngut+1 to ngut+nrch in linkR and we want them in the order they
        // should be processed in. The order is in LL(ngut+1) to LL(ngut+nrch). DGT
        // introduced an additive factor of MAXCHN to avoid confusion between basin
        // and reach numbers and we have to remove it.

        for (jsub = 1; jsub <= Nsub; jsub++) {
            js = jsub;
            prec = bRain(js-1,max(istep, 1)-1);
            jr = ll[js-1];

            // snow
            if (istep > 0) {

                //        Subroutine to compute ET  DGT 17 May 1998
                // parameters for ET - because at least albedo may be fitted, it is not possible
                // to put this loop outside calcts - therefore leave as is for the time being RPI 18/3/2004
                // debugging.  To debug a particular watershed uncomment the lines below and enter its topnetID
                //      if (jsub == 87 .and. istep == 245)then
                //	  js=jsub

                //	snowcontrol(1)=1
                //	else
                //	snowcontrol(1)=0
                //	endif
                albedo = Sp(11,js-1);
                rlapse = Sp(12,js-1);
                elevsb = Sp(13,js-1);
                temper_temp = temper(istep-1);
                dewp_temp = dewp[istep-1];
                trange_temp = tRange[istep-1];
                if (ievap_method != 0)
                    dewp_temp = bTdew[js-1][istep-1]; // for use with penman-Monteith
                if (ievap_method != 0)
                    trange_temp = min(30.0, max(0.0, bTmax[js-1][istep-1] - bTmin[js-1][istep-1])); // check for bad data
                if (ievap_method != 0)
                    temper_temp = (bTmax[js-1][istep-1] + bTmin[js-1][istep-1])/2;
                for (i = 0; i < 12; i++) {
                    for (j = 0; j < maxSlp; j++) {
                        bdtBarR(i,j) = bdtBar[i][j];
                    }
                }
                sHour = 0;
                //  Warning here - this code would have to be changed for time steps other than daily
                //  In general it would be better to inherit sHour from the calling program but coming in at 240000
                //  the integration of radiation across a day fails.  A more general solution to this issue would involve
                //  reprogramming hyri to handle time steps that cross the day break.
                //  Setting sHour=0 also achieves compatibility with the convention that inputs are associated with
                //  measurements recorded at the end of the time step (240000 for daily) but ET and snow computations
                //  integrate from the start time over the time step.
                elev = elevsb; // we are driving this using a basin average temperature, so it must be at basin average elev
                ArrayXd bdtBarR1 = bdtBarR.col(js-1);
                //cout << bdtBarR1.size() << endl;exit(0); // column is size 12
                etall(bXlat[js-1], bXlon[js-1], stdlon, elev, bdtBarR1, PETsngl, temper_temp, dewp_temp,
                      trange_temp, elevsb, albedo, rlapse, sDate, sHour, interval, m, istep, iyear, month, iday, ihr,
                      imm, isec, hour1, bTmin[js-1][istep-1], bTmax[js-1][istep-1], wind2m[istep-1], ievap_method);
                pet = PETsngl;
                ta= (bTmin[js-1][istep-1] + bTmax[js-1][istep-1])/2.0;
                p = bRain(js-1,istep-1)*3600.0/interval/1000.0; //mm/int *3600s/h / (s/int) / (1000mm/m) = m/hr for snowueb
                ws = wind2m[istep-1];
                //		RH=1 !unknown for Nooksack, hope we don't need it!
                qsiobs = 0; // unknown for Nooksack, hope we don't need it!
                qnetob = 0; // unknown for Nooksack, hope we don't need it!
                isnow_method = 2;
                if (isnow_method == 2) {
                    snowcontrol(4) = iyear*10000 + month*100 + iday;
                    ihh = hour1;
                    imm = (hour1-ihh)*60;
                    iss = hour1*3600-ihh*3600-imm*60;
                    snowcontrol(5) = ihh*10000 + imm*100*iss;
                    snowforcing(0) = ta;
                    snowforcing(1) = p;
                    snowforcing(2) = ws;
                    snowforcing(3) = dewp_temp;
                    //	snowforcing(4) = 237.3/(1/(log(RH)/17.27+Ta/(Ta+237.3))-1)  ! from http://williams.best.vwh.net/ftp/avsig/avform.txt
                    snowforcing(4) = bTmax[js-1][istep-1] - bTmin[js-1][istep-1];
                    snowforcing(5) = qsiobs;
                    snowforcing(6) = qnetob;
                    snowsitev(0) = 0;   //  Sp(39,js) !0. !forest cover   DGT 4/1/05
                    snowsitev(1) = Sp(13,js-1); //elevsb
                    snowsitev(2) = bXlat[js-1]; //lat
                    snowsitev(3) = bXlon[js-1]; //lon
                    snowsitev(4) = stdlon;  //stdlong
                    snowsitev(5) = Sp(13,js-1); //elevtg is assumed  =  elevsb
                    snowsitev(6) = Sp(12,js-1); //rlapse
                    snowsitev(7) = 0.0; //slope
                    snowsitev(8) = 0.0; //azimuth
                    ArrayXd snowstatev1 = snowstatev.col(js-1);	// a reference to a slice
                    ArrayXd snowsurfacetemp1 = snowsurfacetemp.col(js-1);
                    ArrayXd snowaveragetemp1 = snowaveragetemp.col(js-1);
                    snowueb(snowsitev, snowstatev1, snowparam, ndepletionpoints, dfc, snowcontrol, bdtBarR1,
                            snowforcing, snowsurfacetemp1, snowaveragetemp1, timestep, nstepday,
                            surfacewaterinput, snowevaporation,  // outputs (both in m/h)
                            areafractionsnow, js);   // added js so that within snow one knows which element one is in for debugging

                    snowstatev.col(js-1) = snowstatev1;	// a reference to a slice
                    snowsurfacetemp.col(js-1) = snowsurfacetemp1;
                    snowaveragetemp.col(js-1) = snowaveragetemp1;
                    cump[js-1]  = cump[js-1]  + p*timestep;
                    cume[js-1]  = cume[js-1]  + snowevaporation*timestep;
                    cummr[js-1] = cummr[js-1] + surfacewaterinput*timestep;
                    errmbal[js-1] = w1[js-1] + cump[js-1] - cummr[js-1] - cume[js-1] - snowstatev(1,js-1);
                    if (snowcontrol(1) >= 2) {
                        snowcontrol3_File << dec << setw(2) << js;
                        snowcontrol3_File << dec << setw(5) << istep;
                        snowcontrol3_File << dec << setw(5) << iyear;
                        snowcontrol3_File << dec << setw(3) << month;
                        snowcontrol3_File << dec << setw(3) << iday;
                        snowcontrol3_File << fixed << setw(5)  << setprecision(1) << hour1;
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << surfacewaterinput;
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << snowevaporation;
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << areafractionsnow;
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << snowstatev(1,js-1);
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << cump[js-1];
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << cume[js-1];
                        snowcontrol3_File << fixed << setw(18) << setprecision(9) << cummr[js-1] << '\n';	//swe
                    }
                }
                else {
                    if (ddf >= 0) {
                        snow(snowoutFile, temper, elevtg[0], elevsb, rlapse, prec, ddf, snowst[js-1], interval, Nsub, m, js, istep, maxSlp, maxInt);
                    }
                }
                if ( modwrt && istep == m ) {
                    if (idebugoutput >= 1) {  // DGT 8/17/05 debugout
                        lunpFile << "\n\n Output for sub-catchment " << dec << setw(3) << js << '\n';
                        lunpFile << " ---------------------------\n";
                    }
                }
            }

            prec = surfacewaterinput*1000.0*interval/3600.0;   //mm/timestep=m/h*1000mm/m*h/timestep
            if (nooksack) {
                n_irrig = 1;   // DGT warning.  If this is changed the logic associated with
                //                averaging depth_irrig_dem will need to be adjusted
                n_drainage = 2;
            }
            else {
                n_irrig    = 0;
                n_drainage = 0;
            }

            wt0 = 1.0;
            kk = 0;
            depth_irrig_dem          = 0.0;
            //	xIRR(:)=0
            baseflow(js-1)           = 0.0;
            totalrunoff(js-1)        = 0.0;
            evap_for_watermgmt[js-1] = 0.0;
            potentialevap(js-1)      = 0.0;   // DGT 10/21/12 initialize for averaging over categories
            surfro(js-1)             = 0.0;
            canstore(js-1)           = 0.0;
            soilstore(js-1)          = 0.0;
            tiled(js-1)              = 0.0;
            ditchd(js-1)             = 0.0;

            zbm_d            = 0;
            acsem_d          = 0;
            s1_d             = 0;
            zbar_d           = 0.0;
            sr_d             = 0;
            cv_d             = 0;
            bal_d            = 0;
            sumr_d           = 0;
            sumq_d           = 0;
            sumae_d          = 0;
            s0_d             = 0;
            sumpe_d          = 0;
            qinst_out_d      = 0;
            dr_out_d         = 0;		// Initializing
            art_drainage_out = 0.0;		// Artificial drainage

#ifdef ZBAR_OUT
            zbarFile << dec << setw(4) << istep;
            zbarFile << dec << setw(3) << js;
            // all nreg regions have the same value at this point in the program
#endif
//#ifndef ZBAR_IN
            // the previous timestep's call to watermgmt calculated groundwater_to_take
            if (istep > 1) { //can't do any pumping on first timestep because we haven't called watermgmt yet: assume zero pumping when istep=1
                dth1 = Sp(3,js-1);  // dgt 11/4/07 added this line to get dth1 from corresponding basin parameter
                for (n = 0; n < nreg; n++) {
                    zbar(js-1,n) += groundwater_to_take[js-1]/(Sp(0,js-1)/1.0e6)/dth1;     //RAW 18-Jul-2005 bug fix: Sp(1,JS)/1e6 was Sp(1,JS)*1e6
                }
                //  m/timestep = m3/timestep/(m^2)/porosity
                //   DGT 11/4/07 changed the above from DTH to DTH1 because topmodel saturated zone works with that.
            }
//#endif
            for (i_irrig = 1; i_irrig <= (n_irrig+1); i_irrig++) {
                if (i_irrig == 1) {
                    wt1 = wt0*(1.0 - Sp(19,js-1)); //unirrigated fraction
                    rate_irrig = 0;
                }
                else {
                    wt1 = wt0*Sp(19,js-1); // irrigated fraction
                    if (istep == 0 || wt1 <= 0.0) {
                        rate_irrig = 0;
                    }
                    else {
                        //use volume_irrig_sup[js-1] calculated in previous timestep
                        //                       mm^3        /      mm^2      / (sec) = mm/s
                        //   DGT 8/18/05.  Commented out the /float(interval).  The units going to topmod need to be mm/ts, the same as precipitation
                        //   DGT 8/18/05.                 mm^3                 /      mm^2 = mm/ts
                        rate_irrig = volume_irrig_sup(js-1)*1.0e9/(wt1*Sp(0,js-1));
                        //   /(float(interval))
                    }
                }
                for (i_drainage = 1; i_drainage <= (n_drainage+1); i_drainage++) {
                    kp = (i_irrig - 1)*3 + i_drainage;   // DGT 11/3/07
                    // kp is a variable to index the combinations of irrigation and drainage to compare with kcase and control
                    // detail output.  There are 3 drainage cases.  The first the associated with i_irrig will be numbered 1,2,3
                    // The next three will be numbered 4, 5, 6 and so on.
                    if (i_drainage == 1) { //Naturally drained
                        wt2 = 1.0 - Sp(15,js-1) - Sp(16,js-1); //16=tilefraction, 17=ditchfraction
                        art_drainage = 0;
                    }
                    else if (i_drainage == 2) {     //Tile drained
                        wt2 = Sp(15,js-1); //16=tilefraction
                        art_drainage = Sp(17,js-1); //watch for units
                    }
                    else if (i_drainage == 3) {      //Ditch drained
                        wt2 = Sp(16,js-1); //17=ditchfraction
                        art_drainage = Sp(18,js-1);
                    }
                    wt12 = wt1*wt2; //wt12 is fraction of basin covered by this drainage-irrigation combo
                    kk++;
                    if (wt12 > 0) {
                        for ( i = 0; i < MAX_NTDH; ++i ) {
                            dr1(i) = dr[js-1][kk-1][i];
                            qinst1(i) = qinst[js-1][kk-1][i];
                        }

                        topmod(Si, Sp, js, Nka, tl[js-1], atb, pka, nd, cl, pd, units, irr_l, modwrt,
                               ipsub, ipatb, stim, prec, pet, interval, art_drainage, rate_irrig, month,
                               m, mps, mpe, qinst_out_0,  dr_out_0, ndump, ntdh, istep, maxC, zbm[js-1][kk-1],
                               maxA, maxSlp, maxInt, sumr[js-1][kk-1], sumq[js-1][kk-1], sumae[js-1][kk-1], s0[js-1][kk-1], q0[js-1],
                               sr[js-1][kk-1], cv[js-1][kk-1], aciem[js-1][kk-1], acsem[js-1][kk-1], sumpe[js-1][kk-1], sumie[js-1][kk-1],
                               sumqb[js-1][kk-1], sumce[js-1][kk-1], sumsle[js-1][kk-1], sumr1[js-1][kk-1], qb[js-1], qinst1,
                               dr1, sumqv[js-1][kk-1], sumse[js-1][kk-1], zbar(js-1,kk-1), zbar_in[istep][js-1][kk-1], tdh, zr[js-1], ak0fzrdt[js-1],
                               logoqm[js-1], qvmin[js-1], dth[js-1], sumad, evap_mm, qlat_mm, ipflag, rirr, js);

                        for ( i = 0; i < MAX_NTDH; ++i ) {
                            dr[js-1][kk-1][i] = dr1[i];
                            qinst[js-1][kk-1][i] = qinst1[i];
                        }
                        //  DGT 11/2/07 added qlat and ipflag.  ipflag for debugging.  qlat to return lateral outflows without
                        //     overland flow delays to facilitate mass balance checks
                        //  Added subscript (JS) to Q0 in the calling statement above so that Q0 differences
                        //  between subbasins are preserved
                        //  Added sumad to retrieve artificial drainage calcs

                        if (i_irrig > 1) {
                            irrigation(sr[js-1][kk-1], Sp, js, interval, maxSlp, dep); //watch for units!!!
                            depth_irrig_dem = depth_irrig_dem + wt2*dep; //dep is in m/timestep
                            // DGT 8/20/05 changed the weight in the above  from wt12 to wt2 because depth only should be weighted
                            // among the different drainage classes within irrigated area
                        }
                        else {
                            dep = 0.0;
                        }
                        zbm_d                    += wt12*zbm[js-1][kk-1];
                        acsem_d                  += wt12*acsem[js-1][kk-1];
                        zbar_d                   += wt12*zbar(js-1,kk-1);
                        sr_d                     += wt12*sr[js-1][kk-1];
                        cv_d                     += wt12*cv[js-1][kk-1];
                        sumr_d                   += wt12*sumr[js-1][kk-1];
                        sumq_d                   += wt12*sumq[js-1][kk-1];
                        sumae_d                  += wt12*sumae[js-1][kk-1];
                        sumpe_d                  += wt12*sumpe[js-1][kk-1];
                        s0_d                     += wt12*s0[js-1][kk-1];
                        qinst_out_d              += wt12*qinst_out_0;
                        dr_out_d                 += wt12*dr_out_0;    // 6/10/05  DGT keeping track of total runoff
                        art_drainage_out         += wt12*sumad;  // 6/28/05   DGT Artificial drainage
                        evap_for_watermgmt[js-1] += wt12*evap_mm;   // DGT 8/17/05 keeping track of evap
                        potentialevap(js-1)      += wt12*rirr(12);  // DGT 10/21/12  keeping track of pet
                        surfro(js-1)             += wt12*rirr(6)/1000.0*(Sp(0,js-1)/1.0e6)/interval;
                        canstore(js-1)           += wt12*rirr(9);
                        soilstore(js-1)          += wt12*rirr(11);
                        if (i_drainage  == 2) {
                            tiled(js-1) += wt12*sumad*1.0e-9; // sumad is mm^3/s - now in m^3/s
                        }
                        if (i_drainage == 3) {
                            ditchd(js-1) += wt12*sumad*1.0e-9;
                        }

                        //  Detail topsbd writes
                        if (idebugoutput >= 1) {
                            if (modwrt && iBout[js-1] == 1) {
                                if (istep >= mi && istep > 0) {
                                    if (js == idebugbasin || idebugbasin == 0) {
                                        if (kp == idebugcase || idebugcase == 0) {
                                            lunmodFile << dec << setw(5) << js << dec << setw(8) << istep;
                                            lunmodFile << " " << dec << setw(3) << kp;
                                            lunmodFile << fixed << setw(11) << setprecision(7) << wt12;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << prec;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << qlat_mm;
                                            for (jj = 2; jj <= 14; jj++) {
                                                lunmodFile << " " << fixed << setw(11) << setprecision(6) << rirr(jj-1);
                                            }
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << rate_irrig;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << groundwater_to_take[js-1]/(Sp(0,js-1)/1.0e6)*1000.0;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << dep*1000.0;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << bRain(js-1,istep-1);
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << snowstatev(1,js-1)*1000.0;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << snowevaporation*timestep*1000.0;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << (bTmax[js-1][istep-1] + bTmin[js-1][istep-1])/2.0;
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << bTdew[js-1][istep-1];
                                            lunmodFile << " " << fixed << setw(11) << setprecision(6) << bTmax[js-1][istep-1] - bTmin[js-1][istep-1];
                                            lunmodFile << " " << fixed << setw(14) << setprecision(8) << rirr(14) << '\n';
                                        }
                                    }
                                }
                            }
                        }
                    }
#ifdef ZBAR_IN
                    else {      // Added for reading zbar from MODFLOW
                        zbar(js-1,kk-1) = zbar_in[istep][js-1][kk-1];
                    }
#endif
                } // i_drainage
            } // i_irrig

            //  The one line below is a substantive change DGT 10/13/12.  It amounts to an assumption of
            //   lumped depth to water table, rather than separate depth to water table for each drainage and
            //   irrigation component.  It was made to fix the issue of groundwater levels declining indefinitely
            //   due to there being no recharge from artificially drained areas
            for (n = 0; n < nreg; n++) {
                zbar(js-1,n) = zbar_d;
#ifdef ZBAR_OUT
                zbarFile << scientific << setw(16) << setprecision(9) << zbar(js-1,n); // all nreg regions now have been set to a single value
#endif
            }
#ifdef ZBAR_OUT
            zbarFile << '\n';
#endif
            if (istep > 0)
                tempave(js-1) = (bTmax[js-1][istep-1] + bTmin[js-1][istep-1])/2.0; // DGT 10/21/12 record ave temperature for output

            //accumulate the fluxes from each part of the sub-basin as m3/timestep
            baseflow(js-1) = qinst_out_d*interval*1.0e-9; // m^3/ts = mm^3/s *s/ts *m3/mm3

            //   Use Total Runoff from the floating point DR variable rather than the integer IRR topmodel output
            //    to avoid some divide by zero's later on in water management
            totalrunoff(js-1) = dr_out_d*interval*1.0e-9; // m^3/ts = mm^3/s *s/ts *m3/mm3
            ArtDrainage(js-1) = art_drainage_out*interval*1.0e-9 ;  // DGT 6/28/05  Keep in same units as totalrunoff

            if (modwrt && istep == m) {
                if (idebugoutput >= 1) {  // DGT 8/17/05 debugout
                    if (modwrt) {
                        lunpFile << "MIN ZBAR=";
                        lunpFile << fixed << setw(22) << setprecision(17) << zbm_d << " MAX CONT AREA=" << acsem_d << '\n';
                    }
                    if (modwrt) {
                        lunpFile << "LAMBDA=";
                        lunpFile << fixed << setw(21) << setprecision(15) << tl[js-1] << '\n';
                    }
                }
                //  Storage at the end
                dth1 = Sp(3,js-1);
                s1_d = - zbar_d*dth1 + sr_d + cv_d;
                bal_d = sumr_d - sumq_d - sumae_d + (s0_d - s1_d);
                if (idebugoutput >= 1) { // dgt 8/17/05 debugout
                    lunpFile << "water balance: bal= " << bal_d << " sumr= " << sumr_d << '\n';
                }
                if ( fabs(bal_d) > (sumr_d + sumq_d)*1.0e-3 ) {
                    // can't write to lunco1, so put it in toperror.txt
                    //			write(21,*)"water balance: bal=",bal_d," sumr=",sumr_d,
                    //    +					" sumq=",sumq_d," it=",istep
                }
                if (idebugoutput >= 1) {  // DGT 8/17/05 debugout
                    lunpFile << '\n';
                    lunpFile << " Rainfall    =" << fixed << setw(9) << setprecision(1) << sumr_d/units;
                    lunpFile << " Total Evap  =" << fixed << setw(9) << setprecision(1) << sumae_d/units << '\n';
                    lunpFile << " Model Runoff=" << fixed << setw(9) << setprecision(1) << sumq_d/units;
                    lunpFile << " Total PET   =" << fixed << setw(9) << setprecision(1) << sumpe_d/units << '\n';
                    lunpFile << " Incr.Storage=" << fixed << setw(9) << setprecision(1) << (s1_d - s0_d)/units << '\n';
                }
            } // finish

            vol_irrig_demand(js-1) = depth_irrig_dem;  // *(Sp(1,JS)/1e6) !m3/timestep = m/timestep * m^2

            if ( modwrt && istep == m ) {
                if (idebugoutput >= 1) {  // DGT 8/17/05 debugout
                    lunpFile << " The parameters are\n";
                    for (i = 1; i <= Nsp; i++) {
                        lunpFile << " " << scientific << setw(15) << setprecision(6) << Sp(i-1,js-1) << '\n';
                    }
                    //lunpFile << '\n';
                    lunpFile << " Initial conditions\n";
                    for (i = 1; i <= Nsi; i++) {
                        lunpFile << " " << scientific << setw(15) << setprecision(6) << Si[i-1][js-1] << '\n';
                    }
                    lunpFile << '\n';
                }
            }
            //     SET UP THE START TIME OF THE OUTPUT TIME SERIES OF EACH GUTTER
        } // subbasin_loop

        // THAT IS THE END OF THE SUBCATCHMENTS

        // now we know runoff and soil moisture, so we can do a water managment step

        if (istep > 0) {
            for (i = 0; i < Nsub; i++) {
                precip_for_watermgmt[i] = bRain(i,istep-1); //mm/timestep
            }
        }
        // DGT 10/21/12  Additional writes for Christina
        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[0], "results/Potential_evapotranspiration_mm.txt", istep, potentialevap, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[1], "results/TemperatureAve_C.txt", istep, tempave, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[2], "results/Surface_runoff_cms.txt", istep, surfro, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[3], "results/Canopy_storage_mm.txt", istep, canstore, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[4], "results/Soil_storage_mm.txt", istep, soilstore, Nsub, scalefactor);

        scalefactor = 1000.0;
        ArrayXd zbar1 = zbar.col(0);
        Write_OutputLine_Eigen(oFile[5], "results/Depth_to_Water_mm.txt", istep, zbar1, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[6], "results/Tile_drainage_cms.txt", istep, tiled, Nsub, scalefactor);

        scalefactor = 1.0;
        Write_OutputLine_Eigen(oFile[7], "results/Ditch_drainage_cms.txt", istep, ditchd, Nsub, scalefactor);

        // Added ArtDrainage to arguments for output
        watermgmt(sDate, sHour, istep, m, totalrunoff, baseflow, ArtDrainage, vol_irrig_demand, maxSlp,
                  evap_for_watermgmt, precip_for_watermgmt, volume_irrig_sup, groundwater_to_take);

        if (istep > 0)
            updatetime(iyear, month, iday, hour1, timestep);

    }	// time_loop
    // End of time loop  for reaches ***************************************************
    // Added ArtDrainage to arguments for output
    //final call to write/close output files
    watermgmt(sDate,sHour, -1, m, totalrunoff, baseflow, ArtDrainage, vol_irrig_demand, maxSlp, //inputs, ignored this time
              evap_for_watermgmt, precip_for_watermgmt, volume_irrig_sup, groundwater_to_take);
    if (idebugoutput >= 1) {
        lunpFile << "Total snow mass balance error = " << errmbal << '\n';	// 'topinfo_v7.txt'
    }

    //  Close files on additional writes for Christina
    istep = -1;
    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[0],"results/Potential_evapotranspiration_mm.txt", istep, potentialevap, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[1],"results/TemperatureAve_C.txt", istep, tempave, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[2],"results/Surface_runoff_cms.txt", istep, surfro, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[3],"results/Canopy_storage_mm.txt", istep, canstore, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[4],"results/Soil_storage_mm.txt", istep, soilstore, Nsub, scalefactor);

    scalefactor = 1000.0;
    ArrayXd zbar1 = zbar.col(0);
    Write_OutputLine_Eigen(oFile[5],"results/Depth_to_Water_mm.txt", istep, zbar1, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[6], "results/Tile_drainage_cms.txt", istep, tiled, Nsub, scalefactor);

    scalefactor = 1.0;
    Write_OutputLine_Eigen(oFile[7], "results/Ditch_drainage_cms.txt", istep, ditchd, Nsub, scalefactor);

    if (idebugoutput >= 1)
        lunmodFile.close();   // debugout: close model writes in case routing fails

    for (js = 0; js < maxSlp; ++js) {
        for (kk = 0; kk < nreg; ++kk) {
            delete [] dr[js][kk];
            delete [] qinst[js][kk];
        }
        delete [] dr[js];
        delete [] qinst[js];
    }
    delete [] dr;
    delete [] qinst;

    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j < Nsub; ++j) {
            delete [] zbar_in[i][j];
        }
        delete [] zbar_in[i];
    }
    delete [] zbar_in;

    delete [] groundwater_to_take;
    delete [] evap_for_watermgmt;
    delete [] precip_for_watermgmt;

#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving calcts(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

    return 0;
}


// *****************************************************************
//     SUBROUTINE SNOW
// *****************************************************************
int snow(ofstream &snowOutFile, const ArrayXd &temper, const double elevTg, const double elevsb, const double rlapse, double &bRain,
         const double ddf, double &snowst, const long int dt, const int Nsub, const int m, const int js, const int it,
         const int maxSlp, const int maxInt)
{
    //SAVE
    double lastsnowst;
    double t, temp, newsnow, melt, snowbal;

#if TRACE
    static int ncalls = 0;
    string save_caller = caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> snow(" << ncalls << ")" << std::endl;
    }
    caller = "snow";
#endif

    if (it == 1) {
        //	if (js == 1) open(unit=78,file="snowout.txt")
        snowOutFile << dec << js << " " << Nsub << " " << m;
        snowOutFile << "it, precip(mm), snow(mm), rain(mm), temper(deg), snowstore(mm), melt(mm), snowbal(mm)\n";
    }

    lastsnowst = snowst;
    t = temper(it-1) - (elevsb - elevTg)*rlapse;
    temp = bRain;
    if (t < 0) {
        newsnow = bRain; //bRain in mm/interval, newsnow in mm
        snowst += newsnow; //snowst and newsnow both in mm
        bRain = 0.0;
        melt = 0.0;
    }
    else {
        newsnow = 0.0;
        melt = min(snowst, ddf*t*(double)dt/86400.0); //ddf in mm/deg/day, t in deg, dt in seconds, melt in mm
        snowst -= melt;	//melt in mm, snowst in mm
        bRain += melt; //melt in mm, bRain in mm/interval
    }
    snowbal = snowst - lastsnowst + melt - newsnow;
    snowOutFile << dec << setw(8) << it;	// unit 78
    snowOutFile << fixed << setw(12) << temp;
    snowOutFile << fixed << setw(12) << newsnow;
    snowOutFile << fixed << setw(12) << bRain;
    snowOutFile << fixed << setw(12) << t;
    snowOutFile << fixed << setw(12) << snowst;
    snowOutFile << fixed << setw(12) << melt;
    snowOutFile << fixed << setw(12) << snowbal;
#if TRACE
    caller = save_caller;
    if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving snow(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

    return 0;
}

