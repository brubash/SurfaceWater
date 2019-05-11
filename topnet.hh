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

    develop_1 branch added May 11, 2019
*/

#ifndef TOPNET_HH
#define TOPNET_HH

#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <Eigen/Dense>
#include "Nlfit.hh"
#include "Tdims.hh"
#include "types.hh"
#include <vector>
#include <valarray>
#include <array>

//#define TRACE 10
#define MAX_TRACE 5
#define ZBAR_OUT
#define ZBAR_IN
//#define WRIA1
//#define LNWB

#if TRACE
    extern std::ofstream traceFile;
    extern std::string caller;
    const int maxCalls = 10000;
#endif
#ifdef ZBAR_OUT
    extern std::ofstream zbarFile;
#endif

//extern int *ifound;
//extern std::ofstream debugFile;
extern std::ofstream topErrorFile;
extern std::vector<int> ifound;
extern std::vector<int> find_counts;
extern std::vector<int> done;
extern int nfound, nfound_last, global_found;
void Find(Eigen::Array<DrainageType,Eigen::Dynamic,1> &iarray1, const int ival1,const int num);

namespace data_array {
	extern Eigen::ArrayXXd real_array;
	extern Eigen::ArrayXXi integer_array;
	extern double **dble_array;
	extern bool real_array_allocated;
	extern bool integer_array_allocated;
	extern bool dble_array_allocated;
	const int real_data = 1;
	const int integer_data = 2;
	const int dble_data = 3;
}

namespace TimeVaryingOutput {
	extern Eigen::ArrayXXd FlowInLinks_cms;			//Timestep,Location
	extern double *FlowAtStreamNodes_cms;		//Location
	extern Eigen::ArrayXXd ReservoirStorage_m3;		//Timestep,Location
	extern Eigen::ArrayXXi DateTime_yyyymmdd_hhmmss;		//Timestep, yyyymmdd/hhmmss
}
namespace input_structures {
	extern RunControlType RunControl;
	extern Eigen::Array<DrainageType,Eigen::Dynamic,1> Drainage;
	extern StreamNodeType *StreamNode;
	extern MeasuredFlowInfoType *MeasuredFlowInfo;
	extern MeasuredFlowDataType *MeasuredFlowData;
	extern RunoffType *Runoff;
	extern RunoffType *Baseflow;
	extern RunoffType *ArtDrainage;
	extern ReservoirType *Reservoir;
	extern UserType *User;
	extern SourceType *Source;
	extern RightsType *Rights;
	extern SourceMixingType *SourceMixing;
	extern SeasonsDefnType *SeasonsDefn;
	extern ReturnFlowType *ReturnFlow;
	//	DemandCoefficientsType *DemandCoefficients;
	extern MonthlyDemandType *MonthlyDemand;
}

namespace other_structures {
	extern NodeType              *Node;
	extern Eigen::Array<LinkType,Eigen::Dynamic,1>      Link;
	extern UserSourceTableType   *UserSourceTable;
	extern StaticOutputTableType StaticOutput;
	extern NodeType              *NodeSave;
	extern Eigen::Array<LinkType,Eigen::Dynamic,1>      LinkSave;
	extern int                   *UserSourceOrder;
	extern int                   *WWTP_list;
}

namespace constant_definitions {

	//Define UseType Codes
	//UserType Descriptions
	const int SoilMoistureIrrigationUseCode = 1;
	const int FixedDemandIrrigationUseCode = 2;
	const int DownstreamReservoirReleaseUseCode = 3;
	const int PWSUseCode = 4;
	const int NonPWSMandIUseCode = 5;
	const int DairyUseCode = 6;
	const int RanchUseCode = 7;
	const int PoultryUseCode = 8;
	const int ParkGolfCemeteryUseCode = 9;
	const int InstreamFlowUseCode = 10;
	const int DiversionUseCode = 11;
	//  The above are user specified.  The below are used internally for automatically created users to handle reservoirs
	const int ReservoirFillUseCode = 12;
	const int InStreamReservoirReleaseUseCode = 13;
	const int OffStreamReservoirReleaseUseCode = 14;

	//Define SourceType Codes
	const int StreamSourceCode = 1;
	const int GroundwaterSourceCode = 2;
	const int ReservoirSourceCode = 3;

	//Define ReservoirType Codes
	const int InStreamReservoirCode = 1;
	const int OffStreamReservoirCode = 2;

	//Define Internal/External Codes
	const int ExternalCode = 1;
	const int InternalCode = 2;

	//Define NodeType Codes
	const int StreamNodeCode       = 1; //Drainage Outlets
	const int GroundwaterNodeCode  = 2;
	const int ReservoirNodeCode    = 3;
	const int DrainageNodeCode     = 4;
	const int UserNodeCode         = 5;
	const int SinkNodeCode         = 6;
	const int MeasuredFlowNodeCode = 7;

	//Define LinkType Nodes
	const int SurfaceRunoffLinkCode = 1;
	const int SubsurfaceRunoffLinkCode = 2;
	const int UnallocatedLinkCode = 3;    // This is used to designate the link at each node that absorbs
	// the extra flow when flow at nodes is balanced.  This is always the link from a stream node to the next downstream node.
	const int UserAbstractionLinkCode = 5;
	const int ReturnFlowLinkCode = 6;
	const int SinkLinkCode = 7;
	const int MeasuredFlowLinkCode = 8;

	//Define AllocationMode Types
	const int WaterRightsAllocationCode = 1;
	const int DemandAllocationCode = 2;
	const int NoAllocationCode = 0;

	//Define StreamNode Types
	const int IsDrainageOutlet = 1;
	const int IsNotDrainageOutlet = 0;

	//Define Units Types
	const int FractionUnits = 1;
	const int VolumeUnits = 2;
	const int FracMinDemandUnits = 3;

} //end module constant_definitions

#define   Ntsg(i,j)   Ntsg[i-1][j-1]
#define    atb(j,i)    atb[j-1][i-1]
#define    pka(j,i)    pka[j-1][i-1]
#define    cl2(j,i)    cl2[j-1][i-1]
#define    pd2(j,i)    pd2[j-1][i-1]
#define     si(j,i)     si[j-1][i-1]
#define     Rp(j,i)     Rp[j-1][i-1]

extern int NumDrainage;
extern int NumStreamNode;
extern int NumMeasuredFlowInfo;
extern int NumMeasuredFlowData;
extern int NumReservoir;
extern int NumUser;
extern int NumSource;
extern int NumRights;
extern int NumSourceMixing;
extern int NumSeasonsDefn;
extern int NumReturnFlow;
extern int NumMonthlyDemand;
extern int NumRunoff;
extern int NumBaseflow;
extern int NumWWTP;
extern std::ofstream oFile[];

namespace mddata {
	extern int maxRno;
}
namespace lakes1 {
	extern int nlakes;
}
namespace ts_save {
	extern double Ts_old;
	extern double Tave_old;
	extern double Ts_Ave;
	extern double Tave_ave;
}
namespace rain {
	extern double rain_factor;
}
namespace ccompr2 {
	extern double t0;
}
namespace model_master {
	extern int iDoIt;
}
namespace model1 {
	extern int nGauge;
	extern int Ns;
	extern int Nrch;
}
namespace model2 {
	extern int sDate;
	extern int sHour;
	extern long int interval;
	extern int m;
	extern int mi;
	extern int mps;
	extern int mpe;
	extern int stim;
}
namespace model4 {
	extern Eigen::ArrayXXi pmap;
	extern int Nout;
	extern int nBout;
}
namespace model5 {
	extern int Neq1;
}
namespace constr {
	extern double minSp[];
	extern double maxSp[];
	extern double minSi[];
	extern double maxSi[];
	extern double minRp[];
	extern double maxRp[];
	extern bool limitC;
}

// from maxvariables.inc
extern int maxInt, maxGauge, maxSlp, maxTGauge;
extern int maxResponse, maxA, maxC, maxChn;
extern int maxSites, maxRchAreas;
extern int max_lakes, max_lheads;
extern int max_of_shifts;

extern int Nrchsav;

extern std::ifstream topinpFile, modelspcFile, rchareasFile, rainFile;
extern std::ofstream lunmodFile;
extern std::ofstream luntopFile;
extern std::ofstream lundatFile;
extern std::ofstream lunpFile;

// watermgmt.cpp
extern int *DrainageOrder;
//extern double *DrainageOutFlow;

//extern double rain_factor;	 common to hydatasn and mddata files
//extern int maxRno;			 common to calv46sn and mddata files

int AllocateWaterToUsers(const int Timestep, const int NumNode, const int NumLink, const int NumUser, const int NumReservoir,
	const int NumSource, const int NumRights, const int NumSourceMixing, int *DrainageOrder, const int NumDrainage,
	const int NumReturnFlow, const int NumUserSource, Eigen::ArrayXd &volume_irrig_sup,
	double *groundwater_to_take, Eigen::ArrayXd &DrainageOutFlow);
int Append_To_Output_Tables(const int Timestep, const int dt,
	const int NumStreamNode, const int yyyymmdd, const int hhmmss,
	const int NumLink, const int NumNode, const int Nsteps, const int NumUser);
int AssignDrainageFlows(const int Timestep, const int NumDrainage, const int NumNode, const int NumLink,
	int *DrainageOrder, const int NumRunoff, const int NumBaseflow, Eigen::ArrayXd &DrainageOutFlow);
int AssignPriorityOrder(const int NumUser, const int NumSource, const int NumRights, const int NumReservoir,
	const int AllocationMode, int *DrainageOrder, const int NumDrainage,
	int &NumUserSource, int &NumUserSourceReturn, const int NumReturnFlow);
int atf(double &atff, const double trange, const int month, const std::vector<double> dtbar, const double a, const double c);

int BalanceFlowsAtReservoirs(const int NumNode, const int NumLink, const int NumUser,
	const int NumReservoir, Eigen::ArrayXd &ReservoirNetStorage);
int BalanceFlowsAtStreamNodes(const int NumNode, const int NumLink, int *DrainageOrder,
	const int NumDrainage, Eigen::ArrayXd &DrainageOutFlow);
int bcparm(std::vector<double> &dtBar, double &bca, double &bcc, const std::string bcFile);
int BuildDrainageOrder(const int NumDrainage, int *DrainageOrder);
int BuildLinkStructure(const int NumDrainage, const int NumUser, const int NumSource,
	const int NumReturnFlow, const int NumReservoir, const int NumMeasuredFlowInfo, const int NumNode, int &NumLink);
int BuildNodeStructure(const int NumDrainage, int &NumUser, const int NumReservoir, int &NumSource,
	int &NumReturnFlow, const int NumMeasuredFlowInfo, int &NumNode);
int calcts( double **Si,            const Eigen::ArrayXXd &Sp,           double **Rp,           Eigen::ArrayXXi &linkR,        const int *ll,
            const int Nsub,                                const int *Nka,        const double *tl,   double **atb,
            double **pka,           const int *nd,         double **cl,           double **pd,        const double units,
            const int ipsub,        const int ipatb,       const bool reinit,     const bool modwrt,  const int stim,
            Eigen::ArrayXXd &bRain,        const long int interval, const int m,        const int mi,
            const int mps,          const int mpe,         bool &ok,              const double *xlat, const double *xlong,
            const double stdlon,    const double *elevtg,  double **bdtBar,       const int sDate,    int &sHour,
            const Eigen::ArrayXd &temper,   const double *dewp,    const double *tRange,  const int Neq,      const int Nout,
            const int nBout,        const int *iBout,      const double *wind2m,  double **bTmin,     double **bTmax,
            double **bTdew,         const double *bXlat,   const double *bXlon,   int *ntdh,    const int ndump,
            const int maxInt,       const int maxSlp,      const int maxA,        const int maxC,     const int maxChn,
            const int idebugoutput, const int idebugbasin, const int idebugcase);
int CalculateDemand(const int ThisMonth, const int doy, const Eigen::ArrayXd &vol_irrig_demand, const int NumUser, const int NumMonthlyDemand);
int chk_sites(int &lrg, const int sites, const int j, int &iMatch);
int cliParam(double *xlat, double *xlong, double &stdlon, double *elevtg, double **dtBar, int &ns_temper, int *temper_id);
int CreateLink(int &NumLink, const std::string Title, const int LinkCode, const int IntExtCode,
	const int USNode, const int DSNode, const double Flow, const int ReturnFlowID);
int etall(const double xlat, const double xlong, const double stdlon, const double elevtg,
	const std::vector<double> dtbar, double &evap, const double temp, const double dewp, const double trange,
	const double elevsb, const double albedo, const double rlapse, const int sdate, const int shour,
	const int dtsec, const int m, const int istep, int &iyear, int &month, int &iday,
	int &ihr, int &imm, int &isec, double &hour1, const double tmin, const double tmax,
	const double wind2m, const int method);
int hyri(const int year, const int month, const int day, const double hour, const double dt,
	const double slope, const double azi, const double lat, double &hri, double &coszen);
int ImposeMeasuredFlows(const int Timestep, const int NumNode, const int NumLink, const int NumRunoff,
	const int NumBaseflow, const int NumDrainage, const int NumMeasuredFlowInfo, const int NumMeasuredFlowData,
	int *DrainageOrder, Eigen::ArrayXd &DrainageOutFlow);
int inputT(int initT[], int iend[], int &Neq, Eigen::ArrayXXd &qact, double actime[],
	std::string &modelid, int &Npar, const int nrx, const int iex);
int irrigation(const double Sr, const Eigen::ArrayXXd &Sp, const int js, const long int interval, const int MAXSLP, double &Dapp);
int read_bdryflow(const std::string fileName, int &expected_numcols, const int ncommentlines, int &nrows);

int Model(const int it, const int iflag, const int iopt, const char prt, int &Neq, int &Npar, const int npx, const int iEx,
	const int nfor, double qfit[], const double par[], const int mfit[], const int ifit, const int ibeale);
int mdData(int &Ngauge, int &Ns, int &Nrch, int *Nka, double *tl, double **atb, double **pka, int *Nd, double **cl2, double **pd2, double &units,
	int *ll, int *Ntr, int *Nts, Eigen::ArrayXXi &linkS, Eigen::ArrayXXi &linkR, double **si, Eigen::ArrayXXd &Sp, double **Rp,  const int iret,
	Eigen::ArrayXXi &pMap, int &Npar, Eigen::ArrayXXi &lrg, Eigen::ArrayXXd &wrg,
	const int iex, int *llout, int &Neq, int &Nout, int &nBout, int *iBout, int &nRchSav, int *qMap, double *rel, int &relFlag,
	double *minSp, double *maxSp,  double *minSi, double *maxSi, double *minRp, double *maxRp, bool &limitC, Eigen::ArrayXXi &ClinkR, Eigen::ArrayXi &kllout,
	int *ishift0, const int maxInt, const int maxGuage, const int maxSlp, const int maxResponse, const int maxA, const int maxC,
	const int maxChn, const int maxRchAreas, const int maxSites, Eigen::ArrayXXd &bp,
	double *bXlat, double *bXlon, Eigen::ArrayXXd &wrg1);
int hyData(int &sDate, int &sHour, long &interval, int &m, int &mi, int &mps, int &mpe, int &Ngauge, int &Neq,
	Eigen::ArrayXXd &bRain, double **flow, int &iret, Eigen::ArrayXd &temper, double *dewp, double *trange, double **dtBar, const int Ns, Eigen::ArrayXXd &wrg,
	Eigen::ArrayXXi &lrg, const double *elevtg, double **bTmax, double **bTmin, double **bTdew, double **bdtBar, const Eigen::ArrayXXd &Sp,
	const int maxGauge, const int maxInt, const int maxSites, const int maxResponse, const int maxTGauge, double *wind2m, Eigen::ArrayXXd &wrg1,
	int &idebugoutput, int &idebugbasin, int &idebugcase);

int PropagateWaterViaUser(const int i, const int j, const double Qtry, const int NumNode,
	const int NumLink, const int NumUser, const int NumReservoir, const int NumSource,
	int *DrainageOrder, const int NumDrainage, const int NumReturnFlow, int &iFeasible, double &Capacity);
int rainfill_doit(Eigen::ArrayXd &train, const int it, const int Ngauge, Eigen::ArrayXXd &rcoeff, Eigen::ArrayXXi &fsite, const int *nfill,
	bool &i_reset_flag, int &it_save_old, const int maxGauge);
int rainfill_read(const std::string rainfill_file, const int Ngauge, const int *rsite,
	Eigen::ArrayXXd &rcoeff, Eigen::ArrayXXi &fsite, int *nfill, int &fillflag, const int maxGauge);
int isrch(const int mSite, const int *rSite, const int nGauge);
bool miss(const double train, const int nfill);

int iposn(const int Nrch, const Eigen::ArrayXXi &linkR, const int ii);
int set_cor_data(const Eigen::ArrayXXi &linkR, const int Neq, const int Nout, const int Nrch, const int Ns,
	const Eigen::Array<int,Eigen::Dynamic,1> &kllout, double **flow, const Eigen::ArrayXXd &Spd, const double *bArea, const int *ll,
	const int *llOut, double *dFlow, double *suma, int *knt, int *iReach, int *lOut, const int maxInt,
	const int maxSlp, const int maxChn, const int maxResponse);
int read_lakes(Eigen::Array<int,Eigen::Dynamic,1> &lake_reach, int *lzero, double *lake_areas, int *lake_beach_slps, int *lk_line, int *num_rat_vals,
	int **lheads, int **loflows, const int max_lakes, const int max_lheads);
int read_lakes_levels(const Eigen::Array<int,Eigen::Dynamic,1> &lake_reach, double *ini_levels, const int max_lakes);
int setup_zbar0s(const int Ns, double **Si, const Eigen::ArrayXXd &Sp, const double *tl, const long int dt,
	const int *lOut, const int *iReach, const double *dFlow, const double *sumA, const int *knt, const int maxResponse, const int maxSlp);
int solve_for_zbar0(double &zbar0, const double *area, const double *f, const double *k,
	const double *Lambda, const double q, const long int dt, const int ns);

int updatetime(int &year, int &month, int &day, double &hour, const double dt);

int datevec(const int yyyymmdd, const int hhmmss, int &yyyy, int &mm, int &dd, int &hh, int &mi, int &ss);
int undatevec(const int yyyy, const int mm, const int dd, const int hh, const int mi, const int ss, int &yyyymmdd, int &hhmmss);

int read_inputs(const std::string dirname, const int dt, const int StartDateTopnet, int &NumDrainage,
	int &NumStreamNode, int &NumMeasuredFlowInfo, int &NumMeasuredFlowData,
	int &NumReservoir, int &NumUser, int &NumSource, int &NumRights, int &NumSourceMixing,
	int &NumSeasonsDefn, int &NumReturnFlow, int &NumMonthlyDemand, const int NumRunoff,
	int &NumBaseflow, int &NumWWTP);
int read_struct_from_text(const std::string fileName, int &expected_numcols, const int ncommentlines, int &nrows);

int Initialise_Output_Tables(const int NumDrainage, const int NumNode, const int NumStreamNode, const int NumLink, const int NumUser,
	const int NumTimesteps, const int NumReservoir, const int NumUserSourceReturn, const int NumReturnFlow);

int snow(std::ofstream &snowOutFile, const Eigen::ArrayXd &temper, const double elevTg, const double elevsb, const double rlapse, double &bRain,
	const double ddf, double &snowst, const long int dt, const int Nsub, const int m, const int js, const int it,
	const int maxSlp, const int maxInt);
int snowueb(const Eigen::ArrayXd &snowsitev, Eigen::ArrayXd &snowstatev1, const Eigen::ArrayXd &snowparam,
    const int ndepletionpoints, double **dfc, const Eigen::Array<int,Eigen::Dynamic,1> &snowcontrol, const std::vector<double> dtbar,
	const Eigen::ArrayXd &snowforcing, Eigen::ArrayXd &snowsurfacetemp1, Eigen::ArrayXd &snowaveragetemp1,
    const double timestep, const int nstepday, double &surfacewaterinput, double &snowevaporation,  //outputs (both in m/h)
    double &areafractionsnow, const int modelelement);
double svp(const double t);
int td8micsec(const int jdatem, int &jhour, long int &isec);
int td81micdh(int &idate, int &ihour, const long int jsec);
int topmod(double **si, const Eigen::ArrayXXd &Sp, const int isub, const int *Nka, const double Lambda,
	double **atb, double **pka, const int *nd, double **cl, double **pd, const double units,
	Eigen::Array<int,Eigen::Dynamic,1> &irr, const bool modwrt, const int ipsub, const int ipatb, const int stim,
	const double r, const double pet, const long int interval, const double art_drainage, const double rate_irrig,
	const int imonth, const int ndata, const int mps, const int mpe, double &qinst_out, double &dr_out,
	const int ndump, int *ntdh, const int istep, const int maxC, double &zbm, const int maxA,
	const int maxSlp, const int maxInt, double &sumr, double &sumq, double &sumae,
	double &s0, double &q0, double &sr, double &cv, double &aciem,
	double &acsem, double &sumpe, double &sumie, double &sumqb, double &sumce,
	double &sumsle, double &sumr1, double &qb, Eigen::ArrayXd &qinst, Eigen::ArrayXd &dr, double &sumqv,
	double &sumse, double &zbar, const double zbar_new, double **tdh, double &zr, double &ak0fzrdt, double &logoqm,
	double &qvmin, double &dth, double &sumad, double &evap_mm, double &qlat_mm,
	const int ipflag, Eigen::ArrayXd &rirr, const int js);
int watermgmt(const int StartDateTopnet, int &StartHourTopnet, const int Timestep, const int NSteps,
	Eigen::ArrayXd &RunoffTopnet, Eigen::ArrayXd &BaseflowTopnet, const Eigen::ArrayXd &ArtDrainageTopnet,
	const Eigen::ArrayXd &vol_irrig_demand, const int maxSlp, const double *evaporation, const double *precipitation,
	Eigen::ArrayXd &volume_irrig_sup, double *groundwater_to_take);
int Write_OutputLine(std::ofstream &oFile, const std::string filenm, const int timestep, const double *Rvariable,
	const int NumDrainage, const double scalefactor);
int Write_OutputLine_Eigen(std::ofstream &oFile, const std::string fileName, const int timestep, const Eigen::ArrayXd &Rvariable,
	const int NumDrainage, const double scalefactor);
int Write_OutputLocalContributions(std::ofstream &o1File, std::ofstream &o2File, const int NumStreamNode, const int NumDrainage,
		const Eigen::ArrayXd &BaseflowTopnet, const Eigen::ArrayXd &RunoffTopnet, const int timestep, const double scalefactor);
int Write_Static_Output_Tables(const std::string dirname, const int NumUserSourceReturn);
int Write_TimeVaryingOutput_Tables(const std::string dirname, const int NumUser, const int NumTimesteps,
	const int NumNode, const int NumLink, const int NumReturnFlow, const int NumWWTP);

#endif

// Fortran indexing notes:
//   LINKR(4,MAXCHN) -- array containing reach network information ,
//     i.e. LINKR(1,I) = arbitrary numbering (must be greater than NS).
//          LINKR(2,I) = slope or upstream reach feeding into this reach.
//          LINKR(3,I) = second slope or upstream reach feeding into this
//          LINKR(4,I) = number of slopes or upstream reaches feeding into

//  LL(MAXCHN) -- pointer giving the order for processing all reaches.

//  Sort out llout array so that it is an index to ll array
// Array LL defines the processing order of reaches. Array LLOUT contains
// the position within LL of reaches with measured values or for which
// output is required.

// Set up the mapping as follows: use an array called qmap. Define
// qmap(i) to be the column in runoff.dat corresponding to i'th reach
// in the list read in from modelspc.dat, e.g. through the reach-site no.
// linkage data for the 2nd reach in the list might be in the 4th column,
// in which case qmap(2)=4

//  LRG(MAXSLP,MAXGAUGE) -- Extended basin-raingauge link data
//    i.e. LRG(1,I) = raingauge numbers to be used for sub-basin I.
// Adjusted for C++ indexing when read from input file in mdData. (lrg[i][k]--;)
