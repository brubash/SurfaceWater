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

#ifndef SNOW_HH
#define SNOW_HH

extern std::ofstream snowcontrol3_File;	// unit 10
double svpi(const double t);
double svpw(const double t);
int qlif(double &qliff, const double ta, const double rh, const double tk, const double sbc, const double cf);
double ablat(const double DimBaw, const int ndfc, double **dfc);
int snowLSub(int &year, int &month, int &day, double &hour, const int dt,
	const int nStep, Eigen::ArrayXXd &inpt, Eigen::ArrayXd &sitev, Eigen::ArrayXd &statev,
	const Eigen::ArrayXd &param, Eigen::ArrayXi &iflag, const std::vector<double> dtbar, const int nstepday,
    double &cump, double &cume, double &cummr, Eigen::ArrayXd &outv,
    Eigen::ArrayXd &tsbackup, Eigen::ArrayXd &tavebackup,
    const int ndepletionpoints, double **dfc, const int modelelement, const int jj);
int snowueb2(const double dt, const int nt, const Eigen::ArrayXXd &input, const Eigen::ArrayXd &sitev, Eigen::ArrayXd &statev,
	Eigen::ArrayXd &tsprevday, Eigen::ArrayXd &taveprevday, const int nstepday, const Eigen::ArrayXd &param, Eigen::Array<int,Eigen::Dynamic,1> &iflag,
	double &cump, double &cume, double &cummr, Eigen::ArrayXd &outv, const Eigen::ArrayXd &mtime, const int modelelement, const int jj);
double daily_ave(const Eigen::ArrayXd &backup, const int n, const double a);
double partsnow(const double p, const double ta, const double tr, const double ts);
double tavg(const double ub, const double w, const double rhow, const double cs,
	const double to, const double rhog, const double de, const double cg, const double hf);
double albedo(const double tausn, const double coszen, const double d,
	const double aep, const double abg, const double avo, const double airo);
int agesn(double &tausn, const double dt, const double ps,
	const double tsurf, const double tk, const double dnews);
int predicorr(const double dt, double &ub, double &w, const double a,
	const double ta, const double pRain, const double ps, const double ws, const double rh,
	const double qsi, const double  qli, const double iradfl, const double rkn,
	const double qnetob, const double rid,
	const Eigen::ArrayXd &param, const Eigen::ArrayXd &sitev, const double iTsMethod,
	const Eigen::ArrayXd &mtime,
	double &qh, double &qe, double &e,  double &mr, double &qm, double &q, double &fm,
	double &tsurf, double &tave, double &qnet, double &refDepth, double &totalRefDepth ,
	double &smelt);
int qfm(const double ub, const double w, const double a, const double ta, const double pRain,
	const double ps, const double ws, const double rh, const double qsi, const double qli,
	const double rkn, const double iradfl, const double qnetob,
	const double rid, const Eigen::ArrayXd &param, const Eigen::ArrayXd &sitev, const int iTsMethod, const Eigen::ArrayXd &mtime,
	double &fm, double &q, double &qm, double &mr, double &qe,
	double &e, double &tsurf, double &tave, double &qh, double &qnet,
	const double dt, double &refDepth, const double totalRefDepth, double &smelt);
int prehelp(const double w1, const double w, const double dt, double &fm,
	const double fm1, const double fac, const double ps, const double pRain,
	double &e, const double rhow, const double hf, double &q, double &qm,
	double &mr, double &qe, const double hneu);
double qpf(const double pr, const double ta, const double to, const double ps,
	const double rhow, const double hf, const double cw, const double cs);
double rkinst(const double rkn, const double ws, const double ta,
	const double ts, const double z, const double g, const double fstab);
double QcEst(const double ts,  const double rkn, const double ws, const double tak,
	const double g, const double qp,
	const double densa, const double cp, const double hneu,
	 const double pr, const double ea, const double tk, const double dens,
	 const double cs, const double rs, const double tavek, const double qsn,
	 const double qli, const double fc, const double sbc,
	 const double qnetob, const int iradfl, const Eigen::ArrayXd &mtime,
	 const Eigen::ArrayXd &param, const int iTsMethod, const double w, const double dt);
int Grad(const double qc1, const double qc2, const double t1, const double t2, double &a, double &b);
double refDep(const double flans, const double a, const double b, const double hf,
	const double rhom, const double dt, const double x1 );
double LanE(const double lans, const double lang, const double zs, const double rho,
	const double rhog, const double cs, const double cg, const double r, double &ze, const double w1day);
double fmelt(const double ub, const double rhow, const double w, const double hf,
	const double lc, const double rid, const double ks, const double pRain);
int turbflux(const double pr, const double ra, const double ta, const double tk,
	const double ts, const double z, const double g, const double cp,
	const double rkn, const double ws, double &ea, const double rhow,
	const double hneu, double &qh, double &qe, double &e, const double fstab);
double srftmp(const double qsi, const double a, const double qli, const double qpin,
		const double ea, const double ta, const double tave, const double tk, const double pr,
		const double ra, const double cp,const double rho, const double rkn, const double hneu,
		const double es, const double sbc, const double cs, const double rs, const double w,
		const double qnetob, const double iradfl, const double ws, const double z, const double g,
		const double fc, const double fstab, const Eigen::ArrayXd &mtime, const Eigen::ArrayXd &param,
        const int iTsMethod, const double dt, const double ub, const double refDepth, double &smelt);
double surfeb(const double ts, const double rkn, const double ws, const double tak, const double g,
	const double qp, const double densa, const double cp, const double hneu, const double pr, const double ea,
	const double tk, const double dens, const double cs, const double rs, const double tavek, const double qsn,
	const double qli, const double fc, const double sbc, const double qnetob, const int iradfl,
	const Eigen::ArrayXd &mtime, const Eigen::ArrayXd &param, const int iTsMethod, const double w, const double dt,
	const double ub, const double refDepth);

#endif
