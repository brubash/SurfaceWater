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

using namespace std;

int get_fx_and_dfx(const double zbar0, const double *area, const double *f, const double *k,
	const double *Lambda, const double q, const long int ns, double &fx, double &dfx, const long int dt);


int solve_for_zbar0(double &zbar0, const double *area, const double *f, const double *k,
	const double *Lambda, const double q, const long int dt, const int ns)
{
	// this program uses the newton raphson method to solve for the roots of eqns
	// The user supplied equation is in the function subprogram called 'func'
	// The user supplied derivative of the equation is located in 'dfunc'.

	int iter, maxit;
	double dif, tol, x, xold, fx, dfx;
	// *****************Variable Dictionary***********************************
	// x = the current best guess for the location of the root.
	// xold = the location of the root from the previous iteration.
	// func = the user supplied function of x
	// dfunc = the derivative of the function (also user supplied)
	// dif = the % difference between x and xold (used for convergence criteria)
	// tol = tolerance for convergence
	// maxit = maximum number of iterations allowed
	// ***********************************************************************
	//  initialize variables
	iter = 0;
	dif = 1.0;
	tol = 0.001;
	maxit = 50;

	//  put in initial guess for root.
	x = zbar0;

	//  start loop to perform newton-raphson method
L15: ;
	xold = x;
	get_fx_and_dfx(x, area, f, k, Lambda, q, ns, fx, dfx, dt);
	x = xold - fx/dfx;

	//  add 1 to # of iterations.  If maximum number of iterations (maxit)
	//  is exceeded then stop the program.
	iter++;
	if(iter > maxit)
		cerr << "Exceeded maximum number of iterations\n";

	//  Check the difference between current guess for the root and the
	//  previous (x vs xold).  If % difference is less than the tolerance,
	//  convergence is reached.  Otherwise, go back up and repeat the process.
	dif = fabs((x - xold)/x);
	if(dif > tol)
		goto L15;


	zbar0 = x;

	return 0;
}

// *****************Subroutine***********************************
//  user suplied function
int get_fx_and_dfx(const double zbar0, const double *area, const double *f, const double *k,
	const double *Lambda, const double q, const long int ns, double &fx, double &dfx, const long int dt)
{
	int i;
	double temp1;

	fx  = 0.0;
	dfx = 0.0;
    for (i = 1; i <= ns; i++) {
		temp1 = area[i-1]*k[i-1]*exp(-Lambda[i-1])*exp(-f[i-1]*zbar0);
		fx += temp1/f[i-1];
		dfx -= temp1;
    }
	// at this point q is in units of mm**3/ts but fx is in units of m*mm**2/s
	// so divide the flow by the interval length to get compatible units
	fx *= 1000.0/3600.0 - q/dt;
	dfx *= 1000.0/3600.0;

	return 0;
}
