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

#ifndef TDIMS_HH
#define TDIMS_HH

const int MAX_NTDH       = 450;		// Maximum number of values for the Time Delay Histogram used by TOPMOD
const int END_STEPS_MAX  = 205;		// Maximum number of time steps used in forecasting
const int INI_ARRAY_SIZE = 100;		// The initial size allocation to runoff routing arrays - these
									//  arrays are self expanding
const int Nsi  =  3;				// Number of initial condition parameters
const int Nsp  = 39;				// Number of basin parameters
const int Nrp  =  4;				// Number of channel parameters
const int Nip1 = 15;				// The size of the IRR array used for topmodel output - DGT.
const int dpm  = 17;				// The number of parameters that can be assigned by nlfit during calibration.
									//  It should correspond to nlfit parameter nmx
const int nooksack = 1;
const int num_basinpars = 46;

// WARNING - make sure nrx in NLFIT.INC is at least the maximum mnumber of time steps
// in a simulation + 200

#endif
