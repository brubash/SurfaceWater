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

#ifndef NLFIT_HH
#define NLFIT_HH

// These parameters set the dimensions for all arrays used in programs
//    NLFIT, RESPONS, PREDICT, COMPAT and EDPMF.

// itb setup for nooksack
const int Nmx   =    17;	// Maximum number of model parameters
const int Nrx   = 22000;	// Number of observations in longest response record
const int iex   =   500;	// Maximum number of responses used in fitting of model
const int mArMa =     2;	// Max number of AR or MA parameters

//      Probailistic search parameters used by NLFIT only:
const int nPopX    = 600;	// Maximum number of populations in GA or random starts
const int maxBitX  =  16;	// Maximum number of bits in GA string
const int maxBestX = 150;	// Maximum number of good solutions to be saved

#endif
