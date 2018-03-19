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

// This version, V2, has been adjusted for inclusion of lake modelling

#include "topnet.hh"

using namespace Eigen;


int iposn(const int Nrch, const ArrayXXi &linkR, const int ii)
{
	// Function that gives the position of the array element which has the reach number ii.
	int i;

	for (i = 0; i < Nrch; i++) {
		if (linkR(0,i) == ii) {
			return i;
		}
	}
	return 0;
}

