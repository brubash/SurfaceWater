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
#include <sys/stat.h>
#include <unistd.h>

using namespace std;
using namespace Eigen;

// *******************************************************************
//          SUBROUTINE  RAINFILL_DOIT
// Code to fill in missing rainfall data using regression on other sites
// This new section of code contains two subroutines and two functions
// Designed to be called after rain.dat has been read (e.g. from hyData),
// but before rain data is used.
// -1 in rain.dat is assumed to be the indicator for missing data
// Regression constants are stored in text file to be developed by model user
// File name is hard-coded 'rainfill.txt'.
// If the file is not present, the code doesn't do any filling

int rainfill_doit(ArrayXd &train, const int it, const int Ngauge, ArrayXXd &rcoeff, ArrayXXi &fsite, const int *nfill,
	bool &i_reset_flag, int &it_save_old, const int maxGauge)
{
	int itry, iok, ig, ig1, ifill;
	int it_save;
	int im;
	double train_last[maxGauge];
	bool canfill;

	if (it == 1) {
		it_save = 0;
		i_reset_flag = false;
	}
	//	do 100 it=1,nt !repeat separately for each timestep
	//keep looping around the gauges until they're all filled, or we've tried too many times
	// write(6,*)'Timestep ',it
	itry = 0; //just to get started
	iok = 0;
	while (iok < Ngauge && itry < Ngauge) {
		cerr << "Try # " << itry << '\n';
		itry++;
		iok = 0;
		for (ig = 1; ig <= Ngauge; ig++) {	//check each gauge to see if it's missing
			cerr << "Gauge " << ig << '\n';
			if (miss(train(ig-1), nfill[ig-1])) {
				ifill = 0;					//if this gauge is missing try to fill it using available eqns
				for (im = 1; im <= abs(nfill[ig-1]); im++) {
					ig1 = fsite(ig-1,im-1); //fill it if it's not already done
					if (ig != ig1) { //first decide whether we have some data to fill with
						canfill = !miss(train(ig1-1), nfill[ig1-1]);
					} else {
						if (it > 1) {
							canfill = !miss(train_last[ig-1], nfill[ig-1]);
						} else {
							canfill = false;
						}
					}
					if ((ifill == 0) && canfill) { //if we haven't filled it yet and there's data available at the fill site then do it
						if (ig != ig1) {
							train(ig-1) *= rcoeff(ig-1,im-1);
						} else {
							train(ig-1) = train_last[ig-1]*rcoeff(ig-1,im-1);
						}
						ifill = 1;
						//write(errlun,*)'ig,im,ig1,train(ig1,it),rcoeff(ig-1,im-1),train(ig,it)'
						//write(errlun,*)ig,im,ig1,train(ig1,it),rcoeff(ig-1,im-1),train(ig,it)
					}
				}
				iok = iok+ifill;
			} else {
				iok++;
				train_last[ig-1] = train(ig-1);
			}
		}
	}
	if (iok == 0) {
		cerr << "Warning: all gauges missing data: step " << it << '\n';
		// introduced the following code to reduce the length of the forecast
		// window when there is missing data at all gauges from time step IT_OLD_SAVE to
		// the end of the rainfall data
		if (it == it_save+1 && i_reset_flag == false) {
			it_save_old = it_save;
			i_reset_flag = true;
		}
		for (ig = 1; ig <= Ngauge; ig++) {
			if (it > 1) {
				train(ig-1) = train_last[ig-1];
			} else {
				train(ig-1) = 0;
			}
		}
	} else {
		i_reset_flag = false;
	}
	it_save = it;  // RPI 3/3/03
	//100	continue
	return 0;
}

bool miss(const double train, const int nfill)
{
	bool result;

	if (train < 0 || (train == 0 && nfill < 0) ) {
		result = true;
	} else {
		result = false;
	}
	return result;
}


// *****************************************************************
//     SUBROUTINE  RAINFILL_READ
// ******************************************************************

int rainfill_read(const string rainfill_file, const int Ngauge, const int *rsite,
	ArrayXXd &rcoeff, ArrayXXi &fsite, int *nfill, int &fillflag, const int maxGauge)
{

	int ig, is[maxGauge], i, ierr, j;
	int nfill0;
	int fsite0[maxGauge];
	int msite;
	double rcoeff0[maxGauge];
	struct stat filestatus;

	fillflag = -1;

	ifstream rainfillFile(rainfill_file.c_str());		// fortran unit lun
	if (!rainfillFile.is_open()) {
		cerr << "Failed to open " << rainfill_file << '\n';
		exit(EXIT_FAILURE);
	}
    stat( rainfill_file.c_str(), &filestatus );
    if (filestatus.st_size < 2) {
        return -1;
    }

	fillflag = 1;
	for (j = 1; j <= Ngauge; j++) {
		ierr = 0;
		rainfillFile >> msite >> nfill0;
		for (i = 1; i <= abs(nfill0); i++) {
			rainfillFile >> rcoeff0[i-1] >> fsite0[i-1];
		}
		ig = isrch(msite, rsite, Ngauge);
		if (ig == 0) {
			cerr << "Rainfill: Site " << msite << " not in rain.dat\n";
			//write(errlun,*)'Rainfill: Site ',msite,' not in rain.dat'
			ierr = 1;
		}
		for (i = 1; i <= abs(nfill0); i++) {
			is[i-1] = isrch(fsite0[i-1], rsite, Ngauge);
			if (is[i-1] == 0) {
				cerr << "Rainfill: Site " << fsite0[i-1] << " not in rain.dat\n";
				//write(errlun,*)'Rainfill: Site ',fsite0(i),' not in rain.dat'
				ierr = 1;
			}
		}

		if (ierr == 0) {
			nfill[ig-1] = nfill0;
			for (i = 1; i <= abs(nfill0); i++) {
				rcoeff(ig-1,i-1) = rcoeff0[i-1];
				fsite(ig-1,i-1)  = is[i-1];
			}
		}
	}
	rainfillFile.close();

	return 0;
}


int isrch(const int mSite, const int *rSite, const int nGauge)
{
	int result, ig;

	result = 0;
	for (ig = 1; ig <= nGauge; ig++) {
		if (result == 0 && mSite == rSite[ig-1]) {
			result = ig;
		}
	}

	return result;
}

