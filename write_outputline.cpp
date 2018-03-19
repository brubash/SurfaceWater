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

using namespace std;
using namespace Eigen;

//  Subroutine to write output at each time step
//  David Tarboton  6/29/05
int Write_OutputLine(ofstream &oFile, const string fileName, const int timestep, const double *Rvariable, const int NumDrainage, const double scalefactor)
{
	string LocationTypeString, str;
	int i, j;

	if (timestep == 0) {  // Initialize
		oFile.open(fileName.c_str());
		if (!oFile.is_open()) {
			cerr << "Failed to open " << fileName << '\n';
			exit(EXIT_FAILURE);
		} //else {
			//cerr << '\'' << fileName << "' opened" << endl;
		//}
		LocationTypeString = "Drainage";
		oFile << setw(9) << "TimeStep     ";
		for (i = 1; i <= NumDrainage; i++) {
			oFile << setw(12) << LocationTypeString << dec << setw(3) << i;
		}
		oFile << '\n';
	} else if (timestep > 0) {
		if (!oFile.is_open()) {
			cerr << fileName << " is not open\n";
			exit(EXIT_FAILURE);
		}
		oFile << dec << setw(12) << timestep;
		for (j = 0; j < NumDrainage; j++) {     // remove branch below (remove fixed) in the release version
            if (fabs(Rvariable[j]*scalefactor) < 1.0e-99) {
                oFile << fixed << setw(15) << setprecision(5) << 0.0;
            } else {
                if (fabs(Rvariable[j]*scalefactor) < 0.1) {
                    oFile << scientific << setw(15) << setprecision(7) << Rvariable[j]*scalefactor;
                } else {
                    oFile << fixed << setw(15) << setprecision(5) << Rvariable[j]*scalefactor;
                }
            }
        }
		oFile << '\n';

	} else {	// timestep < 0
		oFile.close();
	}

	return 0;
}

using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;

//   Subroutine added to output local node drainage contributions
int Write_OutputLocalContributions(ofstream &o1File, ofstream &o2File, const int NumStreamNode, const int NumDrainage,
		const ArrayXd &BaseflowTopnet, const ArrayXd &RunoffTopnet, const int timestep, const double scalefactor)
{
	double LocalSurfaceRunoff, LocalBaseflow, fracrunoff;
	int j, iddrain;

	if (timestep == 0){  // Initialize
		o1File.open("results/NodeLocalSurfaceRunoff_cms.txt");
		if (!o1File.is_open()) {
			cerr << "Failed to open 'NodeLocalSurfaceRunoff_cms.txt'\n";
			exit(EXIT_FAILURE);
		} //else {
		//	cerr << "'NodeLocalSurfaceRunoff_cms.txt' opened" << endl;
		//}
		o1File << "TimeStep  Node_ID  Flow_cms\n";
		o2File.open("results/NodeLocalBaseflow_cms.txt");
		if (!o2File.is_open()) {
			cerr << "Failed to open 'NodeLocalBaseflow_cms.txt'\n";
			exit(EXIT_FAILURE);
		} //else {
			//cerr << "'NodeLocalBaseflow_cms.txt' opened" << endl;
		//}
		o2File << "TimeStep  Node_ID  Flow_cms\n";
	} else if (timestep > 0) {
		for (j = 1; j <= NumStreamNode; j++) {
			iddrain = StreamNode[j-1].DrainageID;

			fracrunoff = StreamNode[j-1].LocalArea/Drainage(iddrain-1).CatchArea*1.0e6;  // Drainage(iddrain)%CATCHAREA is in mm^2
			LocalSurfaceRunoff = (RunoffTopnet(iddrain-1) - BaseflowTopnet(iddrain-1))*fracrunoff*scalefactor;
			LocalBaseflow = BaseflowTopnet(iddrain-1)*fracrunoff*scalefactor;
			o1File << dec << setw(8) << timestep << dec << setw(8) << StreamNode[j-1].NodeID;
			o1File << scientific << setw(15) << setprecision(7) << LocalSurfaceRunoff << '\n';
			o2File << dec << setw(8) << timestep << dec << setw(8) << StreamNode[j-1].NodeID;
			o2File << scientific << setw(15) << setprecision(7) << LocalBaseflow << '\n';
		}
	} else {	// timestep < 0
		o1File.close();
		o1File.close();
	}
	return 0;
}

//========================= Eigen =========================

int Write_OutputLine_Eigen(ofstream &oFile, const string fileName, const int timestep, const ArrayXd &Rvariable, const int NumDrainage, const double scalefactor)
{
	string LocationTypeString, str;
	int i, j;

	if (timestep == 0) {  // Initialize
		oFile.open(fileName.c_str());
		if (!oFile.is_open()) {
			cerr << "Failed to open " << fileName << '\n';
			exit(EXIT_FAILURE);
		} //else {
			//cerr << fileName << " opened" << endl;
		//}
		LocationTypeString = "Drainage";
		oFile << setw(9) << "TimeStep     ";
		for (i = 1; i <= NumDrainage; i++) {
			oFile << setw(12) << LocationTypeString << dec << setw(3) << i;
		}
		oFile << '\n';
	} else if (timestep > 0) {
		if (!oFile.is_open()) {
			cerr << fileName << " is not open\n";
			exit(EXIT_FAILURE);
		}
		oFile << dec << setw(12) << timestep;
		for (j = 0; j < NumDrainage; j++) {     // remove branch below (remove fixed) in the release version
            if (fabs(Rvariable(j)*scalefactor) < 0.1) {
                oFile << scientific << setw(15) << setprecision(7) << Rvariable(j)*scalefactor;
            } else {
                oFile << fixed << setw(15) << setprecision(5) << Rvariable(j)*scalefactor;
            }
		}
		oFile << '\n';

	} else {	// timestep < 0
		oFile.close();
	}

	return 0;
}
