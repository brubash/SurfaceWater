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

using namespace std;
using namespace data_array;

int read_struct_from_text(const string fileName, int &expected_numcols, const int ncommentlines, int &nrows)
{
	string title, headers, inLine, checkLine;
	int i, k;
	ifstream inFile;
	stringstream value;
	bool csvFlag;
	size_t pos;

	inFile.open(fileName.c_str());
	if (!inFile.is_open()) {
		cerr << "Failed to open " << fileName << '\n';
		exit(EXIT_FAILURE);
	}
	if (ncommentlines >= 1)
		getline(inFile, title, '\n');
	if (ncommentlines >= 2)
		getline(inFile, headers, '\n');
	nrows = 0;
	while (inFile.good()) {
		getline(inFile, inLine, '\n');
		if (inLine.length() > 1) {		// We can't use empty lines at the end of a file to look for commas.
			checkLine = inLine;			// so we need to save the last extended line.
			nrows++;
		} else {
			break;
		}
	}

	// Check for commas.
	if ( ( pos=checkLine.find_first_of(',')) != string::npos ) {	// if csv
		csvFlag = true;
	} else {
		csvFlag = false;
	}
	if (nrows > 0) {
		// rewind
		inFile.clear();
		inFile.seekg(0);
		if (ncommentlines >= 2) {
			if (expected_numcols < 0) {
				inFile >> expected_numcols;
			} else {
				getline(inFile, title, '\n');
			}
		}
		real_array.resize(nrows, expected_numcols);

		getline(inFile, headers, '\n');
		for (k = 0; k < nrows; k++) {
			for (i = 0; i < expected_numcols-1; i++) {
				if ( csvFlag ) {
					getline(inFile, inLine, ',');
					value.str(inLine);
					real_array(k,i) = strtod(value.str().c_str(), NULL);
					value.clear();
					value.str("");
				} else {
					inFile >> real_array(k,i);
				}
			}
			if ( csvFlag ) {
				getline(inFile, inLine, '\n');  // look for line end instead of comma
				value.str(inLine);
				real_array(k,expected_numcols-1) = strtod(value.str().c_str(), NULL);
				value.clear();
				value.str("");
			} else {
				inFile >> real_array(k,expected_numcols-1);
				getline(inFile, inLine, '\n');  // read line end or inline comment
			}
		}
	}
	inFile.close();

	return 0;
}

	// ==================================================

int read_bdryflow(const string fileName, int &expected_numcols, const int ncommentlines, int &nrows)
{
	string title, headers, inLine;
	int i, k;
	ifstream inFile;

	inFile.open(fileName.c_str());
	if (!inFile.is_open()) {
		cerr << "Failed to open " << fileName << '\n';
		exit(EXIT_FAILURE);
	}
	getline(inFile, title, '\n'); //"This file provides ..."
	getline(inFile, title, '\n'); //"Flow values are  ..."
	getline(inFile, title, '\n'); //"*12205000*HNW-2052*12208000*"
	integer_array.resize(100,1);

	inFile >> title >> expected_numcols;
	for (i = 0; i < expected_numcols; i++) {
		inFile >> integer_array(i,0);	// "Ver2 3  8  18  21 Date Hour"
	}
	expected_numcols += 2;

	nrows = 0;
	getline(inFile, inLine, '\n');
	while (inLine.length() > 1 && !inFile.eof()) {
		getline(inFile, inLine, '\n');
		if (inLine.length() <= 1) {
			break;
		}
		nrows++;
	}
	if (nrows > 0) {
		// rewind
		inFile.clear();
		inFile.seekg(0);
		getline(inFile, title, '\n');
		getline(inFile, title, '\n');
		getline(inFile, title, '\n');
		getline(inFile, title, '\n');
		dble_array = new double*[nrows];
		for (k = 0; k < nrows; k++) {
			dble_array[k] = new double[expected_numcols];
		}
		for (k = 0; k < nrows; k++) {
			for (i = 0; i < expected_numcols; i++) {
				inFile >> dble_array[k][i];
			}
		}
	}
	inFile.close();

	return 0;
}

