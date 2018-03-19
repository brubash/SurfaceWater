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
#include <vector>
#include <algorithm>

using namespace input_structures;

int SetOrder(int &k, int *DrainageOrder, int &Position, const int NumDrainage);

int BuildDrainageOrder(const int NumDrainage, int *DrainageOrder)
{
	int Position, i;

	for (i = 0; i < NumDrainage; i++) {
		DrainageOrder[i] = 0;
	}
	Position = 0;

    i = 1;
	while (Position < NumDrainage) {
		while (DrainageOrder[i-1] > 0) {	// Seek to the first unordered position
			i++;
		}
		// find first unordered drainage
        SetOrder(i, DrainageOrder, Position, NumDrainage);
	}
    /*
    for (i = 0; i < NumDrainage; i++) {
        cout << dec << setw(4) << DrainageOrder[i];
        if ((i+1)%10 == 0)
            cout << "\n";
    }
    cout << "\n\n";
	for (i = 0; i < NumDrainage; i++) {
        cout << dec << setw(4) << Drainage(i).DSDrainage;
        if ((i+1)%10 == 0)
            cout << "\n";
    }
    cout << endl;
	exit(0); */
	return 0;
}


int SetOrder(int &k, int *DrainageOrder, int &Position, const int NumDrainage)
{
	static int i;
	int n;
	bool skip_step = false;

	if (DrainageOrder[k-1] > 0) {
		std::cerr << "attempt to order a node twice\n";
		return 0;
	} else {
		// find an unlisted entry in the downcatachments column
		i = k;
        do {
            Find(Drainage,i,NumDrainage);
            i++;
        } while (nfound > 0);   // until a drainage has no entry in the downcatchments vector

		// Record it
		if (find_counts.size() > 0 && nfound == 0) {
            find_counts.erase(find_counts.end()-1); // the next-to-last element is the one needed because nfound == 0
            if (find_counts.back() == 0 && find_counts.size() > 0) {
                find_counts.erase(find_counts.end()-1);
            }
            nfound_last = find_counts.back();

            if (find_counts.size() > 0) {
                if ( std::find(done.begin(), done.end(), ifound[global_found-nfound_last]) == done.end()) {   // if not set
                    DrainageOrder[Position] = ifound[global_found-nfound_last];
                    done.push_back(ifound[global_found-nfound_last]);
                    Position++;
                    skip_step = false;
                }

            } else {    //if (find_counts.size() == 0) {
                for (n = 1; n <= NumDrainage; n++) {
                    if ( std::find(done.begin(), done.end(), n) == done.end()) {    // if not set
                        DrainageOrder[Position] = n;
                        done.push_back(n);
                        Position++;
                        if ( std::find(done.begin(), done.end(), n+1) == done.end()) {    // if next not set
                            Find(Drainage,n,NumDrainage);

                            if (find_counts.back() == 0)
                                find_counts.erase(find_counts.end()-1);
                            nfound_last = find_counts.back();

                            break;  // n-loop
                        }
                    }
                }
            }

            if (!skip_step) {
                if (ifound.size() > 0) {
                    ifound.erase (ifound.begin()+global_found-nfound_last);
                    find_counts.back()--;
                    global_found--;
                    if (nfound_last > global_found)
                        nfound_last = find_counts.back();
                }
            }

            if (find_counts.size() > 0 && global_found > 0) {
                if (find_counts.back() == 0)
                    find_counts.erase(find_counts.end()-1);
            }
        }

        if (ifound.size() > 0) {
            k = ifound[global_found-nfound_last+1];
        } else {
            k = ifound[0];
        }
	}

	return 0;
}
