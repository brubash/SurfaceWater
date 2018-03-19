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
#include "types.hh"

using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;
using namespace Eigen;
using namespace std;

int isort(const UserSourceTableType *table, const int n, int *s);

int AssignPriorityOrder(const int NumUser, const int NumSource, const int NumRights, const int NumReservoir,
	const int AllocationMode, int *DrainageOrder, const int NumDrainage,
	int &NumUserSource, int &NumUserSourceReturn, const int NumReturnFlow)
{
	int latestdate, nrfound, i, j, k, m, n, j_rights, j_src, j_res, j_ret, ID = -999;
	int ifoundresfill, ifoundtemp, nfound, *ifound;
	bool any_positive;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> AssignPriorityOrder(" << ncalls << ")" << std::endl;
    }
	caller = "AssignPriorityOrder";
#endif
	//  code to initialize latest date to avoid uninitialized variable
	if(NumRights > 0) {
		latestdate = Rights[0].PriorityDate;
	} else {
		latestdate = 0;
	}

	//  Code to find the latest Priority date so that users without a priority date specified get a later date
	nrfound = 0;
	any_positive = false;
	for (i = 0; i < NumUser; i++) {
		if (User[i-1].RightID > 0) {
			any_positive = true;
		}
	}
	for (i = 1; i <= NumUser; i++) {
		for (j = 1; j <= User[i-1].NumSources; j++) {
			//for each User-Source combination, find the priority date associated with the water right if it exists
			if (NumRights > 0 && any_positive ) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumRights];
				ifound[nfound] = 0;
				for (n = 0; n < NumRights; n++) {
					if (Rights[n].RightID == User[i-1].RightID[j-1]) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				j_rights = ifound[0];	//problem if nfound<>1
				delete [] ifound;
				if (nrfound == 0) {
					latestdate = Rights[j_rights-1].PriorityDate;
				} else {
					if (latestdate < Rights[j_rights-1].PriorityDate)
						latestdate = Rights[j_rights-1].PriorityDate;
				}
				nrfound++;
			}
		}
	}
	//  End DGT modifications to find latestdate


	//Build a UserSource Table
	NumUserSource = 0;
	NumUserSourceReturn = 0;
	for (i = 1; i <= NumUser; i++) {
		for (j = 1; j <= User[i-1].NumSources; j++) {
			//for each User-Source combination, we want to identify the user,
			//the source, the date, and the drainage the water is coming from
			NumUserSource++;
			UserSourceTable[NumUserSource-1].UserID = i;
			UserSourceTable[NumUserSource-1].SourceCounter = j;
			if (NumRights > 0 && any_positive ) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumRights];
				ifound[nfound] = 0;
				for (n = 0; n < NumRights; n++) {
					if (Rights[n].RightID == User[i-1].RightID[j-1]) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				j_rights = ifound[0];	//problem if nfound<>1
				delete [] ifound;
				UserSourceTable[NumUserSource-1].PriorityDate = Rights[j_rights-1].PriorityDate;
			} else {
				UserSourceTable[NumUserSource-1].PriorityDate = latestdate + 1;  //  DGT 9/4/05  Changed from 0.
				// Where priority date does not exist give it a date later than all where one does exist
			}
			if (User[i-1].UsersType == ReservoirFillUseCode) { //RAW 31-Aug-2005 we always want to give top priority to filling instream reservoirs
				UserSourceTable[NumUserSource-1].PriorityDate = -1;
			}
			// find1()
			nfound = 0; //none found
			ifound = new int[NumSource];
			ifound[nfound] = 0;
			for (n = 0; n < NumSource; n++) {
				if (Source[n].SourceID == User[i-1].SourceID[j-1]) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			j_src = ifound[0];	//problem if nfound<>1
			delete [] ifound;

			switch (Source[j_src-1].Type) {
				case StreamSourceCode:
				case GroundwaterSourceCode: //water is coming from stream or groundwater
					ID = Source[j_src-1].SourceLocationID; //DrainageId is stored with the Source
					if (ID == 0)
						ID = User[i-1].POU_ID;
					break;
				case ReservoirSourceCode:	 //water is coming from reservoir
					// find1()
					nfound = 0; //none found
					ifound = new int[NumReservoir];
					ifound[nfound] = 0;
					for (n = 0; n < NumReservoir; n++) {
						if (Reservoir[n].ReservoirID == Source[j_src-1].SourceLocationID) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					j_res = ifound[0];	//problem if nfound<>1
					delete [] ifound;
					ID = Reservoir[j_res-1].DrainageID; //DrainageID is with the reservoir
					break;
				default:;
			}
			UserSourceTable[NumUserSource-1].DrainageID = ID;
			if (User[i-1].ReturnFlowID >= 0) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumReturnFlow];
				ifound[nfound] = 0;
				for (n = 0; n < NumReturnFlow; n++) {
					if (ReturnFlow[n].ReturnFlowID == User[i-1].ReturnFlowID) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				j_ret = ifound[0];	//problem if nfound<>1
				delete [] ifound;
				NumUserSourceReturn = NumUserSourceReturn + ReturnFlow[j_ret-1].NumReturnFlows;
			}
			//count these return flow so we can dimension an array later
		}
	}

	//Order that UserSource Table
	if (AllocationMode == WaterRightsAllocationCode) { //sort by date
		isort(UserSourceTable, NumUserSource, UserSourceOrder);
	} else if (AllocationMode == DemandAllocationCode || AllocationMode == NoAllocationCode) { //sort by drainage order
		k = 0;
		for (i = 1; i <= NumDrainage; i++) { // we already know what order the drainages are in, so just use it
			// find1()
			nfound = 0; //none found
			ifound = new int[NumUserSource];
			ifound[nfound] = 0;
			for (n = 0; n < NumUserSource; n++) {
				if (UserSourceTable[n].DrainageID == DrainageOrder[i-1]) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			// extra magic to make sure that (in-stream) reservoirs are filled as highest priority
			if (nfound != 0) {
				ifoundresfill = 0;
				m = 1;
				//this loop will hunt through all the users in this drainage, looking for the ReservoirFillUseCode
				//when it finds that ReservoirFillUseCode, it will switch it with the first user in ifound, so that it comes before anything else in this drainage
				//the do while loop ends as soon as this condition is satisfied
				while (ifoundresfill == 0 && m < nfound) {
					if (User[UserSourceTable[ifound[m-1]-1].UserID-1].UsersType == ReservoirFillUseCode ) {
						ifoundresfill = 1;
						ifoundtemp  = ifound[m-1];
						ifound[m-1] = ifound[0];
						ifound[0]   = ifoundtemp;
					}
					m++;
				}

				for (m = 1; m <= nfound; m++) {
					UserSourceOrder[k+m-1] = ifound[m-1];
				}
			}
			delete [] ifound;
			k += nfound;
		}
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving AssignPriorityOrder(" << ncalls << ")" << "\n\n";
	}
	ncalls++;
#endif

	return 0;
}

int isort(const UserSourceTableType *table, const int n, int *s)
{
	Array<int,Dynamic,1> i(n);
	int j, k, m, itemp;

	for (m = 0; m < n; m++) {
		i(m) = table[m].PriorityDate;
	}
	for (j = 1; j <= n; j++) {
		s[j-1] = j;
	}

	for (j = 1; j <= n-1; j++) {
		for (k = j+1; k <= n; k++) {
			if (i(j-1) >= i(k-1)) {
				itemp = i(k-1);
				i(k-1) = i(j-1);
				i(j-1) = itemp;
				itemp  = s[k-1];
				s[k-1] = s[j-1];
				s[j-1] = itemp;
			}
		}
	}

	return 0;
}

