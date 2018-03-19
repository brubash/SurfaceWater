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

using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;
using namespace TimeVaryingOutput;
using namespace std;

int Initialise_Output_Tables(const int NumDrainage, const int NumNode, const int NumStreamNode, const int NumLink, const int NumUser,
	const int NumTimesteps, const int NumReservoir, const int NumUserSourceReturn, const int NumReturnFlow)
{
	int n, i, j, k, js, j_src, jr, j_ret, j_ret_ind, nfound, *ifound;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> Initialise_Output_Tables(" << ncalls << ")" << std::endl;
    }
	caller = "Initialise_Output_Tables";
#endif
	FlowInLinks_cms.resize(NumTimesteps,NumLink);
	ReservoirStorage_m3.resize(NumTimesteps,NumReservoir);
	DateTime_yyyymmdd_hhmmss.resize(NumTimesteps, 2);
	FlowAtStreamNodes_cms = new double[NumStreamNode];

	for (i = 1; i <= NumDrainage; i++) {

		// find2()
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == StreamNodeCode && Node[n].SelfID == i) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		j = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		// find2()
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).USNode == j && Link(n).LinkCode == UnallocatedLinkCode) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		k = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		StaticOutput.StreamFlowLinks[i-1].DrainageID = i;
		StaticOutput.StreamFlowLinks[i-1].LinkID = k;
	}

	for (i = 1; i <= NumDrainage; i++) {
		StaticOutput.DrainageID[i-1].TopnetID   = Drainage(i-1).DrainageID;
		StaticOutput.DrainageID[i-1].DrainageID = Drainage(i-1).RealDrainageID;
	}

	k = 0;
	for (i = 1; i <= NumUser; i++) {
		for (js = 1; js <= User[i-1].NumSources; js++) {
			j_src = User[i-1].SourceID[js-1];
			j_ret = User[i-1].ReturnFlowID;
			if (j_ret > 0) {

				// find1()
				nfound = 0; //none found
				ifound = new int[NumReturnFlow];
				ifound[nfound] = 0;
				for (n = 0; n < NumReturnFlow; n++) {
					if (ReturnFlow[n].ReturnFlowID == j_ret) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				j_ret_ind = ifound[0];
				delete [] ifound;

				for (jr = 1; jr <= ReturnFlow[j_ret_ind-1].NumReturnFlows; jr++) {
					k++;
					if (k > NumUserSourceReturn) {
						std::cerr << "too many return flows for the output table\n";
					}
					StaticOutput.DrainageInfo[k-1].UserSrceReturn_ID     = k;
					StaticOutput.DrainageInfo[k-1].User_ID               = i;
					StaticOutput.DrainageInfo[k-1].UsersType             = User[i-1].UsersType;
					StaticOutput.DrainageInfo[k-1].SourceDrainageID      = Source[j_src-1].RealSourceLocationID;
					StaticOutput.DrainageInfo[k-1].DestinationDrainageID = ReturnFlow[j_ret_ind-1].RealReturnFlowsLocn[jr-1];
					StaticOutput.DrainageInfo[k-1].WWTP_ID               = ReturnFlow[j_ret_ind-1].WWTP_ID[jr-1];
					StaticOutput.DrainageInfo[k-1].RF_counter            = jr;
					StaticOutput.DrainageInfo[k-1].SourceID              = j_src;
					StaticOutput.DrainageInfo[k-1].ReturnFlowID          = j_ret;
				}
			}
		}
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving Initialise_Output_Tables(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

