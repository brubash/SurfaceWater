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
using namespace Eigen;
using namespace std;

int AssignDrainageFlows(const int Timestep, const int NumDrainage, const int NumNode, const int NumLink,
	int *DrainageOrder, const int NumRunoff, const int NumBaseflow, ArrayXd &DrainageOutFlow)
{
	int i, n, k, isn, nfound, *ifound;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> AssignDrainageFlows(" << ncalls << ")" << std::endl;
    }
	caller = "AssignDrainageFlows";
#endif
	//Link(k).Flow=Drainage(1:NumDrainage).Flow
	//for each drainage, there is a runoff amount.
	//find the link which connects the drainage node to its stream node
	//assign the runoff to that link
	for (i = 1; i <= NumDrainage; i++) {
		// find2()
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == StreamNodeCode && Node[n].DrainageID == i) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		isn = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		// find2()
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).LinkCode == SurfaceRunoffLinkCode && Link(n).DSNode == isn) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		k = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		Link(k-1).Flow = Runoff[Timestep-1].Rate[Drainage(i-1).DrainageID-1] - Baseflow[Timestep-1].Rate[Drainage(i-1).DrainageID-1];
		// find2()
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).LinkCode == SubsurfaceRunoffLinkCode && Link(n).DSNode == isn) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		k = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		Link(k-1).Flow = Baseflow[Timestep-1].Rate[Drainage(i-1).DrainageID-1];
	}

	BalanceFlowsAtStreamNodes(NumNode, NumLink, DrainageOrder, NumDrainage, DrainageOutFlow);
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving AssignDrainageFlows(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}
