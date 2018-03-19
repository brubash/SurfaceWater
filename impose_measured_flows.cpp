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

int ScaleUpstreamRunoff(const double FlowRatio, const int k, const int Timestep, const int NumNode,
	const int NumRunoff, const int NumBaseflow, const int NumDrainage);

int ImposeMeasuredFlows(const int Timestep, const int NumNode, const int NumLink, const int NumRunoff,
	const int NumBaseflow, const int NumDrainage, const int NumMeasuredFlowInfo, const int NumMeasuredFlowData,
	int *DrainageOrder, ArrayXd &DrainageOutFlow)
{
	double NodeOutFlow, FlowRatio;
	int i, n, k, j, j_out, nfound, *ifound;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> ImposeMeasuredFlows(" << ncalls << ")" << std::endl;
    }
	caller = "ImposeMeasuredFlows";
#endif
	// Impose any external measured flows by
	// identifying the node with NodeType =  3  STREAM (whose ID is NODE-ID) at which the flow was measured,
	// noting the ratio of measured to modelled flow at this node, and then
	// applying that ratio to Runoff and Baseflow rates for this timestep in any
	// drainage which is upstream of this node (including the drainage the node lies in)

	//This is different to the previous method I proposed, which just fiddled
	//the flow from the measuredflownode to this node until the total inflow to
	//this node was as measured
	// I changed methods so that we could get the runoff from all catchments set
	// properly, because the water quality folks need runoff&baseflow from ALL drainages
	for (i = 1; i <= NumMeasuredFlowInfo; i++) {
		// find2()
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == StreamNodeCode && Node[n].DrainageID == MeasuredFlowInfo[i-1].DrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		k = ifound[0]; //problem if nfound<>1 warning if nfound=0
		delete [] ifound;

		// find1()
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).USNode == k) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		NodeOutFlow = 0.0;
		for (j = 1; j <= nfound; j++) {
			j_out = ifound[j-1];
			NodeOutFlow += Link(j_out-1).Flow;
		}
		delete [] ifound;
		//  added the check for NodeOutFlow Greater than 0 to avoid infinity and Nan in the results
		if (MeasuredFlowData[Timestep-1].Flow[MeasuredFlowInfo[i-1].ColInMeasFlow-1] >= 0.0 && NodeOutFlow > 0.0) {
			FlowRatio = MeasuredFlowInfo[i-1].ScalingFactor*MeasuredFlowData[Timestep-1].Flow[MeasuredFlowInfo[i-1].ColInMeasFlow-1]/NodeOutFlow;
		} else {
			FlowRatio = 1.0;
		}
		ScaleUpstreamRunoff(FlowRatio, k, Timestep, NumNode, NumRunoff, NumBaseflow, NumDrainage);
	}

	//now use these new runoffs to calculate flows throughout the network
	AssignDrainageFlows(Timestep, NumDrainage, NumNode, NumLink, DrainageOrder, NumRunoff, NumBaseflow, DrainageOutFlow);
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving ImposeMeasuredFlows(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}


int ScaleUpstreamRunoff(const double FlowRatio, const int k, const int Timestep, const int NumNode,
	const int NumRunoff, const int NumBaseflow, const int NumDrainage)
{
	int i[10], kk[10];
	int j, j_node, m, n, ni, nk, ii, number_processed, list_length, nfound, *ifound;
	Array<int,Dynamic,1> nodelist(NumNode);

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> ScaleUpstreamRunoff(" << ncalls << ")" << std::endl;
    }
	caller = "ScaleUpstreamRunoff";
#endif
	nodelist(0) = k;
	list_length = 1;
	number_processed = 0;
	while (number_processed < list_length) {
		j_node = nodelist(number_processed);
		j = Drainage(Node[j_node-1].DrainageID-1).DrainageID;
		Runoff[Timestep-1].Rate[j-1] = Runoff[Timestep-1].Rate[j-1]*FlowRatio;
		Baseflow[Timestep-1].Rate[j-1] = Baseflow[Timestep-1].Rate[j-1]*FlowRatio;
		number_processed++;
		// find1() Drainage(i) is  just upstream of Node(k)
		nfound = 0; //none found
		ifound = new int[NumDrainage];
		ifound[nfound] = 0;
		for (n = 0; n < NumDrainage; n++) {
			if (Drainage(n).DSDrainage == Node[j_node-1].DrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		if (nfound > 0) {
			for (m = 0; m < nfound; m++) {
				 i[m] = ifound[m];
			}
			delete [] ifound;
			ni = nfound;
			nk = 0;
			for (ii = 1; ii <= ni; ii++) { //kk will contain list of node indices
				// find2()
				nfound = 0; //none found
				ifound = new int[NumNode];
				ifound[nfound] = 0;
				for (n = 0; n < NumNode; n++) {
					if (Node[n].Type == StreamNodeCode && Node[n].DrainageID == Drainage(i[ii-1]-1).DrainageID) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				kk[ii-1] = ifound[0];
				delete [] ifound;
				nk += nfound;
			}
			for (j = 1; j <= nk; j++) {
				nodelist(list_length+j-1) = kk[j-1];
			}
			list_length += nk;
		} else {
			delete [] ifound;
		}
	}	// while
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving ScaleUpstreamRunoff(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

