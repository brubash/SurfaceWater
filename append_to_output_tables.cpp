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

using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;
using namespace TimeVaryingOutput;
using namespace std;

int Append_To_Output_Tables(const int Timestep, const int dt,
	const int NumStreamNode, const int yyyymmdd, const int hhmmss,
	const int NumLink, const int NumNode, const int Nsteps, const int NumUser)
{
	int ius, ids, i, n, j, nfound, *ifound;
	double ThisDrainageOutflow;

	string LocationTypeString;
	string fmtstr;
	string str;
	static ofstream outFile;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> Append_To_Output_Tables(" << ncalls << ")" << std::endl;
    }
	caller = "Append_To_Output_Tables";
#endif
	DateTime_yyyymmdd_hhmmss(Timestep-1, 0) = yyyymmdd; //Timestep	Date	Time
	DateTime_yyyymmdd_hhmmss(Timestep-1, 1) = hhmmss; //Timestep	Date	Time

	for (i = 1; i <= NumLink; i++) {
		FlowInLinks_cms(Timestep-1,i-1) = Link(i-1).Flow/(double)(dt); //m3/s = m3/step / (s/step)
	}

	// considered writing outputs as we go here to avoid having to save all the time steps
	// but FlowInLinks is used in a post processing calculation of WWTP outflows in Write_Output_tables
	// so can not do this without restructuring that.

	// find1()
	nfound = 0; //none found
	ifound = new int[NumNode];
	ifound[nfound] = 0;
	for (n = 0; n < NumNode; n++) {
		if (Node[n].Type == ReservoirNodeCode) {
			nfound++;
			ifound[nfound-1] = n+1;
		}
	}

	for (i = 1; i <= nfound; i++) {
		ReservoirStorage_m3(Timestep-1,Node[ifound[i-1]-1].SelfID-1) = Node[ifound[i-1]-1].Store;
	}
	delete [] ifound;
	// Flow through stream nodes
	// first the drainage outlets

	// new routing method

	// there are two components to flow at any node
	// 1. the local runoff - use FracRunoff (see read_inputs) to calculate this (we do not route this water from node to node)
	// 2. flow(s) from drainages which flow directly into the drainage containing this node, and whose flowpaths do pass through this node

	for (i = 0; i < NumStreamNode; i++) {
		FlowAtStreamNodes_cms[i] = 0.0;
	}
	for (i = 1; i <= NumStreamNode; i++) {
		//   if(i.eq.107)then
		//   FlowAtStreamNodes_cms(1)=FlowAtStreamNodes_cms(1)
		//   endif
		if (StreamNode[i-1].DOutFlag == 1) { //we need to propagate  this flow to any internal stream nodes
			//  DGT 8/21/05 replaced DrainageOutFlow by Link%Flow because DrainageOutFlow is lacking management
			//  effects like diversions
			ThisDrainageOutflow = Link(StreamNode[i-1].DrainageID*3-1).Flow/(double)dt;  // DGT 8/21/05 Rely on the fact that LinkID is 3*DrainageID
			//RAW 29-Aug-2005 look for user that is InStreamReservoirReleaseUseCode
			//This is the user that sends spill flow and environmental releases to the downstream drainage
			//To make this code faster, do this find2 ahead of time in the first call to watermgmt, and store the results
			//or check for the existence of an instream reservoir located in this drainage
			// find2()
			nfound = 0; //none found
			ifound = new int[NumUser];
			ifound[nfound] = 0;
			for (n = 0; n < NumUser; n++) {
				if (User[n].UsersType == InStreamReservoirReleaseUseCode && User[n].POU_ID == StreamNode[i-1].DrainageID) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			if (nfound == 1) {
				ThisDrainageOutflow += User[ifound[0]].Withdrawal[Timestep-1]/(double)dt;
			}
			delete [] ifound;

			// find2()
			nfound = 0; //none found
			ifound = new int[NumUser];
			ifound[nfound] = 0;
			for (n = 0; n < NumUser; n++) {
				if (User[n].UsersType == DownstreamReservoirReleaseUseCode && User[n].POU_ID == StreamNode[i-1].DrainageID) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			for (j = 1; j <= nfound; j++) {
				ThisDrainageOutflow += User[ifound[j-1]-1].Withdrawal[Timestep-1]/(double)dt;
			}
			delete [] ifound;

			// Code to also append Instreamflowuse
			// find2()
			nfound = 0; //none found
			ifound = new int[NumUser];
			ifound[nfound] = 0;
			for (n = 0; n < NumUser; n++) {
				if (User[n].UsersType == InstreamFlowUseCode && User[n].POU_ID == StreamNode[i-1].DrainageID) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			for (j = 1; j <= nfound; j++) {
				ThisDrainageOutflow += User[ifound[j-1]-1].Withdrawal[Timestep-1]/(double)dt;
			}
			delete [] ifound;

				//		FlowAtStreamNodes_cms(i)=DrainageOutFlow(StreamNode[i-1].DrainageID)/float(dt)
			FlowAtStreamNodes_cms[i-1] = ThisDrainageOutflow;
			//now propagate this flow thru stream nodes down to next drainage outlet or end of network
			ius = i;
				//		call find1(StreamNode%NodeId,StreamNode(ius)%DownNodeID,NumStreamNode) !would be more efficient to do this find back in watermgmt
				//		ids=ifound(1) !now start propagating downstream until we reach next drainage outlet, or end of network
			ids = StreamNode[ius-1].DownNodeID;  // find not necessary.  Rely on node numbers in nodelinks.txt being a counting sequence starting at 1
			if (ids >= 1) {  // avoid addressing error
				while (StreamNode[ius-1].DownNodeID != -1 && StreamNode[ids-1].DOutFlag != 1) {
					FlowAtStreamNodes_cms[ids-1] += ThisDrainageOutflow;
					//	DrainageOutFlow(StreamNode[i-1].DrainageID)/float(dt)
					ius = ids;
					//			call find1(StreamNode%NodeId,StreamNode(ius)%DownNodeId,NumStreamNode)  !would be more efficient to do this find back in watermgmt
					//			ids=ifound(1)
					ids = StreamNode[ius-1].DownNodeID;   // find not necessary
					if(ids < 1)
						break;  // avoid addressing error
				} //while
			}
		} else {// this is an internal streamnode, we need to calculate local runoff to it
				j = StreamNode[i-1].DrainageID; // j indicates the drainage which contains this node
				FlowAtStreamNodes_cms[i-1] += StreamNode[i-1].FracRunoff*Runoff[Timestep-1].Rate[j-1]/(double)dt;
		}
	} //i

	// all the water management calculations take place at the drainage outlet
	// this calculation point is downstream of all "internal" stream nodes
	// this model does not know whereabouts inside the drainage any water mgmt activity occurs
	// thus the values of FlowAtStreamNodes_cms for nodes that are not at drainage outlets will not
	// reflect the within-drainage water management
	// if more spatial detail is requried then a finer drainage discretisation must be used

	// Write the outputs as we go to save having to store all the time steps
	if (Timestep == 1) {
		outFile.open("results/FlowAtStreamNodes_cms.txt");
		LocationTypeString = "Node";
		outFile << setw(12) << "TimeStep";
		for (i = 1; i <= NumStreamNode; i++) {
			outFile << setw(15) << LocationTypeString << i;
		}
		outFile << '\n';
	}
	outFile << dec << setw(12) << Timestep;
	for (j = 0; j < NumStreamNode; j++) {
		outFile << scientific << setw(15) << setprecision(7) << FlowAtStreamNodes_cms[j];
	}
	outFile << '\n';
	if(Timestep == Nsteps) {
		outFile.close();
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving Append_To_Output_Tables(" << ncalls << ")\n" << endl;
    }
    ncalls++;
#endif

	return 0;
}


