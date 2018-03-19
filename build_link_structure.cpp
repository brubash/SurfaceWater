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

using namespace std;
using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;
using namespace Eigen;

int BuildLinkStructure(const int NumDrainage, const int NumUser, const int NumSource,
	const int NumReturnFlow, const int NumReservoir, const int NumMeasuredFlowInfo, const int NumNode, int &NumLink)
{
	int USNode = -999, DSNode = -1, SrcLocnID = -999, RFLocn, RFType, nfound, *ifound;
	int i, j, k, m, n, j_ret, idn = -1, isn, ign, iType, isink, NumReturnFlows;
	//int IntExtCode;
	string Title;
	bool found;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> BuildLinkStructure(" << ncalls << ")" << std::endl;
    }
	caller = "BuildLinkStructure";
#endif
	NumLink = 0;
	//Link=[]
	// Connect drainage nodes and stream nodes
	for (i = 1; i <= NumDrainage; i++) {
		//two links for each drainagenode-streamnode connection, one to carry surface runoff and one for baseflow
		// find2()
		found = false;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == DrainageNodeCode && Node[n].DrainageID == i) {
				//nfound++;
				//ifound[nfound-1] = n+1;
				idn = n+1;
				if (!found) {
					found = true;
				} else {
					cerr << "BuildLinkStructure(): Duplicate link found" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

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

		CreateLink(NumLink, "SRFLOW", SurfaceRunoffLinkCode, ExternalCode, idn, isn, 0.0, 0);

		// find2()
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == GroundwaterNodeCode && Node[n].DrainageID == i) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		ign = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		CreateLink(NumLink,"SBFLOW", SubsurfaceRunoffLinkCode, ExternalCode, ign, isn, 0.0, 0);

		//one link for each streamnode-downstreamnode connection (including sink)
		Title = "UNALLC";
		iType = UnallocatedLinkCode;
		if (Drainage(i-1).DSDrainage >= 1) {
			//IntExtCode = InternalCode;
			// find2()
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == StreamNodeCode && Node[n].DrainageID == Drainage(i-1).DSDrainage) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			DSNode = ifound[0]; //problem if nfound<>1
			delete [] ifound;
		} else {
			//IntExtCode = ExternalCode;
			// find1()
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == SinkNodeCode) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			DSNode = ifound[0];  //problem if nfound<>1
			delete [] ifound;
		}
		//cout << " IntExtCode " << IntExtCode << '\n';
		CreateLink(NumLink, "UNALLC", iType, InternalCode, isn, DSNode, 0.0, 0);	// should this be IntExtCode instead of InternalCode?
	}

	//Connect sources of water to users
	for (i = 1; i <= NumUser; i++) {
		for (j = 1; j <= User[i-1].NumSources; j++) { //one link for each sourcenode-usernode connection
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
			k = ifound[0];  //problem if nfound<>1
			delete [] ifound;
			SrcLocnID = Source[k-1].SourceLocationID;
			if (SrcLocnID == 0)
				SrcLocnID = User[i-1].POU_ID;
			switch (Source[k-1].Type) {
				case (StreamSourceCode):
					// find2()
					nfound = 0; //none found
					ifound = new int[NumNode];
					ifound[nfound] = 0;
					for (n = 0; n < NumNode; n++) {
						if (Node[n].Type == StreamNodeCode && Node[n].SelfID == SrcLocnID) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					USNode = ifound[0]; //problem if nfound<>1
					delete [] ifound;
					break;
				case (GroundwaterSourceCode):
					// find2()
					nfound = 0; //none found
					ifound = new int[NumNode];
					ifound[nfound] = 0;
					for (n = 0; n < NumNode; n++) {
						if (Node[n].Type == GroundwaterNodeCode && Node[n].SelfID == SrcLocnID) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					USNode = ifound[0]; //problem if nfound<>1
					delete [] ifound;
					break;
				case (ReservoirSourceCode):
					// find2()
					nfound = 0; //none found
					ifound = new int[NumNode];
					ifound[nfound] = 0;
					for (n = 0; n < NumNode; n++) {
						if (Node[n].Type == ReservoirNodeCode && Node[n].SelfID == SrcLocnID) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					USNode = ifound[0]; //problem if nfound<>1
					delete [] ifound;
					break;
				default:;
			}
			// find2()
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == UserNodeCode && Node[n].SelfID == i) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			DSNode = ifound[0]; //problem if nfound<>1
			delete [] ifound;
			CreateLink(NumLink, "USERAB", UserAbstractionLinkCode, InternalCode, USNode, DSNode, 0.0, 0);
		}
	} //making links from sources to user nodes

	//one link for each return flow pathway, from the user node to the node where return flows go for this use
	// find1()
	nfound = 0; //none found
	ifound = new int[NumNode];
	ifound[nfound] = 0;
	for (n = 0; n < NumNode; n++) {
		if (Node[n].Type == SinkNodeCode) {
			nfound++;
			ifound[nfound-1] = n+1;
		}
	}
	isink = ifound[0];  //problem if nfound<>1
	delete [] ifound;

	for (i = 1; i <= NumUser; i++) {
		if (User[i-1].ReturnFlowID > 0) {
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
			j_ret = ifound[0];  //problem if nfound<>1
			delete [] ifound;
			NumReturnFlows = ReturnFlow[j_ret-1].NumReturnFlows;
			for (m = 1; m <= NumReturnFlows; m++) {
				// find2()
				nfound = 0; //none found
				ifound = new int[NumNode];
				ifound[nfound] = 0;
				for (n = 0; n < NumNode; n++) {
					if (Node[n].Type == UserNodeCode && Node[n].SelfID == i) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				USNode = ifound[0]; //problem if nfound<>1
				delete [] ifound;
				RFLocn = ReturnFlow[j_ret-1].ReturnFlowsLocn[m-1];
				// DGT 9/3/05  Added 3 lines below to define RFLocn for InstreamFlowUseCode where RFLocn is 0.  This is to allow
				// multiple instreamflow uses to reference the same record in returnflow.txt and not require the downstream drainage
				// to be specified.  This was an attempt to overcome the problem with the select case below, but does not work.
				// For now users are required to specify return flow location in returnflow.txt (and get it right)
				if (RFLocn == 0 && User[i-1].UsersType == InstreamFlowUseCode) {
					RFLocn = Drainage(User[i-1].POU_ID-1).DSDrainage;   // This is the downstream drainge
				}
				if (RFLocn == 0)
					RFLocn = User[i-1].POU_ID;
				RFType = ReturnFlow[j_ret-1].ReturnFlowsType[m-1];
				if (RFType == 0)
					RFType = StreamNodeCode;
				if (RFLocn > 0) {
					// find2()
					nfound = 0; //none found
					ifound = new int[NumNode];
					ifound[nfound] = 0;
					for (n = 0; n < NumNode; n++) {
						if (Node[n].Type == RFType && Node[n].SelfID == RFLocn) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					j = ifound[0]; //problem if nfound<>1
					delete [] ifound;
				} else {
					j = isink;
				}
				//   DGT 9/3/05 abandoned this select case structure.  j is an index for nodes, not drainages so the downstream
				//   drainage link is incorrect.  This is handled above.  The user can also specify the RFLocn.
				//				select case (User[i-1].UsersType) !DSNode depends on the type of use
				//				case (InstreamFlowUseCode)
				//					DSNode=Drainage(j)%DSDrainage      ! for instream flow give drainageID d/s of source
				//				case default
				DSNode = j;
				//				end select
				CreateLink(NumLink, "RETURN", ReturnFlowLinkCode, InternalCode, USNode, DSNode, 0.0, m);
			}
		}
	} //making return flow links from sources to nodes

	//Measured Flow
	for (i = 1; i <= NumMeasuredFlowInfo; i++) { //one link for each measured_inflow-streamnode connection
		//find the MeasuredFlow node
		// find2()
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == MeasuredFlowNodeCode && Node[n].DrainageID == MeasuredFlowInfo[i-1].DrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		USNode = ifound[0]; //problem if nfound<>1
		delete [] ifound;

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
		DSNode = ifound[0]; //problem if nfound<>1
		delete [] ifound;
		CreateLink(NumLink, "MEASIN", MeasuredFlowLinkCode, ExternalCode, USNode,DSNode, 0.0, 0);
	}

	//one link for consumptive use, from each source to the sink node
	for (i = 1; i <= NumUser; i++) {
		// find2()
		found = false;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == UserNodeCode && Node[n].SelfID == i) {
				USNode = n+1;
				if (!found) {
					found = true;
				} else {
					cerr << "BuildLinkStructure(WASTE_): Duplicate link found" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		// find1()
		found = false;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == SinkNodeCode) {
				DSNode = n+1;
				if (!found) {
					found = true;
				} else {
					cerr << "BuildLinkStructure(WASTE_): Duplicate link found" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}

		CreateLink(NumLink, "WASTE_", SinkLinkCode, ExternalCode, USNode, DSNode, 0.0, 0);
	} //making links from sources to sink
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving BuildLinkStructure(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}


int CreateLink(int &NumLink, const string Title, const int LinkCode, const int IntExtCode,
	const int USNode, const int DSNode, const double Flow, const int ReturnFlowID)
{
	int k;
#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> CreateLink(" << Title << "), link #" << NumLink+1 << std::endl;
    }
	caller = "CreateLink";
#endif
	NumLink++;
	k = NumLink;
	Link(k-1).Title        = Title;
	Link(k-1).LinkCode     = LinkCode;
	Link(k-1).IntExtCode   = IntExtCode;
	Link(k-1).USNode       = USNode;
	Link(k-1).DSNode       = DSNode;
	Link(k-1).Flow         = Flow;
	Link(k-1).ReturnFlowID = ReturnFlowID; //only used by ReturnFlow links
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving CreateLink(" << Title << "), DSNode " << DSNode << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}
