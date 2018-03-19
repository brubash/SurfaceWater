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

// Local declarations
int CreateNode(int &NumNode, const string Title, const int NodeCode, const int IntExtCode,
    const double StoreMin, const double StoreMax, const double Store, const int DrainageID, const int SelfID);
int CreateUser(int &NumUser, const int UsersType, const int POU_ID, const double DemandVble, const double DemandRate,
	const int InYearDemandType, const int ReturnFlowID, const int SourceMixingID, const int NumSources,
	int SourceID[], int RightID[]);
int CreateSource(int &NumSource, const int SourcesType, const int SourceLocationID, const int RealSourceLocationID,
       const double PhysicalDailyMax, const double PhysicalAnnMax);
int CreateReturnFlow(int &NumReturnFlow, const int NumReturnFlows, const int ReturnFlowsUnits,
       const double ReturnFlowsAmt, const int ReturnFlowsType, const int ReturnFlowsLocn);


int BuildNodeStructure(const int NumDrainage, int &NumUser, const int NumReservoir, int &NumSource,
	int &NumReturnFlow, const int NumMeasuredFlowInfo, int &NumNode)
{
	int ReturnFlowsLocation = -999, ReservoirReleaseUseCode = -999;
	int i, j_ret, j_src, j_source, n, nfound, *ifound;
	double temp;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> BuildNodeStructure(" << ncalls << ")" << std::endl;
    }
	caller = "BuildNodeStructure";
#endif
	// Drainages
	NumNode = 0;
	for (i = 1; i <= NumDrainage; i++) { //Create 3 nodes for each drainage: one for the outlet, one for groundwater, one for runoff-production
		//Stream nodes (one node at the outlet of each model unit drainage)
		CreateNode(NumNode, "STREAM", StreamNodeCode, InternalCode, 0.0, 0.0, 0.0, i, i);
		//Groundwater nodes (one node per model unit drainage)
		CreateNode(NumNode, "GWATER", GroundwaterNodeCode, ExternalCode, -1.0e20, 0.0, 0.0, i, i);
     	   //Model Unit Drainage nodes (one node for the production of each local runoff)
		CreateNode(NumNode, "RUNOFF", DrainageNodeCode, ExternalCode, 0.0, 0.0, 0.0, i, i);
	}

	//Reservoirs
	for (i = 1; i <= NumReservoir; i++) {
		// three nodes per reservoir, a user node to take water, a reservoir node to hold it, and a user node to release water
		// (in-stream and off-stream are identical here, only differ in the links they set up)
		// in this loop we only create one node (reservoir itself), and do the
		// preparatory work of creating users, so that two other user nodes per reservoir
		// can be set up in the next loop. For each user node we will make, we will need to know
		// about the User, including their Source and ReturnFlow details.

		// Reservoir node which can store water
		// ReservoirID	DrainageID	InOffStream	RightID	StoreMax	StoreInitial	StoreMin	MaxInflow	MaxWithdrawal	MinEnvRelease	LossRate
		CreateNode(NumNode, "RSVOIR", ReservoirNodeCode, InternalCode, Reservoir[i-1].StoreMin, Reservoir[i-1].StoreMax,
                     Reservoir[i-1].StoreInitial, Reservoir[i-1].DrainageID,i);

		if (Reservoir[i-1].InOffStream == InStreamReservoirCode) {
                     //Preparation for first new user node
                     //Create an entry in the source table so we can get water for this reservoir
                     //SourceID	Type	SourceLocationID	PhysicalDailyMax	PhysicalAnnMax
			CreateSource(NumSource, StreamSourceCode, Reservoir[i-1].DrainageID,Reservoir[i-1].RealDrainageID,
                            Reservoir[i-1].MaxInflow, BigReal);
			j_src = NumSource;

                     //Create an entry in the return flow table, to send water from the user to the reservoir
                     //ReturnFlowID	NumReturnFlows	ReturnFlowUnits	ReturnFlowAmt   ReturnFlowType	ReturnFlowLocn
			CreateReturnFlow(NumReturnFlow, 1, FractionUnits, 1.0, ReservoirNodeCode, Reservoir[i-1].ReservoirID);
			j_ret = NumReturnFlow;

                     //Create User to take the water for reservoir
			temp = min(Reservoir[i-1].StoreMax - Reservoir[i-1].StoreMin, Reservoir[i-1].MaxInflow);
			CreateUser(NumUser, ReservoirFillUseCode, Reservoir[i-1].DrainageID, temp, 1.0, 0,
				ReturnFlow[j_ret-1].ReturnFlowID, 0, 1, &Source[j_src-1].SourceID, &Reservoir[i-1].RightID);
                     //END of Preparation for first new user node
		}

		//Preparation for second new user node
		//Create an entry in the source table so the "release" user can get their water from somewhere
		CreateSource(NumSource, ReservoirSourceCode, Reservoir[i-1].ReservoirID, Reservoir[i-1].ReservoirID, 1.0e20, 1.0e20);
		j_src = NumSource;

              //Create an entry in the return flow table, to send all water from the reservoir to the stream (possibly d/s)
		switch (Reservoir[i-1].InOffStream) {
			case InStreamReservoirCode:
				// find1()
				nfound = 0; //none found
				ifound = new int[NumDrainage];
				ifound[nfound] = 0;
				for (n = 0; n < NumDrainage; n++) {
					if (Drainage(n).DrainageID == Reservoir[i-1].DrainageID) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				j_source = ifound[0]; //problem if nfound<>1
				delete [] ifound;

				ReturnFlowsLocation = Drainage(j_source-1).DSDrainage;
				ReservoirReleaseUseCode = InStreamReservoirReleaseUseCode;
				break;
			case OffStreamReservoirCode:
				ReturnFlowsLocation = Reservoir[i-1].DrainageID;
				ReservoirReleaseUseCode = OffStreamReservoirReleaseUseCode;
				break;
			default: ;
		}
		CreateReturnFlow(NumReturnFlow, 1, FractionUnits, 1.0, StreamNodeCode, ReturnFlowsLocation);
			j_ret = NumReturnFlow;

		//Create User to take the overflow water and MinEnvRelease from reservoir
		//    CreateUser(NumUser,UsersType,POU_ID,DemandType,DemandVble,DemandRate,InYearDemandType, ...
		//                                    ReturnFlowID,SourceMixingID,NumSources,SourceID,RightID);
		temp = min(Reservoir[i-1].StoreMax-Reservoir[i-1].StoreMin, Reservoir[i-1].MinEnvRelease);
		CreateUser(NumUser, ReservoirReleaseUseCode, Reservoir[i-1].DrainageID, temp, 1.0, 0,
			ReturnFlow[j_ret-1].ReturnFlowID, 0, 1, &Source[j_src-1].SourceID, 0);
		//END of Preparation for second new user node

	} //of loop on Reservoirs

       //Users
	for (i = 1; i <= NumUser; i++) { //create one node for each user
		//one node for each user !!!!Node(k).UsersSourceID=(1:User(i).NumSources);
		CreateNode(NumNode, "USER__", UserNodeCode, InternalCode, 0.0, 0.0, 0.0, User[i-1].POU_ID, i);
		User[i-1].NodeNumber = NumNode; //this saves us having to find the node for this user
	}

       //Measured Flows
	for (i = 1; i <= NumMeasuredFlowInfo; i++) {  //one node per measured inflow, to bring in the required flow difference to get the inflow at the node right
		//one node for each Measured Inflow
		CreateNode(NumNode, "MEASIN", MeasuredFlowNodeCode, ExternalCode, 0.0, 0.0, 0.0, MeasuredFlowInfo[i-1].DrainageID, -1);
	}
	//Sink - one node for all outlets (could split later)
	CreateNode(NumNode, "SINK__", SinkNodeCode, ExternalCode, 0.0, 0.0, 0.0, -1, -1);
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving BuildNodeStructure(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

int CreateNode(int &NumNode, const string Title, const int NodeCode, const int IntExtCode,
    const double StoreMin, const double StoreMax, const double Store, const int DrainageID, const int SelfID)
{
	int k;

	NumNode++;
	k = NumNode;
	Node[k-1].Title      = Title;
	Node[k-1].Type       = NodeCode;
	Node[k-1].IntExt     = IntExtCode;
	Node[k-1].StoreMin   = StoreMin;
	Node[k-1].StoreMax   = StoreMax;
	Node[k-1].Store      = Store;
	Node[k-1].StoreOld   = Store;
	Node[k-1].DrainageID = DrainageID;
	Node[k-1].SelfID     = SelfID; //only the user nodes need this SelfID so far (it's an index back into the table that this type of node comes from)

	return 0;
}

int CreateUser(int &NumUser, const int UsersType, const int POU_ID, const double DemandVble, const double DemandRate,
	const int InYearDemandType, const int ReturnFlowID, const int SourceMixingID, const int NumSources,
	int SourceID[], int RightID[])
{
	int k, i;
	int idmax;

	NumUser++;
	k = NumUser;
	if (k == 1) {
		User[k-1].UserID = 1;
	} else {
		idmax = 0;
		for (i = 1; i <= k-1; i++) {
			idmax = max(idmax, User[i-1].UserID);
		}
		User[k-1].UserID = 1 + idmax;
	}
	User[k-1].UsersType        = UsersType;
	User[k-1].POU_ID           = POU_ID;
	User[k-1].DemandVble       = DemandVble;
	User[k-1].DemandRate       = DemandRate;
	User[k-1].InYearDemandType = InYearDemandType;
	User[k-1].ReturnFlowID     = ReturnFlowID;
	User[k-1].SourceMixingID   = SourceMixingID;
	User[k-1].NumSources       = NumSources;
	if (SourceID != 0) {
		for (i = 0; i < MaxNumSources; i++) {
			User[k-1].SourceID[i] = SourceID[i];
		}
	} else {
		for (i = 0; i < MaxNumSources; i++) {
			User[k-1].SourceID[i] = 0;
		}
	}
	if (RightID != 0) {
		for (i = 0; i < MaxNumSources; i++) {
              User[k-1].RightID[i]  = RightID[i];
		}
	} else {
		for (i = 0; i < MaxNumSources; i++) {
              User[k-1].RightID[i]  = 0;
		}
	}

	return 0;
}

int CreateSource(int &NumSource, const int SourcesType, const int SourceLocationID, const int RealSourceLocationID,
       const double PhysicalDailyMax, const double PhysicalAnnMax)
{
	int k, i;
	int idmax, ids;

	NumSource++;
	k = NumSource;

	if (k <= 1) {
		ids = 1;
	} else {
		idmax = 0;
		for (i = 1; i <= k-1; i++) {
			idmax = max(idmax, Source[i-1].SourceID);
		}
		ids = 1 + idmax;
	}

	Source[k-1].SourceID             = ids;
	Source[k-1].Type                 = SourcesType;
	Source[k-1].SourceLocationID     = SourceLocationID;
	Source[k-1].RealSourceLocationID = RealSourceLocationID;
	Source[k-1].PhysicalDailyMax     = PhysicalDailyMax;
	Source[k-1].PhysicalAnnMax       = PhysicalAnnMax;

	return 0;
}

int CreateReturnFlow(int &NumReturnFlow, const int NumReturnFlows, const int ReturnFlowsUnits,
	const double ReturnFlowsAmt, const int ReturnFlowsType, const int ReturnFlowsLocn)
{
	int k, i;
	int idmax;

	NumReturnFlow++;
	k = NumReturnFlow;
	if (k == 1) {
		ReturnFlow[k-1].ReturnFlowID = 1;
	} else {
		idmax = 0;
		for (i = 1; i <= k-1; i++) {
			idmax = max(idmax, ReturnFlow[i-1].ReturnFlowID);
		}
		ReturnFlow[k-1].ReturnFlowID = 1 + idmax;
       }
	ReturnFlow[k-1].NumReturnFlows     = NumReturnFlows;
	ReturnFlow[k-1].ReturnFlowsUnits   = ReturnFlowsUnits;
	ReturnFlow[k-1].ReturnFlowsAmt[0]  = ReturnFlowsAmt;
	ReturnFlow[k-1].ReturnFlowsType[0] = ReturnFlowsType;
	ReturnFlow[k-1].ReturnFlowsLocn[0] = ReturnFlowsLocn;

	return 0;
}

