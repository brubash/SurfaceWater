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
using namespace std;
using namespace Eigen;

int TrialFlow(const double DemandByUserOnSource, const int i, const int j, const int j_source,
	const int NumRights, double &Qtry);

int AllocateWaterToUsers(const int Timestep, const int NumNode, const int NumLink, const int NumUser, const int NumReservoir,
	const int NumSource, const int NumRights, const int NumSourceMixing, int *DrainageOrder, const int NumDrainage,
	const int NumReturnFlow, const int NumUserSource, ArrayXd &volume_irrig_sup, double *groundwater_to_take, ArrayXd &DrainageOutFlow)
{
	int i_count, i, j, ii, k, n, j_source, ilink, inode, j_srcmx, i_count_max, iFeasible, nfound, *ifound;
	double DemandByUserOnSource, Qtry0, Qtry, Capacity, AllocatedFlow;
#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> AllocateWaterToUsers(" << ncalls << ")" << std::endl;
    }
	caller = "AllocateWaterToUsers";
#endif
	//a user may need a right to take water from a particular source
	//so we process the user-source combinations in order of priority
	//this may mean that for example we allocate some water to one of a user's sources today,
	//but then the next water from that source gets allocated to another user. We might
	//then later come back to allocate water to other sources which the same user has.
	//It is the user-source combinations that are ordered
	//(either by priority date or by position in the river system), not the users.
	for (i = 0; i < NumUser; i++) {
		User[i].Deficit[Timestep-1] = User[i].DemandToday;
		User[i].Demand[Timestep-1]  = User[i].DemandToday;
		User[i].Withdrawal[Timestep-1] = 0.0;
	}
	for (ii = 1; ii <= NumUserSource; ii++) {
		k = UserSourceOrder[ii-1];
		i = UserSourceTable[k-1].UserID;
		if (User[i-1].UsersType == InStreamReservoirReleaseUseCode || User[i-1].UsersType == OffStreamReservoirReleaseUseCode) {
			//  Reservoir release is spill due to fill and spill operation rule.
			//  The demand associated with a reservoir release user is the storage in excess of capacity
			ilink = User[i-1].LinkSourceToUser[0];
			inode = Link(ilink-1).USNode;
			User[i-1].DemandToday = max(Node[inode-1].Store - Node[inode-1].StoreMax, User[i-1].Demand[Timestep-1]);
			User[i-1].Deficit[Timestep-1] = User[i-1].DemandToday;
			User[i-1].Demand[Timestep-1] = User[i-1].DemandToday;
			User[i-1].Withdrawal[Timestep-1] = 0.0;
		} else if (User[i-1].UsersType == ReservoirFillUseCode) {
			// Reservoir fill user requests difference between capacity and storage
			ilink = User[i-1].LinkUserToReturnflow[0];
			inode = Link(ilink-1).DSNode;
			User[i-1].DemandToday = max(Node[inode-1].StoreMax - Node[inode-1].Store, 0.0);
			User[i-1].Deficit[Timestep-1] = User[i-1].DemandToday;
			User[i-1].Demand[Timestep-1] = User[i-1].DemandToday;
			User[i-1].Withdrawal[Timestep-1] = 0.0;
		}
		if (User[i-1].DemandToday > 0) {
			j = UserSourceTable[k-1].SourceCounter;
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
			j_source = ifound[0];	//problem if nfound<>1
			delete [] ifound;
			if (User[i-1].SourceMixingID == 0 || User[i-1].NumSources == 1) {
				DemandByUserOnSource = User[i-1].DemandToday;
			} else {
				j_srcmx = User[i-1].j_srcmx[j-1];
				if (SourceMixing[j_srcmx-1].Units == FractionUnits) {
					DemandByUserOnSource = User[i-1].DemandToday * SourceMixing[j_srcmx-1].Amount;
				} else {
					DemandByUserOnSource = min( User[i-1].DemandToday, SourceMixing[j_srcmx-1].Amount );
				} //SourceMixing(j_srcmx)%Units == FractionUnits
			} //User[i-1].SourceMixingID == 0 || User[i-1].NumSources == 1

			TrialFlow(DemandByUserOnSource, i, j, j_source, NumRights, Qtry0);
			Qtry = Qtry0;
			Capacity = Qtry;  // 10/5/05 Ross and Dave think that Capacity does not need to be initialized here
			AllocatedFlow = 0.0;
			for (n = 0; n < NumLink; n++) {
				LinkSave(n).Title        = Link(n).Title;
				LinkSave(n).LinkCode     = Link(n).LinkCode;
				LinkSave(n).IntExtCode   = Link(n).IntExtCode;
				LinkSave(n).USNode       = Link(n).USNode;
				LinkSave(n).DSNode       = Link(n).DSNode;
				LinkSave(n).Flow         = Link(n).Flow;
				LinkSave(n).ReturnFlowID = Link(n).ReturnFlowID;
			}
			for (n = 0; n < NumNode; n++) {
				NodeSave[n].Title      = Node[n].Title;
				NodeSave[n].Type       = Node[n].Type;
				NodeSave[n].IntExt     = Node[n].IntExt;
				NodeSave[n].StoreMin   = Node[n].StoreMin;
				NodeSave[n].StoreMax   = Node[n].StoreMax;
				NodeSave[n].Store      = Node[n].Store;
				NodeSave[n].StoreOld   = Node[n].StoreOld;
				NodeSave[n].DrainageID = Node[n].DrainageID;
				NodeSave[n].SelfID     = Node[n].SelfID;
			}
//double T0 = (double)clock()/(double)CLOCKS_PER_SEC;cout << "Qtry " << fixed << setw(8) << setprecision(2) << Qtry << '\n';
			i_count = 1;
			i_count_max = 100;
			while (Qtry > 0.01*Qtry0 && AllocatedFlow == 0.0 && i_count <= i_count_max) {

				//give Qtry to the user node and follow it through
				PropagateWaterViaUser(i, j, Qtry, NumNode, NumLink, NumUser, NumReservoir, NumSource,
					DrainageOrder, NumDrainage, NumReturnFlow, iFeasible, Capacity);
				//take Qtry from the source, and follow that through
				//            select case (Source(j_source)%Type)
				//                case (StreamSourceCode)
				//                    call PropagateTakeFromStream(i,j,NumNode, &
				//						NumLink,DrainageOrder,NumDrainage,iFeasible,Capacity);
				//                case (ReservoirSourceCode,GroundwaterSourceCode)
				//                    call PropagateTakeFromStore(Qtry,i,j,NumNode,NumLink, &
				//                        NumUser,NumSource,iFeasible,Capacity);
				//            end select

				if (iFeasible != 0) { //true: this flow is feasible
					AllocatedFlow = Qtry;
				} else {
					//Qtry=Qtry*.99;
					Qtry = max(0.0, min(Qtry+Capacity, Qtry0*(1.0 - (double)(i_count)/(double)(i_count_max))));
					for (n = 0; n < NumLink; n++) {
						Link(n).Title        = LinkSave(n).Title;
						Link(n).LinkCode     = LinkSave(n).LinkCode;
						Link(n).IntExtCode   = LinkSave(n).IntExtCode;
						Link(n).USNode       = LinkSave(n).USNode;
						Link(n).DSNode       = LinkSave(n).DSNode;
						Link(n).Flow         = LinkSave(n).Flow;
						Link(n).ReturnFlowID = LinkSave(n).ReturnFlowID;
					}
					for (n = 0; n < NumNode; n++) {
						Node[n].Title      = NodeSave[n].Title;
						Node[n].Type       = NodeSave[n].Type;
						Node[n].IntExt     = NodeSave[n].IntExt;
						Node[n].StoreMin   = NodeSave[n].StoreMin;
						Node[n].StoreMax   = NodeSave[n].StoreMax;
						Node[n].Store      = NodeSave[n].Store;
						Node[n].StoreOld   = NodeSave[n].StoreOld;
						Node[n].DrainageID = NodeSave[n].DrainageID;
						Node[n].SelfID     = NodeSave[n].SelfID;
					}
					i_count++;
				} //iFeasible == 1
			} // while.  Here we are done iterating.  Assignments below are permanent
/*double T1 = (double)clock()/(double)CLOCKS_PER_SEC;
double deltaT = T1 - T0;
if (deltaT > 1.0e-3) {
    cout << scientific << setw(6) << setprecision(1) << deltaT << " seconds for AllocateWaterToUsers()\n";
} */
			User[i-1].Deficit[Timestep-1] = max(0.0, User[i-1].Deficit[Timestep-1] - Qtry);
			User[i-1].Withdrawal[Timestep-1] = User[i-1].Withdrawal[Timestep-1] + Qtry;
			User[i-1].VolumeToDateSource[j-1] = Qtry + User[i-1].VolumeToDateSource[j-1];
			BalanceFlowsAtStreamNodes(NumNode, NumLink, DrainageOrder, NumDrainage, DrainageOutFlow);
			//for passing back to Topnet  ! DGT 8/16/05.  Type -1 left for compatibility with old inputs
			//  type 7 is irrigation with demand calculated based on soil moisture
			//  type 8 is irrigation with demand specified by user files
			//		if (User[i-1].UsersType == -1 .or. User[i-1].UsersType == 7 .or. User[i-1].UsersType == 8) {
			if (User[i-1].UsersType == SoilMoistureIrrigationUseCode || User[i-1].UsersType == FixedDemandIrrigationUseCode) {
				volume_irrig_sup(User[i-1].POU_ID-1) += Qtry;
			} //User[i-1].UsersType == -1
			if (Source[j_source-1].Type == GroundwaterSourceCode) {
				groundwater_to_take[User[i-1].POU_ID-1] += Qtry;
			} //Source(j_source)%Type == GroundwaterSourceCode
		} // User[i-1].DemandToday.gt.0
	} //ii
	// 4. For each user USER__-ID, in descending order of priority
	// For each source of that user, in no particular order
	// Note the SOURCENODE-ID of the source node for this user-source combination, and find the link LINK-ID with LinkType = 4  ABSTRCT and LinkUSNode=SOURCENODE-ID and LinkDSNode=USER__-ID
	// Note the RETURNNODE-ID of the return node for this user-source combination, and find the link RETURNLINK-ID with LinkType = 5  RETURN and LinkUSNode=USER__-ID and LinkDSNode=RETURNNODE-ID
	// Calculate Qtry = min(demand - Flow in link LINK-ID, daily_max_right, annual_max_take-take_to_date)
	// Determine the amount of returnflow to RETURNNODE-ID for Qtry
	// Determine whether it is feasible to allocate Qtry to user from source: (if Qtry is small skip to 2.4.1.6)
	// increase Flow in the link LINK-ID by Qtry
	// increase Flow in the link RETURNLINK-ID by returnflow
	// increase Flow in the link SINKNODE-ID by Qtry - returnflow
	// if the source node has max-min>0 then
	// if the node storage stays in its bounds, reduce the flow stored at node SOURCENODE-ID by Qtry and this allocation is feasible, else it is not feasible
	// elseif the sourcenode has no storage then recalculate Flow in all links of LinkType = 2  UNALLC, by adjusting the unallocated flow at every node THISNODE-ID so that there is water balance at the node. Inflow is Sum of Flow in all links with LinkDSNode=THISNODE-ID.   AllocatedOutflow is Sum of Flow in all links with LinkUSNode=THISNODE-ID and LinkType<>2  UNALLC.  Flow in LinkType=2  UNALLC is Inflow- AllocatedOutflow. Process the nodes from upstream to downstream. If all unallocated flows are non-negative then this allocation of Qtry is feasible otherwise it is not feasible.
	// If the attempted allocation of Qtry was not feasible then subtract recently-added Flow in links LINK-ID RETURNLINK-ID and SINKNODE-ID, else if it was feasible then update the users stats for the year to date
	// if Qtry is not too small, reduce Qtry and return to step 2.4.1.5, otherwise set Qtry=0
	// Weve now allocated as much water as possible to this user from this source
	// Go on to the next source for this user
	// Go on to the next user
	// Might need to embed another loop in here so that we process users in priority order 12123123412345, rather than 12345
	// For MI=1 to n; For J=1 to M, do steps 2.4.1 and 2.4.2 with user J; next J; next M
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving AllocateWaterToUsers(" << ncalls << ")\n" << endl;
    }
    ncalls++;
#endif

	return 0;
}

int TrialFlow(const double DemandByUserOnSource, const int i, const int j, const int j_source,
	const int NumRights, double &Qtry)
{
	double Qmin_inst, Qmin_ann;
	int j_rights = -1;
	bool found;
	//find a trial flow which tries to meet demand subject to Daily and Annual Water Right Limits
	int ianypositive = 0, ii, n;//, nfound, *ifound;
	double xLegalDailyMax, xLegalAnnMax;
#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> TrialFlow(" << ncalls << ")" << std::endl;
    }
	caller = "TrialFlow";
#endif
	for (ii = 1; ii <= User[i-1].NumSources; ii++) {
		if (User[i-1].RightID[ii-1] > 0) {
			ianypositive = 1;
		}
	}

	if (NumRights > 0 && ianypositive == 1) {
		// find1()
		found = false;
		for (n = 0; n < NumRights; n++) {
			if (Rights[n].RightID == User[i-1].RightID[j-1]) {
				j_rights = n+1;
				if (!found) {
					found = true;
				} else {
					cerr << "TrialFlow(): Duplicate link found" << endl;
					exit(EXIT_FAILURE);
				}
			}
		}
		if (!found) {
			cerr << "TrialFlow() find1() failure\n";
			exit(EXIT_FAILURE);
		} else {
			// reset
			found = false;
		}

		xLegalDailyMax = Rights[j_rights-1].LegalDailyMax;
		xLegalAnnMax   = Rights[j_rights-1].LegalAnnMax;
	} else {
		xLegalDailyMax = 1.0e20;
		xLegalAnnMax   = 1.0e20;
	}
	Qmin_inst = min(xLegalDailyMax, Source[j_source-1].PhysicalDailyMax);
	Qmin_ann  = min(xLegalAnnMax, Source[j_source-1].PhysicalAnnMax) - User[i-1].VolumeToDateSource[j-1];
	Qtry      = min(Qmin_inst, min(Qmin_ann, DemandByUserOnSource));
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving TrialFlow(" << ncalls << ")" << endl;
    }
    ncalls++;
#endif

	return 0;
}
