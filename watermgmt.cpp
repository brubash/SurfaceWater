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
#include <iomanip>

using namespace std;
using namespace constant_definitions;
using namespace input_structures;
using namespace other_structures;
using namespace Eigen;

int *DrainageOrder;
//double *DrainageOutFlow;

int timecalcs(const int Timestep, const int dt, const long int i8startsecs,
	int &doy, int &ThisMonth, int &ThisDay, int &yyyymmdd, int &hhmmss);

int watermgmt(const int StartDateTopnet, int &StartHourTopnet, const int Timestep, const int NSteps,
	ArrayXd &RunoffTopnet, ArrayXd &BaseflowTopnet, const ArrayXd &ArtDrainageTopnet,
	const ArrayXd &vol_irrig_demand, const int maxSlp, const double *evaporation, const double *precipitation,
	ArrayXd &volume_irrig_sup, double *groundwater_to_take)
{
	int i, j, k, m, n, ii, jj, kk, nfound, *ifound;
	int isink, j_source, iSrcLocnID, i_node, j_sink, j_r, j_drainage;
	static int NumUserSource, NumUserSourceReturn, NumReturnFlows;
	//save !so we remember from one Timestep to the next without passing stuff
	double scalefactor;
	static ArrayXd DrainageOutFlow;
	static string dirname;

	static long int i8startsecs;
	static int doy, ThisMonth, ThisDay, dt, yyyymmdd, hhmmss;
	int *rfType, *rfLocn;
	static int nReturnFlows_to_GW = 0;
	static int *ReturnFlows_to_GW_LinkID, *ReturnFlows_to_GW_DrainageID;
	int MaxNodes, MaxLinks, RFNodeID;
	static double t1, t2, t3;
	static int NumLink, NumNode;
	static bool allocateFlag[2] = {false,false};
	bool found;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> watermgmt(" << ncalls << ")" << std::endl;
    }
	caller = "watermgmt";
#endif
	if (Timestep == 0) {
		// initialize
		if (NSteps > MaxNtimeSteps) {  // Checking for too many time steps
			cerr << NSteps << " time steps is more than the maximum allowed: " << MaxNtimeSteps << '\n';
			cerr << "Have a programmer adjust the dimensions in types.f90\n";
			cerr << "Exiting\n";
			exit(EXIT_FAILURE);
		}

		ccompr2::t0 = (double)clock()/(double)CLOCKS_PER_SEC;
		dirname = "."; //current directory
		//dirname = "bethexample";
		//dirname = "daveexample";
		//dirname = "nooksack1";
		dt = 86400;

		read_inputs(dirname, dt, StartDateTopnet, NumDrainage, NumStreamNode, NumMeasuredFlowInfo, NumMeasuredFlowData,
			NumReservoir, NumUser, NumSource, NumRights, NumSourceMixing, NumSeasonsDefn, NumReturnFlow, NumMonthlyDemand,
			NumRunoff, NumBaseflow, NumWWTP);
		RunControl.NumTimesteps = NSteps;

		//Build Node Structure
		MaxNodes = 3*NumDrainage + NumReservoir + (NumUser + 2*NumReservoir) + NumMeasuredFlowData + 1;
		Node = new NodeType[MaxNodes];
		BuildNodeStructure(NumDrainage, NumUser, NumReservoir, NumSource, NumReturnFlow, NumMeasuredFlowInfo, NumNode);
		NodeSave = new NodeType[MaxNodes];
		//NodeSave = Node;
		for (i = 0; i < MaxNodes; i++) {
			NodeSave[i].Title      = Node[i].Title;
			NodeSave[i].Type       = Node[i].Type;
			NodeSave[i].IntExt     = Node[i].IntExt;
			NodeSave[i].StoreMin   = Node[i].StoreMin;
			NodeSave[i].StoreMax   = Node[i].StoreMax;
			NodeSave[i].Store      = Node[i].Store;
			NodeSave[i].StoreOld   = Node[i].StoreOld;
			NodeSave[i].DrainageID = Node[i].DrainageID;
			NodeSave[i].SelfID     = Node[i].SelfID;
		}
		delete [] Node;
		Node = new NodeType[NumNode];
		for (i = 0; i < NumNode; i++) {
			Node[i].Title      = NodeSave[i].Title;
			Node[i].Type       = NodeSave[i].Type;
			Node[i].IntExt     = NodeSave[i].IntExt;
			Node[i].StoreMin   = NodeSave[i].StoreMin;
			Node[i].StoreMax   = NodeSave[i].StoreMax;
			Node[i].Store      = NodeSave[i].Store;
			Node[i].StoreOld   = NodeSave[i].StoreOld;
			Node[i].DrainageID = NodeSave[i].DrainageID;
			Node[i].SelfID     = NodeSave[i].SelfID;
		}
		delete [] NodeSave;
		NodeSave = new NodeType[NumNode];

		//Build Link Structure
		MaxLinks = 3*NumDrainage + 5*(NumUser + 2*NumReservoir) + NumMeasuredFlowData;
		Link.resize(MaxLinks);  // needed for the next step
		BuildLinkStructure(NumDrainage, NumUser, NumSource,	NumReturnFlow, NumReservoir, NumMeasuredFlowInfo, NumNode, NumLink);

		LinkSave.resize(MaxLinks);
		//LinkSave = Link;
		for (i = 0; i < MaxLinks; i++) {
			LinkSave(i).Title        = Link(i).Title;
			LinkSave(i).LinkCode     = Link(i).LinkCode;
			LinkSave(i).IntExtCode   = Link(i).IntExtCode;
			LinkSave(i).USNode       = Link(i).USNode;
			LinkSave(i).DSNode       = Link(i).DSNode;
			LinkSave(i).Flow         = Link(i).Flow;
			LinkSave(i).ReturnFlowID = Link(i).ReturnFlowID;
		}

		Link.resize(NumLink);
		for (i = 0; i < NumLink; i++) {
			Link(i).Title        = LinkSave(i).Title;
			Link(i).LinkCode     = LinkSave(i).LinkCode;
			Link(i).IntExtCode   = LinkSave(i).IntExtCode;
			Link(i).USNode       = LinkSave(i).USNode;
			Link(i).DSNode       = LinkSave(i).DSNode;
			Link(i).Flow         = LinkSave(i).Flow;
			Link(i).ReturnFlowID = LinkSave(i).ReturnFlowID;
		}
		LinkSave.resize(NumLink);

		//Build Drainage Ordering
		DrainageOrder = new int[NumDrainage];
		BuildDrainageOrder(NumDrainage, DrainageOrder);

		//Assign Priority To User-Source Combinations
		UserSourceOrder = new int[5*NumUser];
		UserSourceTable = new UserSourceTableType[5*NumUser];
		AssignPriorityOrder(NumUser, NumSource, NumRights, NumReservoir, RunControl.AllocationMode,
			DrainageOrder, NumDrainage, NumUserSource, NumUserSourceReturn, NumReturnFlow);
			// DGT 9/4/05 Added NumReturnFlow so that where it is used in this subroutine it is not undefined

		//Initialise Volume taken by each user from each source
		for (i = 1; i <= NumUser; i++) {
			for (j = 1; j <= User[i-1].NumSources; j++) {
				User[i-1].VolumeToDateSource[j-1] = 0.0;
			}
		}

		td8micsec(StartDateTopnet, StartHourTopnet, i8startsecs); //seconds to start time since start of 1 Jan 1940;
		StaticOutput.DrainageInfo     = new DrainageInfoType[NumUserSourceReturn];
		StaticOutput.DrainageInfoSize = NumUserSourceReturn;

		StaticOutput.StreamFlowLinks     = new StreamFlowLinksType[NumDrainage];
		StaticOutput.StreamFlowLinksSize = NumDrainage;

		StaticOutput.DrainageID     = new DrainageIDType[NumDrainage];
		StaticOutput.DrainageIDSize = NumDrainage;

		StaticOutput.UserFlowLinks = 0;

		Initialise_Output_Tables(NumDrainage, NumNode, NumStreamNode, NumLink, NumUser,
			RunControl.NumTimesteps, NumReservoir, NumUserSourceReturn, NumReturnFlow);

		//Set up some indices to save time on searching
		for (ii = 1; ii <= NumDrainage; ii++) {
			i = DrainageOrder[ii-1];

			// find2() 6/3 secs
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == StreamNodeCode && Node[n].DrainageID == i) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Drainage(i-1).isn = ifound[0];	//problem if nfound<>1
			delete [] ifound;

			// find1()  4/3 secs
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).DSNode == Drainage(i-1).isn) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Drainage(i-1).n_in = nfound;
			for (n = 0; n < nfound; n++) {
				Drainage(i-1).ifound_in(n) = ifound[n];
			}
			delete [] ifound;

			// find2a() A & (~B)  7/3 secs
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).USNode == Drainage(i-1).isn && Link(n).LinkCode != UnallocatedLinkCode) { 	//A&(~B)
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}

			Drainage(i-1).n_taken = nfound;
			for (n = 0; n < nfound; n++) {
				Drainage(i-1).ifound_taken(n) = ifound[n];
			}
			delete [] ifound;

			// find2() 7/3 secs
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).USNode == Drainage(i-1).isn && Link(n).LinkCode == UnallocatedLinkCode) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Drainage(i-1).j_unallocated = ifound[0]; //problem if nfound<>1
			delete [] ifound;
		}

		for (i = 1; i <= NumReservoir; i++) {
			// find2()
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == ReservoirNodeCode && Node[n].SelfID == i) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Reservoir[i-1].isn = ifound[0]; //problem if nfound<>1
			delete [] ifound;

			// find1()
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).DSNode == Reservoir[i-1].isn) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Reservoir[i-1].n_in = nfound;
			for (n = 0; n < nfound; n++) {
				Reservoir[i-1].ifound_in(n) = ifound[n];
			}
			delete [] ifound;

			// find1()
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).USNode == Reservoir[i-1].isn) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			Reservoir[i-1].n_taken = nfound;
			for (n = 0; n < nfound; n++) {
				Reservoir[i-1].ifound_taken(n) = ifound[n];
			}
			delete [] ifound;
		}

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
		isink = ifound[0];
		delete [] ifound;

		for (ii = 1; ii <= NumUserSource; ii++) {
			kk = UserSourceOrder[ii-1];
			i = UserSourceTable[kk-1].UserID;
			j = UserSourceTable[kk-1].SourceCounter;

			// find1() User(i)%source_ind(1)=j_source
			nfound = 0; //none found
			ifound = new int[NumSource];
			ifound[nfound] = 0;
			for (n = 0; n < NumSource; n++) {
				if (Source[n].SourceID == User[i-1].SourceID[j-1]) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			j_source = ifound[0];  //problem if nfound<>1
			delete [] ifound;

			// find2()
			nfound = 0; //none found
			ifound = new int[NumSourceMixing];
			ifound[nfound] = 0;
			for (n = 0; n < NumSourceMixing; n++) {
				if (SourceMixing[n].SourceMixingID == User[i-1].SourceMixingID && SourceMixing[n].UsersSourceNum == j) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			User[i-1].j_srcmx[j-1] = ifound[0]; //problem if nfound<>1
			delete [] ifound;

			//k is node of this user
			k = User[i-1].NodeNumber;
			iSrcLocnID = Source[j_source-1].SourceLocationID;
			if (iSrcLocnID == 0)
				iSrcLocnID = User[i-1].POU_ID;

			//jj is link from source to user
			// find2()
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == Source[j_source-1].Type && Node[n].SelfID == iSrcLocnID) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			i_node = ifound[0]; //problem if nfound<>1
			delete [] ifound;

			// find3()
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).LinkCode == UserAbstractionLinkCode && Link(n).USNode == i_node  && Link(n).DSNode == k) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			jj = ifound[0]; //problem if nfound<>1
			delete [] ifound;
			User[i-1].LinkSourceToUser[j-1] = jj;

			// find2()
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).LinkCode == SinkLinkCode && Link(n).USNode == k) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			j_sink = ifound[0]; //problem if nfound<>1
			delete [] ifound;
			User[i-1].LinkSourceToSink = j_sink;

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
				j_r = ifound[0];  //problem if nfound<>1
				delete [] ifound;

				NumReturnFlows = ReturnFlow[j_r-1].NumReturnFlows;
				rfType = new int[NumReturnFlows];
				rfLocn = new int[NumReturnFlows];
				for (m = 1; m <= NumReturnFlows; m++) {
					rfType[m-1] = ReturnFlow[j_r-1].ReturnFlowsType[m-1];
					if (rfType[m-1] == 0)
						rfType[m-1] = StreamNodeCode;
					rfLocn[m-1] = ReturnFlow[j_r-1].ReturnFlowsLocn[m-1];
					if (rfLocn[m-1] == 0) {
						if (User[i-1].UsersType == InstreamFlowUseCode) {
							rfLocn[m-1] = Drainage(User[i-1].POU_ID-1).DSDrainage;   // This is the downstream drainge
						} else {
							rfLocn[m-1] = User[i-1].POU_ID;
						}
					}
					if (rfType[m-1] == 1 && rfLocn[m-1]==-1) {

						// find1()
						nfound = 0; //none found
						ifound = new int[NumDrainage];
						ifound[nfound] = 0;
						for (n = 0; n < NumDrainage; n++) {
							if (Drainage(n).DrainageID == User[i-1].POU_ID) {
								nfound++;
								ifound[nfound-1] = n+1;
							}
						}
						j_drainage = ifound[0];  //problem if nfound<>1
						delete [] ifound;
						rfLocn[m-1] = Drainage(j_drainage-1).DSDrainage;
					}
				}
			} else if (User[i-1].ReturnFlowID < 0) {
				NumReturnFlows = 0;
			} else {
				NumReturnFlows = 1;
				rfType = new int[NumReturnFlows];
				rfLocn = new int[NumReturnFlows];
				NumReturnFlows = 1;
				rfType[0] = 1;
				rfLocn[0] = 1;
			}
			for (m = 1; m <= NumReturnFlows; m++) { //length(User[i-1].ReturnFlowFrac) //return flow(s) from user to return flow nodes
				if (rfLocn[m-1] > 0) {
					// find2()
					nfound = 0; //none found
					ifound = new int[NumNode];
					ifound[nfound] = 0;
					for (n = 0; n < NumNode; n++) {
						if (Node[n].Type == rfType[m-1] && Node[n].SelfID == rfLocn[m-1]) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					RFNodeID = ifound[0]; //problem if nfound<>1
					delete [] ifound;
				} else {
					RFNodeID = isink;
				}
				//j_return= link from user node to returnflow node
				// find3()
				found = false;
				for (n = 0; n < NumLink; n++) {
					if (Link(n).LinkCode == ReturnFlowLinkCode && Link(n).USNode == k  && Link(n).DSNode == RFNodeID) {
						User[i-1].LinkUserToReturnflow[m-1] = n+1;
						if (!found) {
							found = true;
						} else {
							cerr << "watermgmt(): Duplicate link found" << endl;
							exit(EXIT_FAILURE);
						}
					}
				}
				if (!found) {
					cerr << "watermgmt() find3() failure\n";
					exit(EXIT_FAILURE);
				} else {
					// reset
					found = false;
				}


			}

		}

		// find1() identify the links which provide return flows
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).LinkCode == ReturnFlowLinkCode) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}

		nReturnFlows_to_GW = 0;
		ReturnFlows_to_GW_DrainageID = new int[nfound];
		ReturnFlows_to_GW_LinkID = new int[nfound-1];
		for (i = 1; i <= nfound; i++) {
			if (Node[Link(ifound[i-1]-1).DSNode-1].Type == GroundwaterNodeCode) { // if this returnflow is going to a G/W node, we need to know about it
				nReturnFlows_to_GW++;
				ReturnFlows_to_GW_LinkID[nReturnFlows_to_GW-1] = ifound[i-1];
				ReturnFlows_to_GW_DrainageID[nReturnFlows_to_GW-1] = Node[Link(ifound[i-1]-1).DSNode-1].DrainageID;
			}
		}
		delete [] ifound;

		t1 = (double)clock()/(double)CLOCKS_PER_SEC;
		cout << t1 - ccompr2::t0 << " seconds to read input files\n";
		//  Stuff for writing files as we go
		scalefactor = 1.0;
		Write_OutputLine_Eigen(oFile[8], "results/Artificial_Drainage_cms.txt",       Timestep, ArtDrainageTopnet, NumDrainage, scalefactor);
		Write_OutputLine(oFile[9], "results/Precipitation_mm.txt",                    Timestep, precipitation, NumDrainage, scalefactor);
		Write_OutputLine(oFile[10], "results/Evaporation_mm.txt",                     Timestep, evaporation, NumDrainage, scalefactor);
		Write_OutputLine_Eigen(oFile[11], "results/Baseflow_cms.txt",                 Timestep, BaseflowTopnet, NumDrainage, scalefactor);
		Write_OutputLine_Eigen(oFile[12], "results/TotalRunoff_noWithdrawal_cms.txt", Timestep, RunoffTopnet, NumDrainage, scalefactor);
		//  Added this function to output Local Contributions to each StreamNode
		Write_OutputLocalContributions(oFile[13], oFile[14], NumStreamNode, NumDrainage, BaseflowTopnet, RunoffTopnet, Timestep, scalefactor);

		// END OF INITIALISE

	} else if (Timestep > 0) {

		// do this Timestep
		timecalcs(Timestep, dt, i8startsecs, doy, ThisMonth, ThisDay, yyyymmdd, hhmmss);

		if (ThisMonth == RunControl.StartOfWaterYearmm && ThisDay == RunControl.StartOfWaterYeardd) {
			for (i = 1; i <= NumUser; i++) {
				for (j = 1; j <= User[i-1].NumSources; j++) {
					User[i-1].VolumeToDateSource[j-1] = 0.0;
				}
			}
		}

		for (j = 1; j <= NumLink; j++) {
			Link(j-1).Flow = 0.0;
		}

		//use the runoff values to assign flows to links
		if (!allocateFlag[1]) {
			DrainageOutFlow.resize(NumDrainage);
			allocateFlag[1] = true;
		}
		//DrainageOutFlow = new double[NumDrainage];
		if (!allocateFlag[0]) {
			Runoff          = new RunoffType[NSteps];
			Baseflow        = new RunoffType[NSteps];
			allocateFlag[0]    = true;
		}

		Runoff[Timestep-1].Timestep = Timestep;
		Runoff[Timestep-1].Date = yyyymmdd;
		Runoff[Timestep-1].Hour = hhmmss;
		for (n = 0; n < NumDrainage; n++) {
			Runoff[Timestep-1].Rate[n] = RunoffTopnet(n); //m3/Timestep
		}
		Baseflow[Timestep-1].Timestep = Timestep;
		Baseflow[Timestep-1].Date = yyyymmdd;
		Baseflow[Timestep-1].Hour = hhmmss;
		for (n = 0; n < NumDrainage; n++) {
			Baseflow[Timestep-1].Rate[n] = BaseflowTopnet(n); //m3/Timestep
		}
		AssignDrainageFlows(Timestep, NumDrainage, NumNode, NumLink, DrainageOrder, NumRunoff, NumBaseflow, DrainageOutFlow);
		//compare the measured runoffs at boundary condition sites to the flow
		//in links there, and adjust the runoff and baseflows so the modelled
		//runoffs match the measured flows, and then recalculate the flows in the links
		ImposeMeasuredFlows(Timestep, NumNode, NumLink, NumRunoff, NumBaseflow, NumDrainage,
			NumMeasuredFlowInfo, NumMeasuredFlowData, DrainageOrder, DrainageOutFlow);

		CalculateDemand(ThisMonth, doy, vol_irrig_demand, NumUser, NumMonthlyDemand);
		for (i = 0; i < maxSlp; i++) {
			volume_irrig_sup(i)    = 0.0;
			groundwater_to_take[i] = 0.0;
		}
		for (n = 0; n < NumNode; n++) {
			Node[n].StoreOld = Node[n].Store;
		}
		if (RunControl.AllocationMode != NoAllocationCode) {
			AllocateWaterToUsers(Timestep, NumNode, NumLink, NumUser, NumReservoir, NumSource, NumRights, NumSourceMixing,
				DrainageOrder, NumDrainage, NumReturnFlow, NumUserSource, volume_irrig_sup, groundwater_to_take, DrainageOutFlow);
			for (i = 1; i <= nReturnFlows_to_GW; i++) { //if any of the return flows deliver water to GW, tell Topnet about it.
				groundwater_to_take[ReturnFlows_to_GW_DrainageID[i-1]-1] -= Link(ReturnFlows_to_GW_LinkID[i-1]-1).Flow;
					//note that the return flows are subtracted because they reduce the amount of water that is to be abstracted from groundwater
					//the array groundwater_to_take is also altered for effects of GW abstraction inside AllocateWaterToUsers
					//thus a value of groundwater_to_take can be either positive (abstraction>returnflow) or negative (returnflow>abstraction)
			}

		}

		// make some entries in the output tables
		Append_To_Output_Tables(Timestep, dt, NumStreamNode, yyyymmdd, hhmmss, NumLink, NumNode, NSteps, NumUser);

		// Artificial drainage and write outputs each step
		scalefactor = 1.0/(double)dt;   // Scale factor to change units.  This was previously done in Append_to_output_tables subroutine
		Write_OutputLine_Eigen(oFile[8],  "results/Artificial_Drainage_cms.txt", Timestep, ArtDrainageTopnet, NumDrainage, scalefactor);
		scalefactor = 1.0;
		Write_OutputLine(oFile[9],  "results/Precipitation_mm.txt",              Timestep, precipitation,     NumDrainage, scalefactor);
		scalefactor = 1.0;
		Write_OutputLine(oFile[10], "results/Evaporation_mm.txt",                Timestep, evaporation,       NumDrainage, scalefactor);
		scalefactor = 1.0/(double)dt;
		// Overwrite BaseflowTopnet for writing.  It is not needed any more
		for (n = 0; n < NumDrainage; n++) {
			BaseflowTopnet(n) = Baseflow[Timestep-1].Rate[n];
		}
		Write_OutputLine_Eigen(oFile[11], "results/Baseflow_cms.txt",            Timestep, BaseflowTopnet,    NumDrainage, scalefactor);
		scalefactor = 1.0/(double)dt;
		for (n = 0; n < NumDrainage; n++) {
			RunoffTopnet(n) = Runoff[Timestep-1].Rate[n];
		}
		Write_OutputLine_Eigen(oFile[12], "results/TotalRunoff_noWithdrawal_cms.txt", Timestep, RunoffTopnet,      NumDrainage, scalefactor);
		Write_OutputLocalContributions(oFile[13], oFile[14], NumStreamNode, NumDrainage, BaseflowTopnet, RunoffTopnet, Timestep, scalefactor);

		//    if (mod(Timestep,10).eq.0) write(6,*)Timestep
		// end do !loop on Timestep

		// END OF do this Timestep !!!!!!

	} else if (Timestep < 0) {

		// FINAL OUTPUTS

		t2 = (double)clock()/(double)CLOCKS_PER_SEC;
		cerr << t2 - t1 << " seconds to calculate flows\n";
		cerr << "writing files\n";

		//now write the output tables
		Write_Static_Output_Tables("results", NumUserSourceReturn);
		//call CreateMoreOutputTables(NumUser,RunControl%NumTimesteps)
		Write_TimeVaryingOutput_Tables("results", NumUser, RunControl.NumTimesteps, NumNode, NumLink, NumReturnFlow, NumWWTP);

		//  Writing as we go close files
		scalefactor = 1.0/(double)dt;   // Scale factor to change units.
		Write_OutputLine_Eigen(oFile[8], "results/Artificial_Drainage_cms.txt", Timestep, ArtDrainageTopnet, NumDrainage, scalefactor);
		scalefactor = 1.0;
		Write_OutputLine(oFile[9], "results/Precipitation_mm.txt",              Timestep, precipitation,     NumDrainage, scalefactor);
		scalefactor = 1.0;
		Write_OutputLine(oFile[10], "results/Evaporation_mm.txt",               Timestep, evaporation,       NumDrainage, scalefactor);
		scalefactor = 1.0/(double)dt;
		Write_OutputLine_Eigen(oFile[11], "results/Baseflow_cms.txt",                 Timestep, BaseflowTopnet,    NumDrainage, scalefactor);
		Write_OutputLine_Eigen(oFile[12], "results/TotalRunoff_noWithdrawal_cms.txt", Timestep, RunoffTopnet,      NumDrainage, scalefactor);
		Write_OutputLocalContributions(oFile[13], oFile[14], NumStreamNode, NumDrainage, BaseflowTopnet, RunoffTopnet, Timestep, scalefactor);
		t3 = (double)clock()/(double)CLOCKS_PER_SEC;
		cerr << fixed << setw(12) << setprecision(6) << t3 - t2 << " seconds to write output files\n";

		// END OF FINAL OUTPUTS
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving watermgmt(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}


int timecalcs(const int Timestep, const int dt, const long int i8startsecs,
	int &doy, int &ThisMonth, int &ThisDay, int &yyyymmdd, int &hhmmss)
{
	int hh, mm, ss, yyyymmdd2, ihour2, ThisYear;
	long int i8nowsecs, i8StartofCalendarYearsecs;

	//inputs: Timestep, dt, i8startsecs
	//outputs: DOY, ThisMonth, ThisDay, yyyymmdd,hhmmss
	//figure out ThisMonth, ThisDate (and DOY in case we have DailyDemandFraction varying by DOY)
	i8nowsecs = i8startsecs + (long int)((Timestep - 1)*dt);
	td81micdh(yyyymmdd, hhmmss, i8nowsecs); //seconds to now since start of 1 Jan 1940
	datevec(yyyymmdd, hhmmss, ThisYear, ThisMonth, ThisDay, hh, mm, ss);
	undatevec(ThisYear, 1, 1, 0, 0, 0, yyyymmdd2, ihour2);
	td8micsec(yyyymmdd2, ihour2, i8StartofCalendarYearsecs);
	doy = (int)((i8nowsecs - i8StartofCalendarYearsecs)/86400);

	return 0;
}


