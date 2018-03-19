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
using namespace constant_definitions; //define names for numeric codes (types of nodes, links, users, etc)
using namespace data_array;
using namespace input_structures;
using namespace other_structures;

// EXAMPLE DATA
int read_inputs(const string dirname, const int dt, const int StartDateTopnet, int &NumDrainage,
	int &NumStreamNode, int &NumMeasuredFlowInfo, int &NumMeasuredFlowData,
	int &NumReservoir, int &NumUser, int &NumSource, int &NumRights, int &NumSourceMixing,
	int &NumSeasonsDefn, int &NumReturnFlow, int &NumMonthlyDemand, const int NumRunoff,
	int &NumBaseflow, int &NumWWTP)
{
	string fileName;
	int nrows, expected_numcols;
	int *ifound_1, *ifound_2, nfound, *ifound;
	double sumarea, area;
	int ncommentlines;
	int i, i_outlet, j, jj, nnode, n, nfound2, k, kk, nSites, ioffset;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> read_inputs(" << ncalls << ")" << std::endl;
    }
	caller = "read_inputs";
#endif
	ncommentlines = 2;
	fileName = dirname + "/" + "WatermgmtControl.txt";
	expected_numcols = 5;
	read_struct_from_text(fileName, expected_numcols, ncommentlines, nrows); //NumTimesteps	AllocationMode	StartDate	StartOfWaterYear

	RunControl.NumTimesteps       = real_array(0,0); //ignore
	RunControl.AllocationMode     = real_array(0,1); //use
	RunControl.StartDate          = real_array(0,2); //ignore
	RunControl.StartOfWaterYearmm = real_array(0,3); //use
	RunControl.StartOfWaterYeardd = real_array(0,4); //use

	//call datevec(RunControl%StartDate,90000,	&
	//			RunControl%StartDateyyyy,RunControl%StartDatemm,RunControl%StartDatedd,hh,mi,ss)

	//call datevec(RunControl%StartOfWaterYear,90000,	&
	//  RunControl%StartOfWaterYearyyyy,RunControl%StartOfWaterYearmm,RunControl%StartOfWaterYeardd,	&
	//  hh,mi,ss)

	fileName = dirname + "/" + "basinpars.txt";
	expected_numcols = 45;
	ncommentlines = 1;
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumDrainage); //CatchID,DrainageID,NodeID,CatchArea
	Drainage.resize(NumDrainage);

	//CatchID, DownCatchID, DrainID, NodeId, Reach_number, Outlet_X, Outlet_Y, direct_area, f, k, dth1, dth2,
	// soilc, c, psif, chv, can_capacity, cr, Albedo, Lapse_rate, average_elevation, lambda, std_dev_of_lambda
	for (i = 0; i < NumDrainage; i++) {
		Drainage(i).DrainageID     = real_array(i,0);
		Drainage(i).DSDrainage     = real_array(i,1);
		Drainage(i).RealDrainageID = real_array(i,2);
		Drainage(i).NodeID         = real_array(i,3);
		Drainage(i).CatchArea      = real_array(i,7);
		Drainage(i).ifound_in.resize(1000);
		Drainage(i).ifound_taken.resize(1000);
	}

	fileName = dirname + "/" + "nodelinks.txt";
	expected_numcols = 10;
	ncommentlines = 1;
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumStreamNode); //NodeID,DownNodeID,DrainageID,NodeCatchID,DoutFlag,Area
	StreamNode = new StreamNodeType[NumStreamNode];
	for (i = 0; i < NumStreamNode; i++) {
		StreamNode[i].NodeID         = real_array(i,0);
		StreamNode[i].DownNodeID     = real_array(i,1);
		StreamNode[i].RealDrainageID = real_array(i,2);
		StreamNode[i].ProjNodeId     = real_array(i,3);
		StreamNode[i].DOutFlag       = real_array(i,4);
		StreamNode[i].LocalArea      = real_array(i,6);
		StreamNode[i].TotalArea      = real_array(i,7);
		StreamNode[i].X              = real_array(i,8);
		StreamNode[i].Y              = real_array(i,9);
		StreamNode[i].AreaInDrainage = StreamNode[i].TotalArea;
	}

	//Drainage(i)%DSDrainage= find the StreamNode which is outlet of this Drainage, go to Downnode, find CatchID
	for (i = 0; i < NumDrainage; i++) {
		// first find the streamnode at the outlet
		// DGT interprets this loop as identifying the Topnet identifier DrainageID for each node
		// from the RealDrainageID that was read in from Nodelinks
		// find1()
		nfound = 0; //none found
		ifound = new int[NumStreamNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumStreamNode; n++) {
			if (StreamNode[n].RealDrainageID == Drainage(i).RealDrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		i_outlet = ifound[0]; //  This line appears redundant
		//i_outlet is the row in the StreamNode table for the outlet of this drainage
		//StreamNode(i_outlet)%DownNodeID is the node which that StreamNode flows to
		for (j = 0; j < nfound; j++) {
			StreamNode[ifound[j]-1].DrainageID = Drainage(i).DrainageID;
		}
		delete [] ifound;
	}

	for (i = 1; i <= NumDrainage; i++) {
		//first find the streamnode at the outlet of this drainage
		// find 2 returns in ifound the list of indices for which firstarg(i)=2ndarg and
		// 3rdarg(i)=4tharg.  The search is from 1 to 5tharg

		// find2()
		nfound = 0; //none found
		ifound = new int[NumStreamNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumStreamNode; n++) {
			if (StreamNode[n].RealDrainageID == Drainage(i-1).RealDrainageID && StreamNode[n].DOutFlag == IsDrainageOutlet) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		i_outlet = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		//now find all streamnodes in this drainage
		// find1()
		nfound = 0; //none found
		ifound = new int[NumStreamNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumStreamNode; n++) {
			if (StreamNode[n].RealDrainageID == Drainage(i-1).RealDrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		nnode = nfound;
		if (nnode == 1) {
			StreamNode[i_outlet-1].FracRunoff = 1; //there's only one node in this drainage, so it was easy
			StreamNode[i_outlet-1].AreaInDrainage = StreamNode[i_outlet-1].LocalArea;
		} else if (nnode > 1) { //if there are other nodes than the one at the outlet, then we have work to do
			//find the fraction of this drainage area which drains through each stream node in this drainage
			ifound_1 = new int[nfound];
			for (n = 0; n < nfound; n++) {
				ifound_1[n] = ifound[n];
			}
			sumarea = 0.0;
			for (j = 1; j <= nnode; j++) {
				sumarea += StreamNode[ifound_1[j-1]-1].LocalArea;
			}
			delete [] ifound;
			for (jj = 1; jj <= nnode; jj++) {
				j = ifound_1[jj-1];
				//is this node getting some area directly from another drainage?
				// find2a()
				nfound = 0; //none found
				ifound = new int[NumStreamNode];
				ifound[nfound] = 0;
				for (n = 0; n < NumStreamNode; n++) {
					if (StreamNode[n].DownNodeID == StreamNode[j-1].NodeID && StreamNode[n].DrainageID != Drainage(i-1).DrainageID) { 	//A&(~B)
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}

				ifound_2 = new int[nfound];
				for (n = 0; n < nfound; n++) {
					ifound_2[n] = ifound[n];
				}
				nfound2 = nfound;
				delete [] ifound;
				if (nfound2 > 0) { //for each of the external areas,
					for (kk = 1; kk <= nfound2; kk++) {
						area = StreamNode[ifound_2[kk-1]-1].TotalArea;
						//now work downstream removing this area from all nodes within this drainage
						k = j;
						while (StreamNode[k-1].DrainageID == Drainage(i-1).DrainageID && nfound > 0) {
							StreamNode[k-1].AreaInDrainage = StreamNode[k-1].AreaInDrainage - area;
							// find1()
							nfound = 0; //none found
							ifound = new int[NumStreamNode];
							ifound[nfound] = 0;
							for (n = 0; n < NumStreamNode; n++) {
								if (StreamNode[n].NodeID == StreamNode[k-1].DownNodeID) {
									nfound++;
									ifound[nfound-1] = n+1;
								}
							}
							if (nfound == 0)
								break;  // DGT 5/27/12 because with nfound = 0 k was undefined and caused subscript out of range in do while above
							k = ifound[0];
							delete [] ifound;
						} //while
					} //kk
				}
				delete [] ifound_2;
			} //jj
			for (jj = 1; jj <= nnode; jj++) {	 // have to delay this operation until all u/s areas have been subtracted
				j = ifound_1[jj-1];
				StreamNode[j-1].FracRunoff = StreamNode[j-1].AreaInDrainage/sumarea;
			} //jj
			delete [] ifound_1;
		}
	} //i

	//	fileName=dirname(1:len_trim(dirname)) // '\' // 'Drainage.txt'
	//	expected_numcols=3
	//	call read_struct_from_text(fileName,expected_numcols,ncommentlines,NumDrainage) !DrainageID	DSDrainage	ColumnInRunoff
	//	allocate (Drainage(NumDrainage))
	//	Drainage(:)%DrainageID=real_array(:,1)
	//	Drainage(:)%DSDrainage=real_array(:,2)
	//	Drainage(:)%ColumnInRunoff=real_array(:,3)
	//	Drainage(:)%FractionOfRunoff=real_array(:,4)
#ifdef LNWB
    fileName = dirname + "/" + "MeasuredFlowInfoLNWB.txt";
#else
	fileName = dirname + "/" + "MeasuredFlowInfo.txt";
#endif
	expected_numcols = 4;
	ncommentlines = 2;

	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumMeasuredFlowInfo); //MeasuredFlowID	DrainageID	ColInMeasFlow
	MeasuredFlowInfo = new MeasuredFlowInfoType[NumMeasuredFlowInfo];
	for (i = 0; i < NumMeasuredFlowInfo; i++) {
		MeasuredFlowInfo[i].MeasuredFlowID	= real_array(i,0);
		MeasuredFlowInfo[i].DrainageID		= real_array(i,1);
		MeasuredFlowInfo[i].ColInMeasFlow	= real_array(i,2);
		MeasuredFlowInfo[i].ScalingFactor	= real_array(i,3);
	}

	//convert any reference to a DrainageID into a CatchID
	for (i = 0; i < NumMeasuredFlowInfo; i++) {
		MeasuredFlowInfo[i].RealDrainageID = MeasuredFlowInfo[i].DrainageID;
		// find1()
		nfound = 0; //none found
		ifound = new int[NumDrainage];
		ifound[nfound] = 0;
		for (n = 0; n < NumDrainage; n++) {
			if (Drainage(n).RealDrainageID == MeasuredFlowInfo[i].RealDrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		MeasuredFlowInfo[i].DrainageID = Drainage(ifound[0]-1).DrainageID;
		delete [] ifound;
	}
#ifdef LNWB
    fileName = dirname + "/" + "bndryflowLNWB.dat";
#else
	fileName = dirname + "/" + "bndryflow.dat";
#endif
	expected_numcols = -1;
	ncommentlines = 4;
	read_bdryflow(fileName, expected_numcols, ncommentlines, NumMeasuredFlowData); //Flow1	Flow2 Date	Time
	nSites = expected_numcols - 2; // nSites = size(dble_array, DIM=2)-2
	// find1()
	nfound = 0; //none found
	ifound = new int[NumMeasuredFlowData];
	ifound[nfound] = 0;
	for (n = 0; n < NumMeasuredFlowData; n++) {
		if ((int)dble_array[n][nSites] == StartDateTopnet) {
			nfound++;
			ifound[nfound-1] = n+1;
		}
	}
	ioffset = ifound[0] - 1;
	delete [] ifound;
	NumMeasuredFlowData = NumMeasuredFlowData - ioffset;
	MeasuredFlowData = new MeasuredFlowDataType[NumMeasuredFlowData];

	for (i = 1; i <= NumMeasuredFlowData; i++) {
		for (j = 1; j <= nSites; j++) {
			MeasuredFlowData[i-1].Flow[j-1]	= dble_array[i+ioffset-1][j-1]*dt; // cumecs to m^3/d
		}
		MeasuredFlowData[i-1].Timestep	= i;
		MeasuredFlowData[i-1].Date		= dble_array[i+ioffset-1][nSites];
		MeasuredFlowData[i-1].Hour		= dble_array[i+ioffset-1][nSites+1];
	}
	//use the elements of integer_array to change the value of MeasuredFlowInfo(:)%ColInMeasFlow from a reference
	//into a column number referring to MeasuredFlowData
	for (i = 1; i <= NumMeasuredFlowInfo; i++) {
		// find1()
		nfound = 0; //none found
		ifound = new int[nSites];
		ifound[nfound] = 0;
		for (n = 0; n < nSites; n++) {
			if (integer_array(n,0) == (int)(MeasuredFlowInfo[i-1].ColInMeasFlow)) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		j = ifound[0];
		delete [] ifound;
		MeasuredFlowInfo[i-1].ColInMeasFlow = j;
	}

	//	fileName=dirname + "/" + "Runoff.txt";
	//	expected_numcols=-1
	//	call read_struct_from_text(fileName,expected_numcols,ncommentlines,NumRunoff) !Flow1	Flow2 Date	Time
	//	allocate (Runoff(NumRunoff))
	//	nSites = size(real_array,DIM=2)-2
	//	for (i = 1; i <= NumRunoff
	//		allocate (Runoff[i-1].Runoff(nSites))
	//		do j = 1,nSites
	//			Runoff[i-1].Rate(j)	=real_array(i,j)
	//		end do
	//		Runoff[i-1].Timestep =i
	//	end do
	//	Runoff(:)%Date		=real_array(:,nSites+1)
	//	Runoff(:)%Hour		=real_array(:,nSites+2)

	fileName = dirname + "/" + "Reservoir.txt";
	expected_numcols = 10;
	ncommentlines = 2;
	//ReservoirID	DrainageID	InOffStream	RightID	MaxStore	InitialStore	MinStore	MaxInflow	MaxWithdrawal	MinEnvRelease	LossRate
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumReservoir);
	Reservoir = new ReservoirType[NumReservoir];
	for (i = 0; i < NumReservoir; i++) {
		Reservoir[i].ReservoirID   = real_array(i,0);
		Reservoir[i].DrainageID    = real_array(i,1);
		Reservoir[i].InOffStream   = real_array(i,2);
		Reservoir[i].RightID       = real_array(i,3);
		Reservoir[i].StoreMax      = real_array(i,4);
		Reservoir[i].StoreInitial  = real_array(i,5);
		Reservoir[i].StoreMin      = real_array(i,6);
		Reservoir[i].MaxInflow     = real_array(i,7);
		Reservoir[i].MaxWithdrawal = real_array(i,8);
		Reservoir[i].MinEnvRelease = real_array(i,9);
		Reservoir[i].ifound_in.resize(1000);
		Reservoir[i].ifound_taken.resize(1000);
		Reservoir[i].ifound_in    = 0;
		Reservoir[i].ifound_taken = 0;
		//	Reservoir(:)%LossRate=real_array(:,11)
	}

	//convert references to a DrainageID into the CatchID
	for (i = 1; i <= NumReservoir; i++) {
		Reservoir[i-1].RealDrainageID = Reservoir[i-1].DrainageID;
		// find1()
		nfound = 0; //none found
		ifound = new int[NumDrainage];
		ifound[nfound] = 0;
		for (n = 0; n < NumDrainage; n++) {
			if (Drainage(n).RealDrainageID == Reservoir[i-1].RealDrainageID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		Reservoir[i-1].DrainageID = Drainage(ifound[0]-1).DrainageID;
		delete [] ifound;
	}

	fileName = dirname + "/" + "User.txt";
	expected_numcols = 15;
	ncommentlines = 2;
      //UserID UserType POU_ID DemandType DemandVble DemandRate ReturnFlowID SourceMixingID NumSources SourceID1 RightID1
      // SourceID2 RightID2 SourceID3 RightID3
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumUser);
	User = new UserType[NumUser+2*NumReservoir];
	//allocate (User(NumUser+2*NumReservoir))
	for (i = 0; i < NumUser; i++) {
		User[i].UserID           = real_array(i,0);
		User[i].UsersType        = real_array(i,1);
		User[i].RealPOU_ID       = real_array(i,2);
		User[i].DemandVble       = real_array(i,3);
		User[i].DemandRate       = real_array(i,4);
		User[i].InYearDemandType = real_array(i,5);
		User[i].ReturnFlowID     = real_array(i,6);
		User[i].SourceMixingID   = real_array(i,7);
		User[i].NumSources       = real_array(i,8);
	}

	for (i = 1; i <= NumUser; i++) {
		for (j = 1; j <= User[i-1].NumSources; j++) {
			User[i-1].SourceID[j-1] = real_array(i-1,9+2*j-2);
			User[i-1].RightID[j-1]  = real_array(i-1,9+2*j-1);
		}
	}

	//convert references to a DrainageID into the CatchID
	for (i = 1; i <= NumUser; i++) {
		// find1()
		nfound = 0; //none found
		ifound = new int[NumDrainage];
		ifound[nfound] = 0;
		for (n = 0; n < NumDrainage; n++) {
			if (Drainage(n).RealDrainageID == User[i-1].RealPOU_ID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		User[i-1].POU_ID = Drainage(ifound[0]-1).DrainageID;
		delete [] ifound;
	}

	fileName = dirname + "/" + "Source.txt";
	expected_numcols = 5;
	ncommentlines = 2;
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumSource); //SourceID	Type	SourceLocationID	PhysicalDailyMax	PhysicalAnnMax
	Source = new SourceType[NumSource+2*NumReservoir];
	for (i = 0; i < NumSource; i++) {
		Source[i].SourceID         = real_array(i,0);
		Source[i].Type             = real_array(i,1);
		Source[i].SourceLocationID = real_array(i,2);
		Source[i].PhysicalDailyMax = real_array(i,3);
		Source[i].PhysicalAnnMax   = real_array(i,4);
	}
	//convert references to a DrainageID into the CatchID
	for (i = 1; i <= NumSource; i++) {
		Source[i-1].RealSourceLocationID = Source[i-1].SourceLocationID;
		// find1()
		nfound = 0; //none found
		ifound = new int[NumDrainage];
		ifound[nfound] = 0;
		for (n = 0; n < NumDrainage; n++) {
			if (Drainage(n).RealDrainageID == Source[i-1].RealSourceLocationID) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		if (nfound == 1) {
			Source[i-1].SourceLocationID = Drainage(ifound[0]-1).DrainageID;
		}
		delete [] ifound;
	}

	fileName = dirname + "/" + "Rights.txt";
	expected_numcols = 4;
	ncommentlines = 2;
	 //RightID	PurposeCode	PriorityDate	LegalDailyMax	LegalAnnMax
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumRights);
	Rights = new RightsType[NumRights];
	for (i = 0; i < NumRights; i++) {
		Rights[i].RightID       = real_array(i,0);
		Rights[i].PriorityDate  = real_array(i,1);
		Rights[i].LegalDailyMax = real_array(i,2);
		Rights[i].LegalAnnMax   = real_array(i,3);
	}

	fileName = dirname + "/" + "SourceMixing.txt";
	expected_numcols = 6;
	ncommentlines = 2;
      //SourceMixingID	UsersSourceNum	Units	Amount	SeasonNumber	SeasonsDefnID
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumSourceMixing);
	SourceMixing = new SourceMixingType[NumSourceMixing];
	for (i = 0; i < NumSourceMixing; i++) {
		SourceMixing[i].SourceMixingID = real_array(i,0);
		SourceMixing[i].UsersSourceNum = real_array(i,1);
		SourceMixing[i].Units          = real_array(i,2);
		SourceMixing[i].Amount         = real_array(i,3);
		SourceMixing[i].SeasonNumber   = real_array(i,4);
		SourceMixing[i].SeasonsDefnID  = real_array(i,5);
	}

	fileName = dirname + "/" + "SeasonsDefn.txt";
	expected_numcols = 5;
	ncommentlines = 2;
	//SourceMixingID	UsersSourceNum	Units	Amount	SeasonNumber	SeasonsDefnID
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumSeasonsDefn);
	SeasonsDefn = new SeasonsDefnType[NumSeasonsDefn];
	for (i = 0; i < NumSeasonsDefn; i++) {
		SeasonsDefn[i].SeasonsDefnID = real_array(i,0);
	}
	for (i = 1; i <= NumSeasonsDefn; i++) {
		for (j = 1; j <= 4; j++) {
			SeasonsDefn[i-1].StartDaySeason[j-1] = real_array(i-1,j);
		}
	}

	fileName = dirname + "/" + "ReturnFlow.txt";
	expected_numcols = 11;
	ncommentlines = 2;
	//ReturnFlowID	NumReturnFlows	ReturnFlowUnits	ReturnFlowAmt1 ReturnFlowType1	ReturnFlowLocn1	ReturnFlowAmt2	ReturnFlowType2	ReturnFlowLocn2
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumReturnFlow);
	ReturnFlow = new ReturnFlowType[NumReturnFlow + 2*NumReservoir];
	WWTP_list = new int[2*(NumReturnFlow + 2*NumReservoir)];
	for (i = 0; i < NumReturnFlow; i++) {
		ReturnFlow[i].ReturnFlowID     = real_array(i,0);
		ReturnFlow[i].NumReturnFlows   = real_array(i,1);
		ReturnFlow[i].ReturnFlowsUnits = real_array(i,2);
	}
	for (i = 1; i <= NumReturnFlow; i++) {
		j = ReturnFlow[i-1].NumReturnFlows;
		for (j = 1; j <= User[i-1].NumSources; j++) {
			ReturnFlow[i-1].ReturnFlowsAmt[j-1]  = real_array(i-1,3+4*j-4);
			ReturnFlow[i-1].ReturnFlowsType[j-1] = real_array(i-1,3+4*j-3);
			ReturnFlow[i-1].ReturnFlowsLocn[j-1] = real_array(i-1,3+4*j-2);
			ReturnFlow[i-1].WWTP_ID[j-1]         = real_array(i-1,3+4*j-1);
			// find1()
			nfound = 0; //none found
			ifound = new int[NumWWTP];
			ifound[nfound] = 0;
			for (n = 0; n < NumWWTP; n++) {
				if (WWTP_list[n] == ReturnFlow[i-1].WWTP_ID[j-1]) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			delete [] ifound;
			if (ReturnFlow[i-1].WWTP_ID[j-1] > 0 && nfound == 0) {
				NumWWTP++;
				WWTP_list[NumWWTP-1] = ReturnFlow[i-1].WWTP_ID[j-1];
			}
			ReturnFlow[i-1].RealReturnFlowsLocn[j-1] = ReturnFlow[i-1].ReturnFlowsLocn[j-1];
			if (ReturnFlow[i-1].ReturnFlowsLocn[j-1] > 0) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumDrainage];
				ifound[nfound] = 0;
				for (n = 0; n < NumDrainage; n++) {
					if (Drainage(n).RealDrainageID == ReturnFlow[i-1].RealReturnFlowsLocn[j-1]) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				for (n = 0; n < MaxNumReturnFlows; n++) {
					ReturnFlow[i-1].ReturnFlowsLocn[n] = Drainage(ifound[0]-1).DrainageID;
				}
				delete [] ifound;
			}
		}
	}
	for (i = (NumReturnFlow+1); i <= (NumReturnFlow+2*NumReservoir); i++) {
		ReturnFlow[i-1].ReturnFlowsAmt[0]      = 0.0;
		ReturnFlow[i-1].ReturnFlowsType[0]     = 0;
		ReturnFlow[i-1].ReturnFlowsLocn[0]     = 0;
		ReturnFlow[i-1].RealReturnFlowsLocn[0] = 0;
		ReturnFlow[i-1].WWTP_ID[0]             = 0;
	}

	//	fileName=dirname + "/" + "DemandCoefficients.txt'
	//	expected_numcols=2
	//	call read_struct_from_text(fileName,expected_numcols,ncommentlines,NumDemandCoefficients) !DemandType	DemandPerUnitPerDay
	//	allocate (DemandCoefficients(NumDemandCoefficients))
	//	DemandCoefficients(:)%DemandType=real_array(:,1)
	//	DemandCoefficients(:)%DemandPerUnitPerDay=real_array(:,2)

	fileName = dirname + "/" + "MonthlyDemandFraction.txt";
	expected_numcols = 13;
	ncommentlines = 2;
	//InYearDemandType	Month1	Month2	Month3	Month4	Month5	Month6	Month7	Month8	Month9	Month10	Month11	Month12
	read_struct_from_text(fileName, expected_numcols, ncommentlines, NumMonthlyDemand);
	MonthlyDemand = new MonthlyDemandType[NumMonthlyDemand];
	for (i = 0; i < NumMonthlyDemand; i++) {
		MonthlyDemand[i].InYearDemandType = real_array(i,0);
	}
	for (i = 1; i <= NumMonthlyDemand; i++) {
		for (j = 1; j <= 12; j++) {
			MonthlyDemand[i-1].Month[j-1] = real_array(i-1,j);
		}
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving read_inputs(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}

int datevec(const int yyyymmdd, const int hhmmss, int &yyyy, int &mm, int &dd, int &hh, int &mi, int &ss)
{
	//split yyyymmdd,hhmmss into components yyyy,mm,dd,hh,mi,ss

	dd = yyyymmdd % 100;
	mm = (yyyymmdd - dd)/100 % 100;
	yyyy = (int)(yyyymmdd/10000);
	ss = hhmmss % 100;
	mi = (hhmmss-ss)/100 % 100;
	hh = int(hhmmss/10000);

	return 0;
}


int undatevec(const int yyyy, const int mm, const int dd, const int hh, const int mi, const int ss, int &yyyymmdd, int &hhmmss)
{
	// combine components yyyy,mm,dd,hh,mi,ss into yyyymmdd,hhmmss

	yyyymmdd = 10000*yyyy + 100*mm + dd;
	hhmmss = 10000*hh + 100*mi + ss;

	return 0;
}

