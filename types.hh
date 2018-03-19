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

#ifndef TYPES_HH
#define TYPES_HH

#include <Eigen/Dense>

const int MaxMeasuredFlowSites = 100;
const int MaxNumTimesteps      = 366*60;
const int MaxNumSources        = 10;
const int MaxNumReturnFlows    = 10;
const int MaxNumDrainage       = 500;

const double BigReal = 1.0e20;
const int MaxNtimeSteps = MaxNumTimesteps;
// MaxNumSources and MaxNumReturnFlows are per each user.  There can be more than these numbers over all users

struct RunControlType {
	int NumTimesteps;
	int AllocationMode;
	int StartDate;
	int StartDateyyyy;
	int StartDatemm;
	int StartDatedd;
	int StartDatehh;
	int StartOfWaterYear;
	int StartOfWaterYearyyyy;
	int StartOfWaterYearmm;
	int StartOfWaterYeardd;
	int StartOfWaterYearhh;
}; // RunControlType

struct DrainageType {
       int DrainageID;			//sequential numbers 1 ... NumDrainage, not necessarily ordered
       int RealDrainageID;		//numbering used by client
       int DSDrainage;			//the CatchID of the drainage d/s of this drainage
       int NodeID;				//the NodeID of the StreamNode at the drainage outlet
       double CatchArea;		//The area of the Drainage (needs to be in same units as StreamNodeType%Area)
       int isn;
       int n_in;
       int n_taken;
       int j_unallocated;
       Eigen::Array<int,Eigen::Dynamic,1> ifound_in;
       Eigen::Array<int,Eigen::Dynamic,1> ifound_taken;
};

struct StreamNodeType {
       int NodeID;				//arbitrary numbering
       int DownNodeID;			//the NodeID of the stream node just downstream of this one
       int RealDrainageID;		//the numbering used by client for the drainage containing this node
       int DrainageID;			//the numbering used by Topnet for the drainage containing this node
       int ProjNodeId;			//the project's externally-set ID for this node
       int DOutFlag;			//1 if the node is drainage outlet, 0 otherwise
       double LocalArea;		//The area draining only to this node (needs to be in same units as DrainageType%CatchArea)
       double TotalArea;		//The area draining to this node (needs to be in same units as DrainageType%CatchArea)
       double FracRunoff;		//The fraction of the runoff in this node's drainage which flows through this node
       double X;				//The X coord of this node
       double Y;				//The Y coord of this node
       double AreaInDrainage; 	//The area drainaing to this node which is in this drainage
};

struct MeasuredFlowInfoType {
       int MeasuredFlowID;
       int DrainageID;
       int RealDrainageID;
       int ColInMeasFlow;
       double	ScalingFactor;
};

struct MeasuredFlowDataType {
       int Timestep;
       int Hour;
       int Date;
       double Flow[MaxMeasuredFlowSites];
};

struct ReservoirType {
       int ReservoirID;
       int DrainageID;
       int RealDrainageID;
       int InOffStream;
       int RightID;
       double StoreMax;
       double StoreInitial;
       double StoreMin;
       double MaxInflow;
       double MaxWithdrawal;
       double MinEnvRelease;
       double LossRate;
       int isn;
       int n_in;
       int n_taken;
       Eigen::Array<int,Eigen::Dynamic,1> ifound_in;
       Eigen::Array<int,Eigen::Dynamic,1> ifound_taken;
};

struct UserType {
       int UserID;
       int UsersType;
       int POU_ID;
       int RealPOU_ID;
       double	DemandVble;
       double	DemandRate;
       int InYearDemandType;
       double	Demand[MaxNumTimesteps];
       double	Deficit[MaxNumTimesteps];
       double	Withdrawal[MaxNumTimesteps];
       double	VolumeToDateSource[MaxNumSources];
       int ReturnFlowID;
       int SourceMixingID;
       int NumSources;
       int SourceID[MaxNumSources];
       int RightID[MaxNumSources];
       double	DemandToday;
       int NodeNumber;
       int LinkSourceToUser[MaxNumSources];
       int LinkSourceToSink;
       int LinkUserToReturnflow[MaxNumReturnFlows];
       int j_srcmx[MaxNumSources];
};

struct SourceType {
		int SourceID;
		int Type;
		int SourceLocationID;
		int RealSourceLocationID;
		double PhysicalDailyMax;
		double PhysicalAnnMax;
};

struct RightsType {
		int RightID;
		int PriorityDate;
		double LegalDailyMax;
		double LegalAnnMax;
};

struct SourceMixingType {
		int SourceMixingID;
		int UsersSourceNum;
		int Units;
		double Amount;
		int SeasonNumber;
		int SeasonsDefnID;
};

struct SeasonsDefnType {
		int SeasonsDefnID;
		int StartDaySeason[4];
};

struct ReturnFlowType {
		int ReturnFlowID;
		int NumReturnFlows;
		int ReturnFlowsUnits;
		double	ReturnFlowsAmt[MaxNumReturnFlows];
		int ReturnFlowsType[MaxNumReturnFlows];
		int ReturnFlowsLocn[MaxNumReturnFlows];
		int RealReturnFlowsLocn[MaxNumReturnFlows];
		int WWTP_ID[MaxNumReturnFlows];
};

struct MonthlyDemandType {
		int InYearDemandType;
		double Month[12];
};

struct RunoffType {
		int Timestep;
		int Date;
		int Hour;
		double Rate[MaxNumDrainage];
};

struct NodeType {
		std::string Title;
		int Type;
		int IntExt;
		double StoreMin;
		double StoreMax;
		double Store;
		double StoreOld;
		int DrainageID;
		int SelfID;
};

struct LinkType {
		std::string Title;
		int LinkCode;
		int IntExtCode;
		int USNode;
		int DSNode;
		double Flow;
		int ReturnFlowID;
};

struct UserSourceTableType {
		int UserID;
		int SourceCounter;
		int PriorityDate;
		int DrainageID;
};

struct StreamFlowLinksType {
		int DrainageID;
		int LinkID;
};

struct DrainageIDType {
		int TopnetID;
		int DrainageID;
};

struct UserFlowLinksType {
		int UserID;
		int SourceDrainage;
		int DestinationDrainage;
		int TakeLinkID;
		int ReturnLinkID;
		int UsersType;
};

struct DrainageInfoType {
		int UserSrceReturn_ID;
		int User_ID;
		int UsersType;
		int SourceDrainageID;
		int ReturnLinkID;
		int DestinationDrainageID;
		int WWTP_ID;
		int RF_counter;
		int SourceID;
		int ReturnFlowID;
};

struct StaticOutputTableType {
	StreamFlowLinksType *StreamFlowLinks;
	int StreamFlowLinksSize;
	UserFlowLinksType   *UserFlowLinks;
	int UserFlowLinksSize;
	DrainageIDType      *DrainageID;
	int DrainageIDSize;
	DrainageInfoType    *DrainageInfo;
	int DrainageInfoSize;
};

#if false

find() prototypes:

		// find1()
		nfound = 0; //none found
		ifound = new int[xxx];
		ifound[nfound] = 0;
		for (n = 0; n < xxx; n++) {
			if (xxx[n].xxx == xxx) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		xxx = ifound[0];	//problem if nfound<>1
		delete [] ifound;


		// find2()
		nfound = 0; //none found
		ifound = new int[xxx];
		ifound[nfound] = 0;
		for (n = 0; n < xxx; n++) {
			if (xxx[n].xxx == xxx && xxx[n].xxx == xxx) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		xxx = ifound[0]; //problem if nfound<>1
		delete [] ifound;

		// find3()
		nfound = 0; //none found
		ifound = new int[xxx];
		ifound[nfound] = 0;
		for (n = 0; n < xxx; n++) {
			if (xxx[n].xxx == xxx && xxx[n].xxx == xxx && xxx[n].xxx == xxx) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		xxx = ifound[0]; //problem if nfound<>1
		delete [] ifound;

//=================================
				found = false;
				for (n = 0; n < NNN; n++) {
					if (XXX == XXX && XXX == XXX  && XXX == XXX) {
						User[i-1].LinkUserToReturnflow[m-1] = n+1;
						if (!found) {
							found = true;
						} else {
							cerr << "___(): Duplicate link found" << endl;
							exit(EXIT_FAILURE);
						}
					}
				}
				if (!found) {
					cerr << "___() find3() failure\n";
					exit(EXIT_FAILURE);
				} else {
					// reset
					found = false;
				}
//=================================
#endif

#endif
