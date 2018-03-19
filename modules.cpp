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

using namespace Eigen;

int NumDrainage;
int NumStreamNode;
int NumMeasuredFlowInfo;
int NumMeasuredFlowData;
int NumReservoir;
int NumUser;
int NumSource;
int NumRights;
int NumSourceMixing;
int NumSeasonsDefn;
int NumReturnFlow;
int NumMonthlyDemand;
int NumRunoff;
int NumBaseflow;
int NumWWTP;

namespace data_array {
	ArrayXXd real_array;
	ArrayXXi integer_array;
	double **dble_array;
	bool real_array_allocated    = false;
	bool integer_array_allocated = false;
	bool dble_array_allocated    = false;
}

namespace TimeVaryingOutput {
	//	real, allocatable :: Precipitation_mm(:,:) //Timestep,Location
	//	real, allocatable :: Evaporation_mm(:,:) //Timestep,Location
	//	real, allocatable :: TotalRunoff_cms(:,:) //Timestep,Location
	//	real, allocatable :: Baseflow_cms(:,:) //Timestep,Location
	//	real, allocatable :: Artificial_Drainage(:,:) //Timestep,Location
	ArrayXXd FlowInLinks_cms;			//Timestep,Location
	double *FlowAtStreamNodes_cms;		//Location
	ArrayXXd ReservoirStorage_m3;		//Timestep,Location
	ArrayXXi DateTime_yyyymmdd_hhmmss;		//Timestep, yyyymmdd/hhmmss
}

namespace input_structures {
	RunControlType RunControl;
	Array<DrainageType,Dynamic,1> Drainage;
	StreamNodeType *StreamNode;
	MeasuredFlowInfoType *MeasuredFlowInfo;
	MeasuredFlowDataType *MeasuredFlowData;
	RunoffType *Runoff;
	RunoffType *Baseflow;
	RunoffType *ArtDrainage;
	ReservoirType *Reservoir;
	UserType *User;
	SourceType *Source;
	RightsType *Rights;
	SourceMixingType *SourceMixing;
	SeasonsDefnType *SeasonsDefn;
	ReturnFlowType *ReturnFlow;
	//	DemandCoefficientsType *DemandCoefficients;
	MonthlyDemandType *MonthlyDemand;
}

namespace other_structures {
	NodeType              *Node;
	Array<LinkType,Dynamic,1> Link;
	UserSourceTableType   *UserSourceTable;
	StaticOutputTableType StaticOutput;
	NodeType              *NodeSave;
	Array<LinkType,Dynamic,1> LinkSave;
	int                   *UserSourceOrder;
	int                   *WWTP_list;
}
