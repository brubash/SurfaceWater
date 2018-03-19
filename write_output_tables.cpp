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
#include <sstream>

using namespace data_array;
using namespace constant_definitions; //define names for numeric codes (types of nodes, links, users, etc)
using namespace other_structures;
using namespace Eigen;
using namespace std;

int write_integers_to_text(const string filenm, const Array<string,Dynamic,1> &headers,
	const int nrows, const int ncols);
int write_struct_to_text(const string filenm, const Array<string,Dynamic,1> &headers, const int nrows, const int ncols,
	const int icall);

int Write_Static_Output_Tables(const string dirname, const int NumUserSourceReturn)
{
	// tables to write: StreamFlowLinks, DrainageID, DrainageInfo and UserFlowLinks
	Array<string,Dynamic,1> headers(10);
	string filenm;
	int ncols, nrows, n;

	//this StreamFlowLinks table has a row for every drainage, telling us which column of the LinkFlow
	//table will tell us the outflow for this drainage
	headers(0) = "DrainageID";
	headers(1)= "LinkID";
	ncols = 2;
    filenm = dirname + "/" + "StreamFlowLinks.txt";
	nrows = StaticOutput.StreamFlowLinksSize; // size(StaticOutput%StreamFlowLinks)

	integer_array.resize(nrows, ncols);
	for (n = 0; n < nrows; n++) {
		integer_array(n,0) = StaticOutput.StreamFlowLinks[n].DrainageID;
		integer_array(n,1) = StaticOutput.StreamFlowLinks[n].LinkID;
	}
	write_integers_to_text(filenm, headers, nrows, ncols);

	//this DrainageID table has a row for every drainage, telling us the connection between the TopnetID of a
	//drainage, and the externally specified Drainage Number
	headers(0) = "TopnetID";
	headers(1) = "DrainageID";
	ncols = 2;
    filenm = dirname + "/" + "DrainageID.txt";
	nrows = StaticOutput.DrainageIDSize; // size(StaticOutput%DrainageID)
	integer_array.resize(nrows, ncols);
	for (n = 0; n < nrows; n++) {
		integer_array(n,0) = StaticOutput.DrainageID[n].TopnetID;
		integer_array(n,1) = StaticOutput.DrainageID[n].DrainageID;
	}
    write_integers_to_text(filenm, headers, nrows, ncols);

	//this DrainageInfo table provides the connections back from each UserSourceReturnFlow combo
	// to the SourceDrainageFlows and DestinationDrainageFlows tables
	headers(0) = "UsrSrcRtn_ID";
	headers(1) = "User_ID";
	headers(2) = "UsersType";
	headers(3) = "SrcDrainID";
	headers(4) = "DestDrainID";
	headers(5) = "WWTP_ID";
	ncols = 6;
    filenm = dirname + "/" + "DrainageInfo.txt";
	nrows = StaticOutput.DrainageInfoSize; // size(StaticOutput%DrainageInfo)
	integer_array.resize(nrows, ncols);
	for (n = 0; n < nrows; n++) {
		integer_array(n,0) = StaticOutput.DrainageInfo[n].UserSrceReturn_ID;
		integer_array(n,1) = StaticOutput.DrainageInfo[n].User_ID;
		integer_array(n,2) = StaticOutput.DrainageInfo[n].UsersType;
		integer_array(n,3) = StaticOutput.DrainageInfo[n].SourceDrainageID;
		integer_array(n,4) = StaticOutput.DrainageInfo[n].DestinationDrainageID;
		integer_array(n,5) = StaticOutput.DrainageInfo[n].WWTP_ID;
	}
    write_integers_to_text(filenm, headers, nrows, ncols);

	//this UserFlowLinks table provides the connections back from each User to the LinkFlow table
	headers(0) = "UserID";
	headers(1) = "SourceDrainage";
	headers(2) = "DestinationDrainage";
	headers(3) = "TakeLinkID";
	headers(4) = "ReturnLinkID";
	headers(5) = "UsersType";
	ncols = 6;
    filenm = dirname + "/" + "UserFlowLinks.txt";
	nrows = StaticOutput.UserFlowLinksSize;	// size(StaticOutput%UserFlowLinks)	This wasn't allocated and used in the Fortran code.
	integer_array.resize(nrows, ncols);
	for (n = 0; n < nrows; n++) {
		integer_array(n,0) = StaticOutput.UserFlowLinks[n].UserID;
		integer_array(n,1) = StaticOutput.UserFlowLinks[n].SourceDrainage;
		integer_array(n,2) = StaticOutput.UserFlowLinks[n].DestinationDrainage;
		integer_array(n,3) = StaticOutput.UserFlowLinks[n].TakeLinkID;
		integer_array(n,4) = StaticOutput.UserFlowLinks[n].ReturnLinkID;
		integer_array(n,5) = StaticOutput.UserFlowLinks[n].UsersType;
	}
    write_integers_to_text(filenm, headers, nrows, ncols);

	return 0;
}


int write_integers_to_text(const string filenm, const Array<string,Dynamic,1> &headers, const int nrows, const int ncols)
{
	//write out a data array that has  columns of integer data
	string str;
	int i, j;
	ofstream outFile(filenm.c_str());
	if (!outFile.is_open()) {
			cerr << "Failed to open '" << filenm << "'\n";
			exit(EXIT_FAILURE);
	} else {
		cout << "Writing to '" << filenm << "'\n";
	}
	for (j = 0; j < ncols; j++) {
		outFile << headers(j) << " ";
	}
	outFile << '\n';
	for (i = 0; i < nrows; i++) {
		for (j = 0; j < ncols; j++) {
			outFile << " " << dec << setw(9) << integer_array(i,j);
		}
		outFile << '\n';
	}
	outFile.close();

	return 0;
}


using namespace input_structures;
using namespace TimeVaryingOutput;

int Write_TimeVaryingOutput_Tables(const string dirname, const int NumUser, const int NumTimesteps,
	const int NumNode, const int NumLink, const int NumReturnFlow, const int NumWWTP)
{
	// five tables to write out with real-valued time series:
	// TotalRunoff_cms,Baseflow_cms,FlowInLinks_cms,FlowAtStreamNodes_cms,ReservoirStorage_m3
	//  and one with integer-valued time series: DateTime_yyyymmdd_hhmmss
	int i, irf_ind, i_src, iuslink, icall, it, irf;
	int j, jr, k, id, j_ret, jrf, js, j_src, jsrc;
	int n, nusr, nrf, nrows, nv, ncols, nsrc, nfound, *ifound;
	double sumreturnflow, rf_to_this_dest, proportion_flow_from_this_source, flow;
	string fmtstr;
	Array<string,Dynamic,1> headers;

	string filenm, LocationTypeString;
	ostringstream location;
	ArrayXXd WWTP_Outflows_cms = ArrayXXd::Zero(NumTimesteps*NumWWTP, 3);
	double sumsourceflow;
	int usernode, k0;

	// 3 tables to create: 1 with WWTP flows, 1 with all flows from sources to users, 1 with all flows from users to destinations
	// build WWTP times series from FlowInLinks
	// for every returnflow of every user, if there's a WWTP on the rf path, add this flow
	//WWTP_Outflows_cms = 0.0;
	for (i = 1; i <= NumUser; i++) {
		irf = User[i-1].ReturnFlowID;
		if (irf > 0) {
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
			usernode = ifound[0];
			delete [] ifound;
			// find1()
			nfound = 0; //none found
			ifound = new int[NumReturnFlow];
			ifound[nfound] = 0;
			for (n = 0; n < NumReturnFlow; n++) {
				if (ReturnFlow[n].ReturnFlowID == irf) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			irf_ind = ifound[0];
			delete [] ifound;

			for (j = 1; j <= ReturnFlow[irf_ind-1].NumReturnFlows; j++) {
				if (ReturnFlow[irf_ind-1].WWTP_ID[j-1] > 0) {
					id = ReturnFlow[irf_ind-1].WWTP_ID[j-1];
					// find3()
					nfound = 0; //none found
					ifound = new int[NumLink];
					ifound[nfound] = 0;
					for (n = 0; n < NumLink; n++) {
						//cout << Link(n).LinkCode << " " << ReturnFlowLinkCode << " && ";
						//cout << Link(n).USNode << " " << usernode << " && ";
						//cout << Link(n).ReturnFlowID << " " << j << '\n';
						if (Link(n).LinkCode == ReturnFlowLinkCode && Link(n).USNode == usernode && Link(n).ReturnFlowID == j) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}
					j_ret = ifound[0]; //problem if nfound<>1
					delete [] ifound;

					// find1()
					nfound = 0; //none found
					ifound = new int[NumWWTP];
					ifound[nfound] = 0;
					for (n = 0; n < NumWWTP; n++) {
						if (WWTP_list[n] == id) {
							nfound++;
							ifound[nfound-1] = n+1;
						}
					}

					for (k = 1; k <= NumTimesteps; k++) {
						k0 = (ifound[0]-1)*NumTimesteps + k;
						WWTP_Outflows_cms(k0-1,0) = k;
						WWTP_Outflows_cms(k0-1,1) = id;
						if (j_ret >= 1) {
                            WWTP_Outflows_cms(k0-1,2) += FlowInLinks_cms(k-1,j_ret-1);
						} else {
						    WWTP_Outflows_cms(k0-1,2) = 0.0;
						    //cerr << "write_output_tables: k = " << k << " no link found\n";
						}
					}
					delete [] ifound;
				}
			}
		}
	}

	// build DrainageFlows time series from FlowInLinks
	// for every path from a source to a user to a returnflow destination
	// each return flow amount is comprised of water from perhaps more than one source
	// estimate the proportion of that retrun flow which comes from each source, by using the proportion of water taken from each source
	k = 0;
	nusr = StaticOutput.DrainageInfoSize; // size(StaticOutput%DrainageInfo,1)
	ArrayXXd SourceDrainageFlows(nusr*NumTimesteps, 3);
	SourceDrainageFlows = 0.0;
	for (i = 1; i <= NumUser; i++) {
		// find2() find user node
		nfound = 0; //none found
		ifound = new int[NumNode];
		ifound[nfound] = 0;
		for (n = 0; n < NumNode; n++) {
			if (Node[n].Type == UserNodeCode && Node[n].SelfID == i) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		usernode = ifound[0];
		delete [] ifound;

		// find2() find all return flows links from this user
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).LinkCode == ReturnFlowLinkCode && Link(n).USNode == usernode) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		Array<int,Dynamic,1> irf_list(nfound);
		for (n = 0; n < nfound; n++) {
			irf_list(n) = ifound[n];
		}
		nrf = nfound;    //RAW 4-Jul-2005 changed ifound(1) to ifound(1:nfound) - will only matter if a user has >1 return flow
		delete [] ifound;

		// find2()  //find all return flows links from this user
		nfound = 0; //none found
		ifound = new int[NumLink];
		ifound[nfound] = 0;
		for (n = 0; n < NumLink; n++) {
			if (Link(n).LinkCode == UserAbstractionLinkCode && Link(n).DSNode == usernode) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		Array<int,Dynamic,1> isrc_list(nfound);
		for (n = 0; n < nfound; n++) {
			isrc_list(n) = ifound[n];
		}
		nsrc = nfound;
		delete [] ifound;

		// only one of these return flows goes to the destination we are interested in
		for (js = 1; js <= User[i-1].NumSources; js++) {
			j_src = User[i-1].SourceID[js-1];
			// find2()  //find the source node
			nfound = 0; //none found
			ifound = new int[NumNode];
			ifound[nfound] = 0;
			for (n = 0; n < NumNode; n++) {
				if (Node[n].Type == Source[j_src-1].Type && Node[n].SelfID == Source[j_src-1].SourceLocationID) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			i_src = ifound[0];
			delete [] ifound;

			// RAW 4-Jul-2005 changed next line from find2 to find3 to include "Link%DSNode,usernode," in our search
			// before this fix we were getting the wrong iuslink

			// find3() //find the link from that node to the user
			nfound = 0; //none found
			ifound = new int[NumLink];
			ifound[nfound] = 0;
			for (n = 0; n < NumLink; n++) {
				if (Link(n).LinkCode == UserAbstractionLinkCode && Link(n).USNode == i_src && Link(n).DSNode == usernode) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			iuslink = ifound[0]; //problem if nfound<>1
			delete [] ifound;

			j_ret=User[i-1].ReturnFlowID;
			if (j_ret > 0) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumReturnFlow];
				ifound[nfound] = 0;
				for (n = 0; n < NumReturnFlow; n++) {
					if (ReturnFlow[n].ReturnFlowID == j_ret) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				irf_ind = ifound[0];	//problem if nfound<>1
				delete [] ifound;

				for (jr = 1; jr <= ReturnFlow[irf_ind-1].NumReturnFlows; jr++) {
					k++;
					for (it = 1; it <= NumTimesteps; it++) {
						sumreturnflow = 0;
						for (jrf = 1; jrf <= nrf; jrf++) {
							if (jrf == StaticOutput.DrainageInfo[k-1].RF_counter)
								rf_to_this_dest = FlowInLinks_cms(it-1,irf_list(jrf-1)-1);
							sumreturnflow += FlowInLinks_cms(it-1,irf_list(jrf-1)-1);
						} //jrf
						sumsourceflow = 0;
						for (jsrc = 1; jsrc <= nsrc; jsrc++) {
							sumsourceflow += FlowInLinks_cms(it-1,isrc_list(jsrc-1)-1);
						} //jsrc
						if (sumsourceflow > 0) {
							proportion_flow_from_this_source = FlowInLinks_cms(it-1,iuslink-1)/sumsourceflow;
						} else {
							proportion_flow_from_this_source = 1.0;
						}
						if (sumreturnflow > 0) {
							// frac_of_all_rf_via_this_rfpath=rf_to_this_dest/sumreturnflow
							// flow=flow_from_this_source*frac_of_all_rf_via_this_rfpath
							flow = rf_to_this_dest*proportion_flow_from_this_source;
						} else {
							flow = 0.0;
						}
						k0 = (k-1)*NumTimesteps + it;
						SourceDrainageFlows(k0-1,0) = it;
						SourceDrainageFlows(k0-1,1) = k;
						SourceDrainageFlows(k0-1,2) = flow;
					} //it
				} //jr
			} //j_ret
		} //js
	} //i

	icall = 1;

	nv = DateTime_yyyymmdd_hhmmss.cols();//.extent(secondDim);//second dimension is location
	nrows = DateTime_yyyymmdd_hhmmss.rows();//.extent(firstDim);
	LocationTypeString = "Drainage";
	filenm = dirname + "/" + "DateTime_yyyymmdd_hhmmss.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	headers(1) = "yyyymmdd";
	headers(2) = "hhmmss";
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
    	for (j = 1; j <= nrows; j++) {
			integer_array(j-1,i) = DateTime_yyyymmdd_hhmmss(j-1,i-1);
    	}
    }
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	icall++;
	// 	nv=size(TotalRunoff_cms,2) !second dimension is location
	// 	nrows=size(TotalRunoff_cms,1)
	// 	LocationTypeString="Drainage"
	// 	filenm=dirname(1:len_trim(dirname))  // '\' // 'TotalRunoff_cms' // '.txt';
	// 	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	// 	headers(1)='TimeStep'
	// 	lenLTS=len(trim(LocationTypeString))
	// 	do i=1,nv
	// 		n=1+int(log10(real(i)))
	// 		write(fmtstr,'(''(A,i'',i1,'')'')')n
	// 		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	// 	end do
	// 	ncols=1+nv
	// 	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	// 	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	// 	do j=1,nrows
	// 		integer_array(j,1)=j
	// 	end do
	//    do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=TotalRunoff_cms(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)

	//  Artificial drainage DGT 6/28/05
	//	icall=icall+1
	//	nv=size(Artificial_Drainage,2) !second dimension is location
	//	nrows=size(Artificial_Drainage,1)
	//	LocationTypeString='Drainage'
	//	filenm=dirname(1:len_trim(dirname))  // '\' // 'Artificial_Drainage' // '.txt';
	//	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	//	headers(1)='TimeStep'
	//	lenLTS=len(trim(LocationTypeString))
	//	do i=1,nv
	//		n=1+int(log10(real(i)))
	//		write(fmtstr,'(''(A,i'',i1,'')'')')n
	//		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	//	end do
	//	ncols=1+nv
	//	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	//	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	//	do j=1,nrows
	//		integer_array(j,1)=j
	//	end do
	//   do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=Artificial_Drainage(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)


	icall++;
	//	nv=size(Baseflow_cms,2) !second dimension is location
	//	nrows=size(Baseflow_cms,1)
	//	LocationTypeString='Drainage'
	//	filenm=dirname(1:len_trim(dirname))  // '\' // 'Baseflow_cms' // '.txt';
	//	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	//	headers(1)='TimeStep'
	//	lenLTS=len(trim(LocationTypeString))
	//	do i=1,nv
	//		n=1+int(log10(real(i)))
	//		write(fmtstr,'(''(A,i'',i1,'')'')')n
	//		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	//	end do
	//	ncols=1+nv
	//	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	//	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	//	do j=1,nrows
	//		integer_array(j,1)=j
	//	end do
 	//   do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=Baseflow_cms(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)

	icall++;
	nv = FlowInLinks_cms.cols();//.extent(secondDim); 	//second dimension is location
	nrows = FlowInLinks_cms.rows();//.extent(firstDim);
	LocationTypeString = "LinkID";
	filenm = dirname + "/" + "FlowInLinks_cms.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	for (i = 1; i <= nv; i++) {
		location << dec << setw(3) << i;
		headers(i) = LocationTypeString + location.str();
		location.clear();
		location.str("");
	}
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
		for (j = 1; j <= nrows; j++) {
			real_array(j-1,i) = FlowInLinks_cms(j-1,i-1);
		}
	};
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	icall++;
	//	nv=size(FlowAtStreamNodes_cms,2) !second dimension is location
	//	nrows=size(FlowAtStreamNodes_cms,1)
	//	LocationTypeString='Node'
	//	filenm=dirname(1:len_trim(dirname))  // '\' // 'FlowAtStreamNodes_cms' // '.txt';
	//	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	//	headers(1)='TimeStep'
	//	lenLTS=len(trim(LocationTypeString))
	//	do i=1,nv
	//		n=1+int(log10(real(i)))
	//		write(fmtstr,'(''(A,i'',i1,'')'')')n
	//		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	//	end do
	//	ncols=1+nv
	//	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	//	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	//	do j=1,nrows
	//		integer_array(j,1)=j
	//	end do
	//    do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=FlowAtStreamNodes_cms(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)

	icall++;
	nv = ReservoirStorage_m3.cols();//.extent(secondDim);	//second dimension is location
	nrows = ReservoirStorage_m3.rows();//.extent(firstDim);
	LocationTypeString = "Reservoir";
	filenm = dirname + "/" + "ReservoirStorage_m3.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	for (i = 1; i <= nv; i++) {
		location << dec << setw(3) << i;
		headers(i) = LocationTypeString + location.str();
		location.clear();
		location.str("");
	}
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
		for (j = 1; j <= nrows; j++) {
			real_array(j-1,i) = ReservoirStorage_m3(j-1,i-1);
		}
	};
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	icall++;
	// DGT 6/29/05  Want to handle precipitation and evaporation on a line by line basis to be more efficient about memory
	//	nv=size(Precipitation_mm,2) !second dimension is location
	//	nrows=size(Precipitation_mm,1)
	//	LocationTypeString='Drainage'
	//	filenm=dirname(1:len_trim(dirname))  // '\' // 'Precipitation_mm' // '.txt';
	//	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	//	headers(1)='TimeStep'
	//	lenLTS=len(trim(LocationTypeString))
	//	do i=1,nv
	//		n=1+int(log10(real(i)))
	//		write(fmtstr,'(''(A,i'',i1,'')'')')n
	//		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	//	end do
	//	ncols=1+nv
	//	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	//	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	//	do j=1,nrows
	//		integer_array(j,1)=j
	//	end do
 	//   do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=Precipitation_mm(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)

	icall++;
	//	nv=size(Evaporation_mm,2) !second dimension is location
	//	nrows=size(Evaporation_mm,1)
	//	LocationTypeString='Drainage'
	//	filenm=dirname(1:len_trim(dirname))  // '\' // 'Evaporation_mm' // '.txt';
	//	if (allocated(headers)) deallocate (headers); allocate(headers(nv+1))
	//	headers(1)='TimeStep'
	//	lenLTS=len(trim(LocationTypeString))
	//	do i=1,nv
	//		n=1+int(log10(real(i)))
	//		write(fmtstr,'(''(A,i'',i1,'')'')')n
	//		write(headers(i+1),fmt=fmtstr)LocationTypeString(1:lenLTS),i
	//	end do
	//	ncols=1+nv
	//	if (allocated(integer_array)) deallocate (integer_array); allocate(integer_array(nrows,ncols))
	//	if (allocated(real_array)) deallocate (real_array); allocate(real_array(nrows,ncols))
	//	do j=1,nrows
	//		integer_array(j,1)=j
	//	end do
 	//   do i=1,nv;
	//		do j=1,nrows; real_array(j,i+1)=Evaporation_mm(j,i); end do
	//	end do;
	//	call write_struct_to_text(filenm,headers,nrows,ncols,icall)

	icall++;
	nv = NumUser; 	//second dimension is user
	nrows = NumTimesteps;
	LocationTypeString = "       User";
	filenm = dirname + "/" + "UserDemand_cms.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	for (i = 1; i <= nv; i++) {
		location << dec << setw(3) << i;
		headers(i) = LocationTypeString + location.str();
		location.clear();
		location.str("");
	}
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
		for (j = 1; j <= nrows; j++) {
			real_array(j-1,i) = User[i-1].Demand[j-1]/86400.0;
		}
	};
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	icall++;
	nv = NumUser; 	//second dimension is user
	nrows = NumTimesteps;
	LocationTypeString = "       User";
	filenm = dirname + "/" + "UserDeficit_cms.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	for (i = 1; i <= nv; i++) {
		location << dec << setw(3) << i;
		headers(i) = LocationTypeString + location.str();
		location.clear();
		location.str("");
	}
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
		for (j = 1; j <= nrows; j++) {
			real_array(j-1,i) = User[i-1].Deficit[j-1]/86400.0;
		}
	};
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	icall++;
	nv = NumUser; 	//second dimension is user
	nrows = NumTimesteps;
	LocationTypeString = "       User";
	filenm = dirname + "/" + "UserWithdrawal_cms.txt";
	headers.resize(nv+1);
	headers(0) = "TimeStep";
	for (i = 1; i <= nv; i++) {
		location << dec << setw(3) << i;
		headers(i) = LocationTypeString + location.str();
		location.clear();
		location.str("");
	}
	ncols = 1 + nv;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = j;
	}
    for (i = 1; i <= nv; i++) {
		for (j = 1; j <= nrows; j++) {
			real_array(j-1,i) = User[i-1].Withdrawal[j-1]/86400.0;
		}
	};
	write_struct_to_text(filenm, headers, nrows, ncols, icall);


	icall = -1;
	nv = WWTP_Outflows_cms.cols();//.extent(secondDim); 	//second dimension is NOT location
	nrows = WWTP_Outflows_cms.rows();//.extent(firstDim);
	filenm = dirname + "/" + "WWTP_Outflows_cms.txt";
	headers.resize(3);
	headers(0) = "TimeStep";
	headers(1) = "WWTP_ID";
	headers(2) = "Flow_cms";
	ncols = 3;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = WWTP_Outflows_cms(j-1,0);
		integer_array(j-1,1) = WWTP_Outflows_cms(j-1,1);
	}
	for (j = 1; j <= nrows; j++) {
		real_array(j-1,2) = WWTP_Outflows_cms(j-1,2);
	}
	write_struct_to_text(filenm, headers, nrows, ncols, icall);


	icall = -1;
	nv = SourceDrainageFlows.cols();//.extent(secondDim); 	//second dimension is NOT location
	nrows = SourceDrainageFlows.rows();//.extent(firstDim);
	filenm = dirname + "/" + "SourceDrainageFlows.txt";
	headers.resize(3);
	headers(0) = "TimeStep";
	headers(1) = "UsrSrcRtn_ID";
	headers(2) = " Flow_cms";
	ncols = 3;
	integer_array.resize(nrows, ncols);
	real_array.resize(nrows, ncols);
	for (j = 1; j <= nrows; j++) {
		integer_array(j-1,0) = SourceDrainageFlows(j-1,0);
		integer_array(j-1,1) = SourceDrainageFlows(j-1,1);
	}
	for (j = 1; j <= nrows; j++) {
		real_array(j-1,2) = SourceDrainageFlows(j-1,2);
	}
	write_struct_to_text(filenm, headers, nrows, ncols, icall);

	return 0;
}


int write_struct_to_text(const string filenm, const Array<string,Dynamic,1> &headers, const int nrows, const int ncols,
	const int icall)
{
	// write out a data array that has 1 columns of integer data, and then columns of real data
	string str;
	string realfmt;

	int i, j;
	ofstream outFile(filenm.c_str());
	cout << "Writing to '" << filenm << "'\n";

	for (j = 0; j < ncols; j++) {
		outFile << setw(18) << headers(j);
	}
	outFile << '\n';
	if (icall > 1) {
		for (i = 0; i < nrows; i++) {
			outFile << dec << setw(18) << integer_array(i,0);
			for (j = 1; j < ncols; j++) {
				outFile << fixed << setw(18) << setprecision(7) << real_array(i,j);
			}
			outFile << '\n';
		}
	} else if (icall == 1) {
		for (i = 0; i < nrows; i++) {
			outFile << dec << setw(18) << integer_array(i,0);
			for (j = 1; j < ncols; j++) {
				outFile << dec << setw(18) << integer_array(i,j);
			}
			outFile << '\n';
		}
	} else if (icall < 0) {
		for (i = 0; i < nrows; i++) {
			outFile << dec << setw(12) << integer_array(i,0);
			outFile << dec << setw(12) << integer_array(i,1);
			//for (j = 1; j < ncols; j++) {
				outFile << scientific << setw(18) << setprecision(7) << real_array(i,2);
			//}
			outFile << '\n';

		}
	} else {
		cerr << "illegal call to write_struct_to_text with icall= " << icall << '\n';
	}
	outFile.close();

return 0;
}


