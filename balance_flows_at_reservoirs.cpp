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

int remove_element(Array<int,Dynamic,1> &j_out, int &nj_out, const int iFound);

int BalanceFlowsAtReservoirs(const int NumNode, const int NumLink, const int NumUser, const int NumReservoir, ArrayXd &ReservoirNetStorage)
{
	Array<int,Dynamic,1> j_in(1000);
	int m, n, isn, nj_in, nj_out, i, k, j_release, j_found, iuser, jret;
	int nfound;// *ifound;
	valarray<int> ifound;
	double flow_in, flow_taken, flow_released;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> BalanceFlowsAtReservoirs(" << ncalls << ")" << std::endl;
    }
	caller = "BalanceFlowsAtReservoirs";
#endif

	//Calculate storage in all reservoir nodes, and if there is overflow, update
	//flows in outflow link as required

	for (i = 1; i <= NumReservoir; i++) {
		//    call find2(Node%Type,ReservoirNodeCode,Node%SelfID,i,NumNode);
		//	isn=ifound(1) !problem if nfound<>1
		//    call find1(Link%DSNode,isn,NumLink); allocate (j_in(nfound));
		//	j_in=ifound; nj_in=nfound
		//    call find1(Link%USNode,isn,NumLink); allocate (j_out(nfound));
		//	j_out=ifound; nj_out=nfound

		isn = Reservoir[i-1].isn;
		nj_in = Reservoir[i-1].n_in;
		for (n = 0; n < nj_in; n++) {
			j_in(n) = Reservoir[i-1].ifound_in(n);
		}

		nj_out = Reservoir[i-1].n_taken;
		Array<int,Dynamic,1> j_out(nj_out);
		for (n = 0; n < nj_out; n++) {
			j_out(n) = Reservoir[i-1].ifound_taken(n);
		}

		//find the link going to the user node of Type=ReservoirReleaseUseCode that takes reservoir overflows from this reservoir
		ifound.resize(NumUser*nj_out,0);   // initialized to 0
		if (Reservoir[i-1].InOffStream == InStreamReservoirCode) {
			// find1()
			nfound = 0; //none found
			for (n = 0; n < NumUser; n++) {
				for (m = 0; m < nj_out; m++) {
					if (User[Node[Link(j_out(m)-1).DSNode-1].SelfID-1].UsersType == InStreamReservoirReleaseUseCode) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
			}
		} else if (Reservoir[i-1].InOffStream == OffStreamReservoirCode) {
			// find1()
			nfound = 0; //none found
			for (n = 0; n < NumUser; n++) {
				for (m = 0; m < nj_out; m++) {
					if (User[Node[Link(j_out(m)-1).DSNode-1].SelfID-1].UsersType == OffStreamReservoirReleaseUseCode) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
			}
		}
		j_release = j_out(ifound[0]-1); //problem if nfound<>1

		// find1()
		nfound = 0; //none found
		ifound.resize(nj_out,0);   // initialized to 0
		for (n = 0; n < nj_out; n++) {
			if (j_out(n) == j_release) {
				nfound++;
				ifound[nfound-1] = n+1;
			}
		}
		j_found = ifound[0];	//problem if nfound<>1
		remove_element(j_out, nj_out, j_found);
		flow_in = 0.0;
		for (k = 1; k <= nj_in; k++) {
			flow_in += Link(j_in(k-1)-1).Flow;
		}
		flow_taken = 0.0;
		for (k = 1; k <= nj_out; k++) {
			flow_taken += Link(j_out(k-1)-1).Flow;
		}
		flow_released = Link(j_release-1).Flow;

		Node[isn-1].Store = Node[isn-1].StoreOld + flow_in - flow_taken - flow_released;
		if (Node[isn-1].Store > Reservoir[i-1].StoreMax) {
			Link(j_release-1).Flow += Node[isn-1].Store-Reservoir[i-1].StoreMax;
			Node[isn-1].Store = Reservoir[i-1].StoreMax;
			iuser = Link(j_release-1).DSNode;  // This is the user identifier of the dummy user that effects reservoir release
			// Need to find the return flow link associated with this user so that the corresponding return flow can also be
			// adjusted by the additional flow in this link.
			// find2()
			nfound = 0; //none found
			ifound.resize(NumLink,0);   // initialized to 0
			for (n = 0; n < NumLink; n++) {
				if (Link(n).LinkCode == ReturnFlowLinkCode && Link(n).USNode == iuser) {
					nfound++;
					ifound[nfound-1] = n+1;
				}
			}
			jret = ifound[0]; //problem if nfound<>1
			Link(jret-1).Flow += Link(j_release-1).Flow; //no consumption at this usernode

			// Logically the line above would be better implemented by
			//  Link(jret)%Flow=Link(j_release)%Flow
			//  Did not do this because DGT did not want to change code and Ross says will not change result.
		}
		ReservoirNetStorage(i-1) = Node[isn-1].Store-Reservoir[i-1].StoreMin;
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving BalanceFlowsAtReservoirs(" << ncalls << ")" << endl;
    }
    ncalls++;
#endif

	return 0;
}

int remove_element(Array<int,Dynamic,1> &j_out, int &nj_out, const int iFound)
{
	int i;

	for (i = iFound; i <= nj_out-1; i++) {
		j_out(i-1) = j_out(i);
	}
	nj_out--;

	return 0;
}

