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
#include <vector>

using namespace Eigen;

//int *ifound;
std::vector<int> ifound;            // record of all individual finds
std::vector<int> find_counts;       // record of find counts
std::vector<int> done;              // tally of finished assignments
int nfound, nfound_last;
int global_found = 0;

void Find(Array<DrainageType,Dynamic,1> &iarray1, const int ival1, const int num)
{
    int i;
    nfound = 0;
    for (i = 1; i <= num; i++) {
        if (iarray1(i-1).DSDrainage == ival1) {
            nfound++;
            global_found++;
            ifound.push_back(i);
        }
    }
    find_counts.push_back(nfound);
}

/*	subroutine find1(iarray1,ival1,num)
	use findmodule
	integer num,ival1,i
	integer iarray1(num)

	if (allocated(ifound)) deallocate (ifound)
	allocate (ifound(num))

	nfound=0 !none found
	do i =1,num
		if (iarray1(i).eq.ival1) then
			nfound=nfound+1
			ifound(nfound)=i
		endif
	end do
	end

	subroutine find2(iarray1,ival1,iarray2,ival2,num)
! this returns in ifound the list of indices for which iarray1(i)=ival1 and
! iarray2(i)=ival2.  The search is from 1 to num
	use findmodule
	integer num,iarray1(num),ival1,iarray2(num),ival2,i

	if (allocated(ifound)) deallocate (ifound)
	allocate (ifound(num))

	nfound=0 !none found
	do i =1,num
		if (iarray1(i).eq.ival1 .and. iarray2(i).eq.ival2) then
			nfound=nfound+1
			ifound(nfound)=i
		endif
	end do
	end

	subroutine find3(iarray1,ival1,iarray2,ival2,iarray3,ival3,num)
	use findmodule
	integer num,iarray1(num),ival1,iarray2(num),ival2,iarray3(num),ival3,i
	if (allocated(ifound)) deallocate (ifound)
	allocate (ifound(num))

	!A & B & C
	nfound=0 !none found
	do i =1,num
		if ((iarray1(i).eq.ival1) .and. (iarray2(i).eq.ival2) .and. (iarray3(i).eq.ival3)) then
			nfound=nfound+1
			ifound(nfound)=i
		endif
	end do
	end

	subroutine find2a(iarray1,ival1,iarray2,ival2,num)
	use findmodule
	integer num,iarray1(num),ival1,iarray2(num),ival2,iarray3(num),ival3,i
	if (allocated(ifound)) deallocate (ifound)
	allocate (ifound(num))

	!A & (~B)
	nfound=0 !none found
	do i =1,num
		if ((iarray1(i).eq.ival1) .and. (iarray2(i).ne.ival2)) then
			nfound=nfound+1
			ifound(nfound)=i
		endif
	end do
	end

	subroutine find3b(iarray1,ival1,iarray2,ival2,iarray3,ival3,num)
	use findmodule
	integer num,iarray1(num),ival1,iarray2(num),ival2,iarray3(num),ival3,i
	if (allocated(ifound)) deallocate (ifound)
	allocate (ifound(num))

	!A & (B | C)
	nfound=0 !none found
	do i =1,num
		if ((iarray1(i).eq.ival1) .and. (iarray2(i).eq.ival2) .or. (iarray3(i).eq.ival3)) then
			nfound=nfound+1
			ifound(nfound)=i
		endif
	end do
	end
*/
