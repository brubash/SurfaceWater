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

int CalculateDemand(const int ThisMonth, const int doy, const ArrayXd &vol_irrig_demand, const int NumUser, const int NumMonthlyDemand)
{
	int i, n, idt, nfound, *ifound;
	double DDF = 0.0;

#if TRACE
	static int ncalls = 0;
	string save_caller = caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " -> CalculateDemand(" << ncalls << ")" << std::endl;
    }
	caller = "CalculateDemand";
#endif
	for (i = 1; i <= NumUser; i++) {
		//  Debugging code to stop at a specific user
		//	if (i .eq. 616)then
		//	  DDF=DDF
		//	endif
		if (User[i-1].UsersType == SoilMoistureIrrigationUseCode) {
			// DGT 5/18/05.  SoilMoistureIrrigationUseCode is irrigation demand from soil moisture
			//			User(i)%DemandVble=vol_irrig_demand(User(i)%POU_ID) // DGT 8/20/05 Reinterpret vol_irrig_demand
			//           as depth of irrig demand.  Therefore use it as the rate variable.
			//           The demand variable is the irrigated area in square meters from user.txt
			User[i-1].DemandRate = vol_irrig_demand(User[i-1].POU_ID-1); // 1.0
			DDF = 1.0;
		} else {
			//		if (User(i)%UsersType.ne.-1) then
			if (User[i-1].InYearDemandType == 0) {
				DDF = 1.0;
			} else if (User[i-1].InYearDemandType > 0) {
				// find1()
				nfound = 0; //none found
				ifound = new int[NumMonthlyDemand];
				ifound[nfound] = 0;
				for (n = 0; n < NumMonthlyDemand; n++) {
					if (MonthlyDemand[n].InYearDemandType == User[i-1].InYearDemandType) {
						nfound++;
						ifound[nfound-1] = n+1;
					}
				}
				idt = ifound[0];	//problem if nfound<>1
				delete [] ifound;
				//               idt=User(i).InYearDemandType  ! DGT 8/20/05 relying on InYearDemandType
				// in MonthlyDemandFraction.txt being a counting sequence starting at 1
				DDF = MonthlyDemand[idt-1].Month[ThisMonth-1];
			} else {
	            std::cout << "Need to implement daily demand fractions for every day of year.\n";
				//call find1(DailyDemand.InYearDemandType,User(i).InYearDemandType, NumDailyDemand);
				//idt=ifound(1) !problem if nfound<>1
				//DDF=DailyDemand(idt).Day(DOY);
			}
		}

		User[i-1].DemandToday = User[i-1].DemandVble*User[i-1].DemandRate*DDF;
	}
#if TRACE
	caller = save_caller;
	if (ncalls < MAX_TRACE) {
        traceFile << setw(30) << caller << " <- Leaving CalculateDemand(" << ncalls << ")" << "\n\n";
    }
    ncalls++;
#endif

	return 0;
}
