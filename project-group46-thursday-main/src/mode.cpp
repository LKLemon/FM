#include "mode.h"

parameter set_parameter(int mode) {
      parameter param;
      switch(mode) {
        case 0:
          param.RF = 2400e3;
          param.IF = 240e3;
          param.AF = 48e3;
          param.SPS = 19;
          break;
        case 1:
          param.RF = 960e3;
          param.IF = 240e3;
          param.AF = 48e3;
          param.SPS = 0;
          break;
        case 2:
          param.RF = 2400e3;
          param.IF = 240e3;
          param.AF = 44.1e3;
          param.SPS = 42;
          break;
        case 3:
          param.RF = 1440e3;
          param.IF = 360e3;
          param.AF = 44.1e3;
          param.SPS = 0;
          break;
      }
      param.U = param.AF / std::__gcd(param.IF,param.AF);
      param.D0 = param.RF / param.IF;
      param.D1 = param.IF * param.U / param.AF;
      param.GCD_RDS = std::__gcd(param.IF,param.SPS*2375);
      param.U_RDS = param.SPS*2375 / param.GCD_RDS;
      param.D_RDS = param.IF / param.GCD_RDS;
      return param;
  }
  // int main()
  // {
  //     int mode = 0;
  //     parameter param = findDDU(mode);
  //     std::cout<< "the RF:\t\t"<< param.RF<<"\nthe IF:\t\t"<< param.IF << "\nthe Audio_Fs:\t"
  //     << param.Audio_Fs <<"\nthe U:\t\t"<< param.U <<"\nthe D0:\t\t"<< param.D0 <<"\nthe D1:\t\t"<< param.D1;
  //
  //     return 0;
  // }
//     return;
// }
