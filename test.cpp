#include <sstream>

#include "EnergyLoss.h"

int main() {
    printf("Test and example of EnergyLoss Library\n\n");

    EnergyLoss *carbon = new EnergyLoss("CarbonCsI.dat", false);
    EnergyLoss *carbonUseFactor = new EnergyLoss("CarbonCsIReduced.dat", 4.5099E+02);

    // CalcRemainder
    printf("Calculating energy remaining of 12C through CsI:\n");
    printf("\t 150 MeV 12C through 0.1 mm CsI. Remaining Energy: %12.9f MeV\n", carbon->CalcRemainder(150., .1));
    printf("\t\t LISE++ Value: 118.1100 MeV\n");
    printf("\t Using carbonUseFactor. 150 MeV 12C through 0.1 mm CsI. Remaining Energy: %12.9f MeV\n", carbonUseFactor->CalcRemainder(150., .1));
    printf("\t\t LISE++ Value: 118.1100 MeV\n");
    printf("\t 72 MeV 12C through 0.1 mm CsI. Remaining Energy: %12.9f MeV\n", carbon->CalcRemainder(72., .1));
    printf("\t\t LISE++ Value: 3.2985 MeV\n");
    printf("\t 15 MeV 12C through 0.005 mm CsI. Remaining Energy: %12.9f MeV\n", carbon->CalcRemainder(15., .005));
    printf("\t\t LISE++ Value: 9.8093 MeV\n");
    printf("\t 1 MeV 12C through 0.001 mm CsI. Remaining Energy: %12.9f MeV\n", carbon->CalcRemainder(1., 0.001));
    printf("\t\t LISE++ Value: 0.3027 MeV\n");

    //AddBack
    printf("\nCalculating initial energy of 12C when it goes through CsI:\n");
    printf("\t 0.5 MeV 12C after 0.001 mm CsI. Initial Energy: %8.4f MeV\n", carbon->AddBack(0.5, 0.001));
    printf("\t\t LISE++ Value: 1.33 MeV\n");
    printf("\t 1 MeV 12C after 0.1 mm CsI. Initial Energy: %8.4f MeV\n", carbon->AddBack(1., 0.1));
    printf("\t\t LISE++ Value: 70.99 MeV\n");

    // CalcRange
    printf("\nCalculating the range 12C can go with initial and final energy given:\n");
    printf("\t 12C 10 MeV to 0.5 MeV in CsI. Range: %8.4f mm\n", carbon->CalcRange(10., 0.5));
    printf("\t\t LISE++ Value: 0.0088 mm\n");
    printf("\t 12C 50 MeV to 0. MeV in CsI. Range: %8.4f mm\n", carbon->CalcRange(50., 0.));
    printf("\t\t LISE++ Value: 0.0624 mm\n");
    printf("\t 12C 100 MeV to 0. MeV in CsI. Range: %8.4f mm\n", carbon->CalcRange(100., 0.));
    printf("\t\t LISE++ Value: 0.1699 mm\n");

    // Test Precision
    printf("\nTesting different number of GL Points for 1000 MeV 12C through 7.75 mm CsI:\n");
    carbon->UseGL16();
    printf("\t16-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL32();
    printf("\t32-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL64();
    printf("\t64-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL128();
    printf("\t128-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL256();
    printf("\t256-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL512();
    printf("\t512-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    carbon->UseGL1024();
    printf("\t1024-point GL. Remaining Energy: %18.16g MeV\n", carbon->CalcRemainder(1000., 7.75));
    printf("\t\t LISE++ Value: 0 MeV\n");

    delete carbon;
    delete carbonUseFactor;

    printf("\n Testing read of basic dEdx file\n\n");

    EnergyLoss* proton = new EnergyLoss();
    proton->ReadBasicdEdx("ProtonD2.dat");

    printf("Calculating energy remaining of Proton through solid D2 target:\n");
    printf("\t 10 MeV Proton through 1.0 mm of solid D2. Remaining Energy: %12.9f MeV\n", proton->CalcRemainder(10., 1.0));
    delete proton;

    return 0;
}
