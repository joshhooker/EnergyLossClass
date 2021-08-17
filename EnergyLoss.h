#ifndef EnergyLoss_h
#define EnergyLoss_h

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iostream>
#include <string>
#include <vector>

#define ASSERT_WITH_MESSAGE_EL(condition, message)                                \
  do {                                                                         \
    if (!(condition)) {                                                        \
      printf((message));                                                       \
    }                                                                          \
    assert((condition));                                                       \
  } while (false)

class CubicSpline_EL {
public:
    CubicSpline_EL();
    template <typename T> CubicSpline_EL(const std::vector<T> &x, const std::vector<T> &y);
    template <typename T, int N, int M> CubicSpline_EL(const T (&x) [N], const T (&y) [M]);
    template <typename T, std::size_t N, std::size_t M> CubicSpline_EL(const std::array<T, N>& x, const std::array<T, M>& y);
    template <typename T> void SetPoints(const std::vector<T> &x, const std::vector<T> &y);
    template <typename T, int N, int M> void SetPoints(const T (&x) [N], const T (&y) [M]);
    template <typename T, std::size_t N, std::size_t M> void SetPoints(const std::array<T, N>& x, const std::array<T, M>& y);
    template <typename T> double operator()(T x) const;
    ~CubicSpline_EL();

private:
    size_t size;
    std::vector<double> xVec, yVec;
    std::vector<double> bVec, cVec, dVec;

    void SetSpline();
};


inline CubicSpline_EL::CubicSpline_EL() {}

template<typename T> inline CubicSpline_EL::CubicSpline_EL(const std::vector<T> &x, const std::vector<T> &y) {
    ASSERT_WITH_MESSAGE_EL(x.size() == y.size(),
                        "In CubicSpline initialization, x vector size != y vector size\n");
    assert(x.size() == y.size());
    size = x.size();
    xVec = x; yVec = y;
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

template <typename T, int N, int M> inline CubicSpline_EL::CubicSpline_EL(const T (&x) [N], const T (&y) [M]) {
    ASSERT_WITH_MESSAGE_EL(N == M,
                        "In CubicSpline initialization, x array size != y array size\n");
    assert(N == M);
    size = N;
    xVec.assign(x, x+N); yVec.assign(y, y+M);
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

template <typename T, std::size_t N, std::size_t M> inline CubicSpline_EL::CubicSpline_EL(const std::array<T, N>& x, const std::array<T, M>& y) {
    ASSERT_WITH_MESSAGE_EL(N == M,
                        "In CubicSpline initialization, x array size != y array size\n");
    size = N;
    xVec.resize(size); yVec.resize(size);
    std::copy(x.begin(), x.begin()+size, xVec.begin());
    std::copy(y.begin(), y.begin()+size, yVec.begin());
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

template <typename T> inline void CubicSpline_EL::SetPoints(const std::vector<T> &x, const std::vector<T> &y) {
    ASSERT_WITH_MESSAGE_EL(x.size() == y.size(),
                        "In CubicSpline SetPoints, x vector size != y vector size\n");
    size = x.size();
    xVec = x; yVec = y;
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

template <typename T, int N, int M> inline void CubicSpline_EL::SetPoints(const T (&x) [N], const T (&y) [M]) {
    ASSERT_WITH_MESSAGE_EL(N == M,
                        "In CubicSpline SetPoints, x array size != y array size\n");
    size = N;
    xVec.assign(x, x + N); yVec.assign(y, y + M);
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

template <typename T, std::size_t N, std::size_t M> inline void CubicSpline_EL::SetPoints(const std::array<T, N>& x, const std::array<T, M>& y) {
    ASSERT_WITH_MESSAGE_EL(N == M,
                        "In CubicSpline SetPoints, x array size != y array size\n");
    size = N;
    xVec.resize(size); yVec.resize(size);
    std::copy(x.begin(), x.begin()+size, xVec.begin());
    std::copy(y.begin(), y.begin()+size, yVec.begin());
    bVec.resize(size); cVec.resize(size); dVec.resize(size);

    SetSpline();
}

void inline CubicSpline_EL::SetSpline() {
    std::vector<double> h(size), alpha(size), l(size), z(size), u(size);

    l[0] = 1.;
    u[0] = 0.;
    z[0] = 0.;
    l[size - 1] = 1.;
    u[size - 1] = 0.;
    cVec[size - 1] = 0.;
    for(unsigned int i = 0; i < size - 1; i++) {
        ASSERT_WITH_MESSAGE_EL(xVec[i + 1] > xVec[i],
                            "In CubicSpline SetSpline, x array is not sorted from smallest to largest\n");
        assert(xVec[i + 1] > xVec[i]);
        h[i] = xVec[i + 1] - xVec[i];
        if(i > 0) {
            alpha[i] = (3./h[i])*(yVec[i + 1] - yVec[i]) - (3./h[i - 1])*(yVec[i] - yVec[i - 1]);
            l[i] = 2.*(xVec[i + 1] - xVec[i - 1]) - h[i - 1]*u[i - 1];
            u[i] = h[i]/l[i];
            z[i] = (alpha[i] - h[i - 1]*z[i - 1])/l[i];
        }
    }
    for(int i = size - 2; i > -1; i--) {
        cVec[i] = z[i] - u[i]*cVec[i + 1];
        bVec[i] = (yVec[i + 1] - yVec[i])/h[i] - h[i]*(cVec[i + 1] + 2.*cVec[i])/3.;
        dVec[i] = (cVec[i + 1] - cVec[i])/(3.*h[i]);
    }
}

template <typename T> inline double CubicSpline_EL::operator()(T x) const{
    double xs = static_cast<double>(x);

    int l = 0;
    int h = size;
    while(l < h) {
        int mid = (l + h)/2;
        if(xs <= xVec[mid]) {
            h = mid;
        } else {
            l = mid + 1;
        }
    }

    size_t idx = (l == 0) ? 0 : l - 1;

    double xi = xs-xVec[idx];
    double result;
    if(idx == 0) result = yVec[0] + bVec[0]*xi + cVec[0]*xi*xi;
    else if(idx == size - 1) result = yVec[size - 1] + bVec[size - 1]*xi + cVec[size - 1]*xi*xi;
    else result = yVec[idx] + bVec[idx]*xi + cVec[idx]*xi*xi + dVec[idx]*xi*xi*xi;
    return result;
}

inline CubicSpline_EL::~CubicSpline_EL() = default;

class EnergyLoss {
public:
    EnergyLoss();
    EnergyLoss(const char*);
    EnergyLoss(const char*, double);
    EnergyLoss(const char*, bool);
    EnergyLoss(const char*, double , bool);

    void ReadBasicdEdx(const char*);

    void ReadLISEFile(const char*);
    void ReadLISEdEdx(const char*, double);
    void ReadLISEdEdx(const char*, double, int);
    void ReadLISEdEdx(const char*, double, double);
    void ReadLISEdEdx(const char*, double, double, int);

    void SetDebug(bool);

    void ReadInitParams();
    double CalcRemainder(double, double);
    double AddBack(double, double);
    double CalcRange(double, double);
    void CalcRemainderError(double);
    void AddBackError(double);
    void AddBackHigh(double);

    friend double GetdEdx(const EnergyLoss&, double);

    void UseGL16();
    void UseGL32();
    void UseGL64();
    void UseGL128();
    void UseGL256();
    void UseGL512();
    void UseGL1024();

    void UseCalcRangeDebug();

private:
    bool debug = false;
    bool calcRangeDebug = false;

    void ReadFile(const char*);
    void ReadFileSimple(const char*, double);

    std::vector<double> energy_;
    std::vector<double> dEdx_;
    CubicSpline_EL energySpline;

    double x16[8], w16[8], x32[16], w32[16], x64[32], w64[32];
    double x128[64], w128[64], x256[128], w256[128];
    double x512[256], w512[256], x1024[512], w1024[512];

    double CalcRemainderErr;
    double AddBackErr;
    double AddBackHighPoint;

    bool GL16, GL32, GL64, GL128, GL256, GL512, GL1024;

    double CalcRangeGL16(double, double);
    double CalcRangeGL32(double, double);
    double CalcRangeGL64(double, double);
    double CalcRangeGL128(double, double);
    double CalcRangeGL256(double, double);
    double CalcRangeGL512(double, double);
    double CalcRangeGL1024(double, double);

    // Variables to read SRIM File
    const int MAX_CHARS_PER_LINE = 1024000;
    const char* const DELIMITER = " ";
    bool beforeMultiply = true;
    std::vector<char*> token;
    double stoppingConversionPower;
    std::vector<double> stoppingConversionPower_;
    std::vector<std::string> stoppingConversionUnit1_;
    std::vector<std::string> stoppingConversionUnit2_;

    // Variables to read LISE File
    std::vector<double> lise_energy_mevu_; // in MeV/u
    std::vector<double> lise_dedx_columns_[6]; // in MeV/um or MeV/mg/cm2

    std::string PrintOutput(std::string Output, std::string Color);
};

inline EnergyLoss::EnergyLoss() {
    debug = false;
}

inline EnergyLoss::EnergyLoss(const char* srimFile) {
    debug = false;
    ReadFile(srimFile);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, bool debugger) {
    debug = debugger;
    ReadFile(srimFile);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, double stoppingPower) {
    debug = false;
    ReadFileSimple(srimFile, stoppingPower);
}

inline EnergyLoss::EnergyLoss(const char* srimFile, double stoppingPower, bool debugger) {
    debug = debugger;
    ReadFileSimple(srimFile, stoppingPower);
}

inline void EnergyLoss::ReadFile(const char* srimFile) {
    std::ifstream inFile(srimFile);
    ASSERT_WITH_MESSAGE_EL(inFile.is_open(), "Cannot find SRIM file!\n");

    if(debug) std::cout << PrintOutput("DEBUGGING: READING IN SRIM FILE", "red") << std::endl;

    while(!inFile.eof()) {
        char buf[1024000];
        inFile.getline(buf, MAX_CHARS_PER_LINE);
        char* segment = strtok(buf, DELIMITER);
        token.push_back(segment);

        if(token[0]) { //check if the line is good
            if(atof(token[0]) && beforeMultiply) { // check if number

                double energy = atof(token[0]);
                std::string energyUnit = strtok(0, DELIMITER);
                double dEElec = atof(strtok(0, DELIMITER));
                double dENuclear  = atof(strtok(0, DELIMITER));

                if(debug) {
                    char* buffer[256];
                    sprintf(reinterpret_cast<char *>(buffer), "Energy: %f %s dEElec: %f dENuclear: %f", energy, energyUnit.c_str(), dEElec, dENuclear);
                    std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
                }

                if(energyUnit.compare("eV") == 0) energy_.push_back(energy*1.e-6);
                else if(energyUnit.compare("keV") == 0) energy_.push_back(energy*1.e-3);
                else if(energyUnit.compare("MeV") == 0) energy_.push_back(energy);
                else if(energyUnit.compare("GeV") == 0) energy_.push_back(energy*1.e3);
                dEdx_.push_back(dEElec+dENuclear);
            }

            if(!strcmp(token[0], "Multiply")) {
                beforeMultiply = false;
            }

            if(atof(token[0]) && !beforeMultiply) {
                stoppingConversionPower_.push_back(atof(token[0]));
                stoppingConversionUnit1_.push_back(strtok(0, DELIMITER));
                strtok(0, DELIMITER);
                stoppingConversionUnit2_.push_back(strtok(0, DELIMITER));
            }
        }

        if(!token.empty()) {
            token.clear(); // clear token vector
        }
    }

    for(size_t i = 0; i < stoppingConversionPower_.size(); i++) {
        if(stoppingConversionUnit1_[i].compare("MeV") == 0 && stoppingConversionUnit2_[i].compare("mm") == 0) {
            stoppingConversionPower = stoppingConversionPower_[i];
            if(debug) {
                char* buffer[256];
                sprintf(reinterpret_cast<char *>(buffer), "Conversion: %f %s/%s", stoppingConversionPower, stoppingConversionUnit1_[i].c_str(), stoppingConversionUnit2_[i].c_str());
                std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
            }
        }
    }
    std::transform(dEdx_.begin(), dEdx_.end(), dEdx_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, stoppingConversionPower));
    energySpline.SetPoints(energy_, dEdx_);

    GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
    GL512 = true;

    CalcRemainderErr = 1.e-8;
    AddBackErr = 1.e-8;
    AddBackHighPoint = 250.;

    ReadInitParams();
}

inline void EnergyLoss::ReadFileSimple(const char* srimFile, double stoppingPower) {
    std::ifstream inFile(srimFile);
    ASSERT_WITH_MESSAGE_EL(inFile.is_open(), "Cannot find SRIM file!\n");

    if(debug) std::cout << PrintOutput("DEBUGGING: READING IN SRIM FILE", "red") << std::endl;

    while(!inFile.eof()) {
        char buf[MAX_CHARS_PER_LINE];
        inFile.getline(buf, MAX_CHARS_PER_LINE);
        char* segment = strtok(buf, DELIMITER);
        token.push_back(segment);

        if(token[0]) { //check if the line is good
            if(atof(token[0]) && beforeMultiply) { // check if number
                double energy = atof(token[0]);
                std::string energyUnit = strtok(0, DELIMITER);
                double dEElec = atof(strtok(0, DELIMITER));
                double dENuclear = atof(strtok(0, DELIMITER));

                if(debug) {
                    char* buffer[256];
                    sprintf(reinterpret_cast<char *>(buffer), "Energy: %f %s dEElec: %f dENuclear: %f", energy, energyUnit.c_str(), dEElec, dENuclear);
                    std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
                }

                if(energyUnit.compare("eV") == 0) energy_.push_back(energy*1.e-6);
                else if(energyUnit.compare("keV") == 0) energy_.push_back(energy*1.e-3);
                else if(energyUnit.compare("MeV") == 0) energy_.push_back(energy);
                else if(energyUnit.compare("GeV") == 0) energy_.push_back(energy*1.e3);
                dEdx_.push_back(dEElec+dENuclear);
            }
        }

        if(!token.empty()) {
            token.clear(); // clear token vector
        }
    }

    std::transform(dEdx_.begin(), dEdx_.end(), dEdx_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, stoppingPower));
    energySpline.SetPoints(energy_, dEdx_);

    GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
    GL512 = true;

    CalcRemainderErr = 1.e-8;
    AddBackErr = 1.e-8;
    AddBackHighPoint = 250.;

    ReadInitParams();
}

inline void EnergyLoss::ReadBasicdEdx(const char* inputFile) {
    std::ifstream inFile(inputFile);
    ASSERT_WITH_MESSAGE_EL(inFile.is_open(), "Cannot find input file!\n");

    double energy, dEdx;
    while (inFile >> energy >> dEdx) {
        energy_.push_back(energy);
        dEdx_.push_back(dEdx);
    }

    energySpline.SetPoints(energy_, dEdx_);

    GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
    GL512 = true;

    CalcRemainderErr = 1.e-8;
    AddBackErr = 1.e-8;
    AddBackHighPoint = 250.;

    ReadInitParams();
}

inline void EnergyLoss::ReadLISEFile(const char* inputFile) {
    std::ifstream inFile(inputFile);
    ASSERT_WITH_MESSAGE_EL(inFile.is_open(), "Cannot find input file!\n");

    std::string str;
    std::getline(inFile, str); // skip the first line

    double energy, dEdx0, dEdx1, dEdx2, dEdx3, dEdx4, dEdx5, dummy;
    while (inFile >> energy >> dEdx0 >> dummy >> dEdx1 >> dummy >> dEdx2 >> dummy >> dEdx3 >> dummy >> dEdx4 >> dummy >> dEdx5) {
        lise_energy_mevu_.push_back(energy);
        lise_dedx_columns_[0].push_back(dEdx0);
        lise_dedx_columns_[1].push_back(dEdx1);
        lise_dedx_columns_[2].push_back(dEdx2);
        lise_dedx_columns_[3].push_back(dEdx3);
        lise_dedx_columns_[4].push_back(dEdx4);
        lise_dedx_columns_[5].push_back(dEdx5);
    }

    inFile.close();
}

inline void EnergyLoss::ReadLISEdEdx(const char* inputFile, double mass) {
    ReadLISEdEdx(inputFile, mass, 10.0, 1);
}

inline void EnergyLoss::ReadLISEdEdx(const char* inputFile, double mass, int column) {
    ReadLISEdEdx(inputFile, mass, 10.0, column);
}

inline void EnergyLoss::ReadLISEdEdx(const char* inputFile, double mass, double density) {
    ReadLISEdEdx(inputFile, mass, density, 1);
}

inline void EnergyLoss::ReadLISEdEdx(const char* inputFile, double mass, double density, int column) {
    ReadLISEFile(inputFile);

    energy_ = lise_energy_mevu_;
    std::transform(energy_.begin(), energy_.end(), energy_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mass));

    dEdx_ = lise_dedx_columns_[column];
    std::transform(dEdx_.begin(), dEdx_.end(), dEdx_.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, density*100.));

    energySpline.SetPoints(energy_, dEdx_);

    GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
    GL512 = true;

    CalcRemainderErr = 1.e-8;
    AddBackErr = 1.e-8;
    AddBackHighPoint = 250.;

    ReadInitParams();
}

inline void EnergyLoss::SetDebug(bool flag) {
    debug = flag;
}

inline double EnergyLoss::CalcRemainder(double initialEnergy, double distance) {
    if(distance == 0.) return initialEnergy;

    distance = fabs(distance);

    if(initialEnergy < CalcRemainderErr) return 0.;

    if(debug) std::cout << PrintOutput("DEBUGGING: CalcRemainder", "red") << std::endl;

    distance = fabs(distance);

    double maxRange = CalcRange(initialEnergy, 0.);
    if(debug) std::cout << PrintOutput("Max Range of particle in the material (mm): ", "green") << maxRange << std::endl;

    if(distance > maxRange) return 0.;

    double lowEnergy = 0.;
    double highEnergy = initialEnergy;
    double guessEnergy = (highEnergy + lowEnergy)/2.0;

    double range = CalcRange(initialEnergy, guessEnergy);
    while(fabs(range - distance) > CalcRemainderErr) {
        if(debug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "Low Energy: %f; Guess Energy: %f; High Energy: %f Range: %f; Distance: %f",
                    lowEnergy, guessEnergy, highEnergy, range, distance);
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        if(range > distance) {
            lowEnergy = guessEnergy;
            guessEnergy = (lowEnergy+highEnergy)/2.0;
        } else {
            highEnergy = guessEnergy;
            guessEnergy = (lowEnergy+highEnergy)/2.0;
        }
        range = CalcRange(initialEnergy, guessEnergy);
    }
    return guessEnergy;
}

inline double EnergyLoss::AddBack(double finalEnergy, double distance) {
    if(distance == 0.) return finalEnergy;

    distance = fabs(distance);

    if(debug) std::cout << PrintOutput("DEBUGGING: AddBack", "red") << std::endl;

    double lowEnergy = finalEnergy;
    double highEnergy = AddBackHighPoint;
    double guessEnergy = (highEnergy + lowEnergy)/2.0;

    double range = CalcRange(guessEnergy, finalEnergy);
    while(fabs(range - distance) > AddBackErr) {
        if(debug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "Low Energy: %f; Guess Energy: %f; High Energy: %f Range: %f; Distance: %f",
                    lowEnergy, guessEnergy, highEnergy, range, distance);
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        if(range < distance) {
            lowEnergy = guessEnergy;
            guessEnergy = (lowEnergy + highEnergy)/2.0;
        } else {
            highEnergy = guessEnergy;
            guessEnergy = (lowEnergy + highEnergy)/2.0;
        }
        range = CalcRange(guessEnergy, finalEnergy);

        if(guessEnergy > (AddBackHighPoint-0.1)) {
            printf("Error: EnergyLoss::AddBack above starting high guess!\n");
            printf("Use the AddBackHigh to set the starting high guess higher!\n");
            return guessEnergy;
        }
    }
    return guessEnergy;
}

inline double EnergyLoss::CalcRange(double initialEnergy, double remainder) {
    if(initialEnergy < 0.) return 0.;
    if(remainder < 0.) return 0.;

    if(remainder > initialEnergy) return 0;

    if(initialEnergy == remainder) return 0;
    double distance;

    if(GL16) distance = CalcRangeGL16(remainder, initialEnergy);
    else if(GL32) distance = CalcRangeGL32(remainder, initialEnergy);
    else if(GL64) distance = CalcRangeGL64(remainder, initialEnergy);
    else if(GL128) distance = CalcRangeGL128(remainder, initialEnergy);
    else if(GL256) distance = CalcRangeGL256(remainder, initialEnergy);
    else if(GL512) distance = CalcRangeGL512(remainder, initialEnergy);
    else if(GL1024) distance = CalcRangeGL1024(remainder, initialEnergy);
    else distance = CalcRangeGL512(remainder, initialEnergy);

    return distance;
}

inline void EnergyLoss::CalcRemainderError(double err) {
    CalcRemainderErr = err;
}

inline void EnergyLoss::AddBackError(double err) {
    AddBackErr = err;
}

inline void EnergyLoss::AddBackHigh(double high) {
    AddBackHighPoint = high;
}

inline void EnergyLoss::UseGL16() {
    GL32 = GL64 = GL128 = GL256 = GL512 = GL1024 = false;
    GL16 = true;
}

inline void EnergyLoss::UseGL32() {
    GL16 = GL64 = GL128 = GL256 = GL512 = GL1024 = false;
    GL32 = true;
}

inline void EnergyLoss::UseGL64() {
    GL16 = GL32 = GL128 = GL256 = GL512 = GL1024 = false;
    GL64 = true;
}

inline void EnergyLoss::UseGL128() {
    GL16 = GL32 = GL64 = GL256 = GL512 = GL1024 = false;
    GL128 = true;
}

inline void EnergyLoss::UseGL256() {
    GL16 = GL32 = GL64 = GL128 = GL512 = GL1024 = false;
    GL256 = true;
}

inline void EnergyLoss::UseGL512() {
    GL16 = GL32 = GL64 = GL128 = GL256 = GL1024 = false;
    GL512 = true;
}

inline void EnergyLoss::UseGL1024() {
    GL16 = GL32 = GL64 = GL128 = GL256 = GL512 = false;
    GL1024 = true;
}

inline double EnergyLoss::CalcRangeGL16(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    if(calcRangeDebug) std::cout << PrintOutput("DEBUGGING: CalcRange", "red") << std::endl;
    for(int i = 7; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w16[i]: %f; -x16[i]: %f energySpline(t0 + dt*(-x16[i])): %f",
                    i, w16[i], -x16[i], energySpline(t0 + dt*(-x16[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w16[i]*(1./energySpline(t0 + dt*(-x16[i])));
    }
    for(int i = 0; i < 8; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w16[i]: %f; x16[i]: %f energySpline(t0 + dt*(x16[i])): %f",
                    i, w16[i], x16[i], energySpline(t0 + dt*(x16[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w16[i]*(1./energySpline(t0 + dt*x16[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL32(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 15; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w32[i]: %f; -x32[i]: %f energySpline(t0 + dt*(-x32[i])): %f",
                    i, w32[i], -x32[i], energySpline(t0 + dt*(-x32[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w32[i]*(1./energySpline(t0 + dt*(-x32[i])));
    }
    for(int i = 0; i < 16; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w32[i]: %f; x32[i]: %f energySpline(t0 + dt*(x32[i])): %f",
                    i, w32[i], x32[i], energySpline(t0 + dt*(x32[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w32[i]*(1./energySpline(t0 + dt*x32[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL64(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 31; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w64[i]: %f; -x64[i]: %f energySpline(t0 + dt*(-x64[i])): %f",
                    i, w64[i], -x64[i], energySpline(t0 + dt*(-x64[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w64[i]*(1./energySpline(t0 + dt*(-x64[i])));
    }
    for(int i = 0; i < 32; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w64[i]: %f; x64[i]: %f energySpline(t0 + dt*(x64[i])): %f",
                    i, w64[i], x64[i], energySpline(t0 + dt*(x64[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w64[i]*(1./energySpline(t0 + dt*x64[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL128(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 63; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w128[i]: %f; -x128[i]: %f energySpline(t0 + dt*(-x128[i])): %f",
                    i, w128[i], -x128[i], energySpline(t0 + dt*(-x128[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w128[i]*(1./energySpline(t0 + dt*(-x128[i])));
    }
    for(int i = 0; i < 64; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w128[i]: %f; x128[i]: %f energySpline(t0 + dt*(x128[i])): %f",
                    i, w128[i], x128[i], energySpline(t0 + dt*(x128[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w128[i]*(1./energySpline(t0 + dt*x128[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL256(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 127; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w256[i]: %f; -x256[i]: %f energySpline(t0 + dt*(-x256[i])): %f",
                    i, w256[i], -x256[i], energySpline(t0 + dt*(-x256[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w256[i]*(1./energySpline(t0 + dt*(-x256[i])));
    }
    for (int i = 0; i < 128; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w256[i]: %f; x256[i]: %f energySpline(t0 + dt*(x256[i])): %f",
                    i, w256[i], x256[i], energySpline(t0 + dt*(x256[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w256[i]*(1./energySpline(t0 + dt*x256[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL512(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 255; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w512[i]: %f; -x512[i]: %f energySpline(t0 + dt*(-x512[i])): %f",
                    i, w512[i], -x512[i], energySpline(t0 + dt*(-x512[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w512[i]*(1./energySpline(t0 + dt*(-x512[i])));
    }
    for(int i = 0; i < 256; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w512[i]: %f; x512[i]: %f energySpline(t0 + dt*(x512[i])): %f",
                    i, w512[i], x512[i], energySpline(t0 + dt*(x512[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w512[i]*(1./energySpline(t0 + dt*x512[i]));
    }
    return result*dt;
}

inline double EnergyLoss::CalcRangeGL1024(double a, double b) {
    double t0 = (a + b)/2.0;
    double dt = (b - a)/2.0;
    double result = 0.;
    for(int i = 511; i > -1; i--) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w1024[i]: %f; -x1024[i]: %f energySpline(t0 + dt*(-x1024[i])): %f",
                    i, w1024[i], -x1024[i], energySpline(t0 + dt*(-x1024[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w1024[i]*(1./energySpline(t0 + dt*(-x1024[i])));
    }
    for(int i = 0; i < 512; i++) {
        if(calcRangeDebug) {
            char* buffer[256];
            sprintf(reinterpret_cast<char *>(buffer), "i: %d; w1024[i]: %f; x1024[i]: %f energySpline(t0 + dt*(x1024[i])): %f",
                    i, w1024[i], x1024[i], energySpline(t0 + dt*(x1024[i])));
            std::cout << PrintOutput(reinterpret_cast<const char *>(buffer), "green") << std::endl;
        }
        result += w1024[i]*(1./energySpline(t0 + dt*x1024[i]));
    }
    return result * dt;
}

inline void EnergyLoss::ReadInitParams() {
    std::ifstream inFile("EnergyLoss.dat", std::ios::in | std::ifstream::binary);
    ASSERT_WITH_MESSAGE_EL(inFile.is_open(), "Cannot find EnergyLoss.dat file!\n");

    inFile.read(reinterpret_cast<char*> (x16), sizeof(x16));
    inFile.read(reinterpret_cast<char*> (w16), sizeof(w16));
    inFile.read(reinterpret_cast<char*> (x32), sizeof(x32));
    inFile.read(reinterpret_cast<char*> (w32), sizeof(w32));
    inFile.read(reinterpret_cast<char*> (x64), sizeof(x64));
    inFile.read(reinterpret_cast<char*> (w64), sizeof(w64));
    inFile.read(reinterpret_cast<char*> (x128), sizeof(x128));
    inFile.read(reinterpret_cast<char*> (w128), sizeof(w128));
    inFile.read(reinterpret_cast<char*> (x256), sizeof(x256));
    inFile.read(reinterpret_cast<char*> (w256), sizeof(w256));
    inFile.read(reinterpret_cast<char*> (x512), sizeof(x512));
    inFile.read(reinterpret_cast<char*> (w512), sizeof(w512));
    inFile.read(reinterpret_cast<char*> (x1024), sizeof(x1024));
    inFile.read(reinterpret_cast<char*> (w1024), sizeof(w1024));
    inFile.close();
}

inline std::string EnergyLoss::PrintOutput(std::string Output, std::string Color) {
    int ColorCode = 0;
    if(Color.compare("red") == 0){
        ColorCode = 31;
    } else if(Color.compare("green") == 0){
        ColorCode = 32;
    } else if(Color.compare("yellow") == 0){
        ColorCode = 33;
    } else if(Color.compare("blue") == 0) {
        ColorCode = 34;
    } else if(Color.compare("magenta") == 0) {
        ColorCode = 35;
    } else if(Color.compare("cyan") == 0) {
        ColorCode = 36;
    } else {
        return Output;
    }
    char buffer[256];
    sprintf(buffer, "\033[1;%dm%s\033[0m", ColorCode, Output.c_str());
    return buffer;
}

inline void EnergyLoss::UseCalcRangeDebug() {
    calcRangeDebug = true;
}

inline double GetdEdx(const EnergyLoss& elClass, double energy) {
    return elClass.energySpline(energy);
}

#endif
