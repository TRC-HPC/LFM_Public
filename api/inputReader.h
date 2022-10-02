#ifndef INPUT_READER_H
#define INPUT_READER_H

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <getopt.h>

class CInputReader
{
public:
    // Constructor
    CInputReader() {}

    // Initialize
    bool initialize(int argc, char *argv[]);

protected:
    // Show Usage
    static void showUsage(char *argv[]);

    // Parse Halo Communication Type
    bool parseHaloCommType(const char* sCommType);

    // Read Runtime Arguments
    bool readRuntimeArgs(int argc, char *argv[]);

    // Extract actual time
    double extractTime(std::string sTime) const;

public:
    // Parallel
    bool m_bParallel;

    // ---------------- //
    // controlDict Info //
    // ---------------- //
    int m_nCommType;
    int m_nHaloCommType;
    bool m_bIsDoublePrecision;
    bool m_bHaveProbes;
    bool m_bHaveSampling;
    bool m_bHaveAverage;
    bool m_bHaveForces;
    bool m_bHaveResidual;
    int m_nSaveForcesSteps;
    int m_nPrintInfoFrequency;
    bool m_bSaveResiduals;
    bool m_bSaveBlendFactor;
    bool m_bSaveRank;
    bool m_bSaveCellOrder;
    double m_dCFLMax;
    double m_tStartAverage;

    //---------------//
    // fvScheme Info //
    //---------------//
    int m_nSolver;
    int m_nDimension;
    int m_nRkStepCount;
    bool m_bMinmodInterpolation;
    bool m_bConstantTimeStep;

    //-----------------//
    // spongeDict Info //
    //-----------------//
    double m_dPInf;
    double m_dTInf;    
    double m_dUInf[3];
    double m_dLs;
    double m_dMach;
    double m_dK;

    // ----------------------------- //
    // thermophysicalProperties Info //
    // ----------------------------- //
    double m_dCp;
    double m_dMolWeight;
    double m_dMu0;
    double m_dPrandtl;

    // ----------------------------- //
    //    turbulenceProperties Info  //
    // ----------------------------- //
    bool m_bisLaminar;
};


#endif // INPUT_READER_H