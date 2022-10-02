#include "runTimeManagerOF.h"
#include "fvCFD.H"

// -------------------------------------------------------------------------- //
CRunTimeManagerOF::CRunTimeManagerOF(const int nRank)
{
    // Case Name
    std::string sCaseName = ".";
    if (nRank >= 0) sCaseName = "processor" + std::to_string(nRank);

    // RunTime
    m_pRunTime = new Foam::Time(Foam::Time::controlDictName, ".", sCaseName);
}

// -------------------------------------------------------------------------- //
CRunTimeManagerOF::~CRunTimeManagerOF()
{    
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    if (m_pRunTime)
        delete pRunTime;
}

// -------------------------------------------------------------------------- //
double CRunTimeManagerOF::getStartTime() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.startTime().value();
}

// -------------------------------------------------------------------------- //
double CRunTimeManagerOF::getEndTime() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;    
    return runTime.endTime().value();
}

// -------------------------------------------------------------------------- //
double CRunTimeManagerOF::getCurrentTime() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.value();
}

// -------------------------------------------------------------------------- //
double CRunTimeManagerOF::getDeltaTime() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.deltaT().value();
}

// -------------------------------------------------------------------------- //
void CRunTimeManagerOF::setDeltaTime(const double dDeltaTime)
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    runTime.setDeltaT(dDeltaTime);
}

// -------------------------------------------------------------------------- //
bool CRunTimeManagerOF::isAdjustDeltaTime() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.isAdjustTimeStep();
}

// -------------------------------------------------------------------------- //
int CRunTimeManagerOF::getTimeStep() const
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.timeIndex();
}

// -------------------------------------------------------------------------- //
void CRunTimeManagerOF::advanceTime()
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    runTime++;
}

// -------------------------------------------------------------------------- //
bool CRunTimeManagerOF::isRunning()
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.run();
}

// -------------------------------------------------------------------------- //
bool CRunTimeManagerOF::isWriteTime()
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    return runTime.writeTime();
}

// -------------------------------------------------------------------------- //
void CRunTimeManagerOF::writeResults()
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;
    runTime.write();
}
