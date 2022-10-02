#ifndef CRUNTIMEMANAGEROF_H
#define CRUNTIMEMANAGEROF_H

#include <ostream>
#include <vector>
#include <map>
#include "runTimeManager.h"

class CRunTimeManagerOF : public IRunTimeManager
{
public:
    // Constructor
    CRunTimeManagerOF(const int nRank);
    virtual ~CRunTimeManagerOF();

    // --------- //
    // Time Info //
    // --------- //        
    // Start / End / Curr Time
    virtual double getStartTime() const;
    virtual double getEndTime() const;
    virtual double getCurrentTime() const;

    // Delta Time    
    virtual double getDeltaTime() const;
    virtual void setDeltaTime(const double dDeltaTime);
    virtual bool isAdjustDeltaTime() const;

    // Time Step
    virtual int getTimeStep() const;

    // Advance Time
    virtual void advanceTime();

    // ----------- //
    // Run Control //
    // ----------- //
    // Check if running
    virtual bool isRunning();

    // Write Results
    virtual bool isWriteTime();
    virtual void writeResults();

    // Get Data Pointer
    virtual void* getRunPointer() {return m_pRunTime;}

protected:
    // RunTime
    void* m_pRunTime;
};

#endif //CRUNTIMEMANAGEROF_H
