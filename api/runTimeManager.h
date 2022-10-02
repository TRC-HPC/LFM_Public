#ifndef IRUNTIMEMANAGER_H
#define IRUNTIMEMANAGER_H

#include <ostream>
#include <vector>

class IRunTimeManager
{
public:
    // Destructor
    virtual ~IRunTimeManager() {}

    // --------- //
    // Time Info //
    // --------- //        
    // Start / End / Curr Time
    virtual double getStartTime() const = 0;
    virtual double getEndTime() const = 0;
    virtual double getCurrentTime() const = 0;

    // Delta Time    
    virtual double getDeltaTime() const = 0;
    virtual void setDeltaTime(const double dDeltaTime) = 0;
    virtual bool isAdjustDeltaTime() const = 0;

    // Time Step
    virtual int getTimeStep() const = 0;

    // Advance Time
    virtual void advanceTime() = 0;

    // ----------- //
    // Run Control //
    // ----------- //
    // Check if running
    virtual bool isRunning() = 0;

    // Write Results
    virtual bool isWriteTime() = 0;
    virtual void writeResults() = 0;

    // Get Data Pointer
    virtual void* getRunPointer() = 0;
};

#endif // IRUNTIMEMANAGER_H
