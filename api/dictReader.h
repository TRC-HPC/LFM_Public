#ifndef DICTREADER_H
#define DICTREADER_H

#include <ostream>
#include <vector>

class IDictReader
{
public:
    // Destructor
    virtual ~IDictReader() {}

    // Read Word
    virtual std::string readDictWord(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false) = 0;
    virtual std::string readDictWord(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) = 0;

    // Read Filename
    virtual std::string readDictFilename(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false) = 0;
    virtual std::string readDictFilename(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) = 0;

    // Read Bool
    virtual bool readDictBool(std::string paramName, std::vector<std::string> subDictList, const bool bDefault = true, const bool bAllowDefault = false) = 0;
    virtual bool readDictBool(std::string paramName, std::string subDict = "", const bool bDefault = true, const bool bAllowDefault = false) = 0;

    // Read scalar
    virtual double readDictScalar(std::string paramName, std::vector<std::string> subDictList, const double dDefault = 0, const bool bAllowDefault = false) = 0;
    virtual double readDictScalar(std::string paramName, std::string subDict = "", const double dDefault = 0, const bool bAllowDefault = false) = 0;

    // Read String
    virtual std::string readDictString(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false) = 0;
    virtual std::string readDictString(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) = 0;

    // Read vector
    virtual std::vector<double> readDictVector(std::string paramName, std::vector<std::string> subDictList) = 0;
    virtual std::vector<double> readDictVector(std::string paramName, std::string subDict = "") = 0;

    // Get Keys
    virtual std::vector<std::string> getKeyList(std::vector<std::string> subDictList) = 0;
    virtual std::vector<std::string> getKeyList(std::string subDict = "") = 0;

    // Get Dicts
    virtual std::vector<std::string> getDictList(std::vector<std::string> subDictList) = 0;
    virtual std::vector<std::string> getDictList(std::string subDict = "") = 0;
};

#endif
