#ifndef CDICTREADEROF_H
#define CDICTREADEROF_H

#include <ostream>
#include <vector>

#include "dictReader.h"

class CDictReaderOF : public IDictReader
{
public:
    // Constructor / Destructor
    CDictReaderOF(std::string dictName, std::string dictDirectory);
    virtual ~CDictReaderOF();

    // Read Word    
    virtual std::string readDictWord(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false);
    virtual std::string readDictWord(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) {return readDictWord(paramName, extractSubDictList(subDict), sDefault, bAllowDefault);}

    // Read Filename
    virtual std::string readDictFilename(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false);
    virtual std::string readDictFilename(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) {return readDictFilename(paramName, extractSubDictList(subDict), sDefault, bAllowDefault);}

    // Read Bool
    virtual bool readDictBool(std::string paramName, std::vector<std::string> subDictList, const bool bDefault = true, const bool bAllowDefault = false);
    virtual bool readDictBool(std::string paramName, std::string subDict = "", const bool bDefault = true, const bool bAllowDefault = false) {return readDictBool(paramName, extractSubDictList(subDict), bDefault, bAllowDefault);}

    // Read scalar
    virtual double readDictScalar(std::string paramName, std::vector<std::string> subDictList, const double dDefault = 0, const bool bAllowDefault = false);
    virtual double readDictScalar(std::string paramName, std::string subDict = "", const double dDefault = 0, const bool bAllowDefault = false) {return readDictScalar(paramName, extractSubDictList(subDict), dDefault, bAllowDefault);}

    // Read String
    virtual std::string readDictString(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault = "", const bool bAllowDefault = false);
    virtual std::string readDictString(std::string paramName, std::string subDict = "", const std::string sDefault = "", const bool bAllowDefault = false) {return readDictString(paramName, extractSubDictList(subDict), sDefault, bAllowDefault);}

    // Read vector
    virtual std::vector<double> readDictVector(std::string paramName, std::vector<std::string> subDictList);
    virtual std::vector<double> readDictVector(std::string paramName, std::string subDict = "") {return readDictVector(paramName, extractSubDictList(subDict));}

    // Get Keys
    virtual std::vector<std::string> getKeyList(std::vector<std::string> subDictList) {return getKeyOrDictList(true, subDictList);}
    virtual std::vector<std::string> getKeyList(std::string subDict = "") {return getKeyOrDictList(true, subDict);}

    // Get Dicts
    virtual std::vector<std::string> getDictList(std::vector<std::string> subDictList) {return getKeyOrDictList(false, subDictList);}
    virtual std::vector<std::string> getDictList(std::string subDict = "") {return getKeyOrDictList(false, subDict);}

protected:
    // Extract subDict Hirarchy
    std::vector<std::string> extractSubDictList(std::string subDict);

    // Get Dict
    const void* getDict(std::vector<std::string> subDictList);
    const void* getDict(std::string subDict) {return getDict(extractSubDictList(subDict));}
    
    // Get KeyOrSubDict
    std::vector<std::string> getKeyOrDictList(const bool bKey, std::vector<std::string> subDictList);
    std::vector<std::string> getKeyOrDictList(const bool bKey, std::string subDict = "") {return getKeyOrDictList(bKey, extractSubDictList(subDict));}

    // Get Class (not default)
    template<typename retT, typename reqT> retT getData(std::string paramName, std::vector<std::string> subDictList);
    template<typename retT, typename reqT> retT getData(std::string paramName, std::string subDict) {return getData<retT,reqT>(paramName, extractSubDictList(subDict));}

    // Get Class (with default)
    template<typename retT, typename reqT> retT getDataOrDefault(std::string paramName, std::vector<std::string> subDictList, const retT defaultValue);
    template<typename retT, typename reqT> retT getDataOrDefault(std::string paramName, std::string subDict, const retT defaultValue) {return getDataOrDefault<retT,reqT>(paramName, extractSubDictList(subDict), defaultValue);}
protected:
    // Dictionary 
    void* m_pIODict;
};

#endif
