
#include "dictReaderOF.h"
#include "fvCFD.H"


template std::string    CDictReaderOF::getData<std::string, word>(std::string , std::vector<std::string>);
template std::string    CDictReaderOF::getData<std::string, fileName>(std::string , std::vector<std::string>);
template bool           CDictReaderOF::getData<bool, bool>(std::string , std::vector<std::string>);
template double         CDictReaderOF::getData<double, scalar>(std::string , std::vector<std::string>);

template std::string  CDictReaderOF::getDataOrDefault<std::string, word>(std::string , std::vector<std::string>, std::string);
template std::string  CDictReaderOF::getDataOrDefault<std::string, fileName>(std::string , std::vector<std::string>, std::string);
template bool         CDictReaderOF::getDataOrDefault<bool, bool>(std::string , std::vector<std::string>, bool);
template double       CDictReaderOF::getDataOrDefault<double, scalar>(std::string , std::vector<std::string>, double);



// -------------------------------------------------------------------------- //
CDictReaderOF::CDictReaderOF(std::string dictName, std::string dictDirectory)
{   
    Foam::Time runTime(".", ".");

    m_pIODict = new IOdictionary(
        IOobject
        (
            dictName,    // dictionary name
            dictDirectory,     // dict is found in "constant"
            runTime,                   // registry for the dict
            IOobject::MUST_READ,    // must exist, otherwise failure
            IOobject::NO_WRITE,      // dict is only read by the solver
            false,
            false
        )
    );
}

// -------------------------------------------------------------------------- //
CDictReaderOF::~CDictReaderOF()
{
    IOdictionary* pIODict = static_cast<IOdictionary *>(m_pIODict);
    if (pIODict)
        delete pIODict;
    m_pIODict = nullptr;
}

// -------------------------------------------------------------------------- //
std::string CDictReaderOF::readDictWord(std::string paramName, std::vector<std::string> subDictList, std::string sDefault, const bool bAllowDefault)
{
    return bAllowDefault ? getDataOrDefault<std::string, word>(paramName, subDictList, sDefault) : getData<std::string, word>(paramName, subDictList);
}

// -------------------------------------------------------------------------- //
std::string CDictReaderOF::readDictFilename(std::string paramName, std::vector<std::string> subDictList, std::string sDefault, const bool bAllowDefault)
{
    return bAllowDefault ? getDataOrDefault<std::string, fileName>(paramName, subDictList, sDefault) : getData<std::string, fileName>(paramName, subDictList);
}

// -------------------------------------------------------------------------- //
bool CDictReaderOF::readDictBool(std::string paramName, std::vector<std::string> subDictList, bool bDefault, const bool bAllowDefault)
{
    return bAllowDefault ? getDataOrDefault<bool, bool>(paramName, subDictList, bDefault) : getData<bool, bool>(paramName, subDictList);
}

// -------------------------------------------------------------------------- //
double CDictReaderOF::readDictScalar(std::string paramName, std::vector<std::string> subDictList, const double dDefault, const bool bAllowDefault)
{
    return bAllowDefault ? getDataOrDefault<double, scalar>(paramName, subDictList, dDefault) : getData<double, scalar>(paramName, subDictList);
}

// -------------------------------------------------------------------------- //
std::string CDictReaderOF::readDictString(std::string paramName, std::vector<std::string> subDictList, const std::string sDefault, const bool bAllowDefault)
{
    const dictionary* pDict = static_cast<const dictionary *>(getDict(subDictList));
    return pDict->lookup(paramName).toString();
}

// -------------------------------------------------------------------------- //
std::vector<double> CDictReaderOF::readDictVector(std::string paramName, std::vector<std::string> subDictList)
{
    const dictionary* pDict = static_cast<const dictionary *>(getDict(subDictList));
    vector v = vector(pDict->lookup(paramName));
    std::vector<double> ret;
    ret.push_back(v.x());
    ret.push_back(v.y());
    ret.push_back(v.z());
    return ret;
}

// -------------------------------------------------------------------------- //
std::vector<std::string> CDictReaderOF::extractSubDictList(std::string subDict)
{
    std::vector<std::string> subDictList;
    if (!subDict.empty())
    {
        std::string seperator = "/";
        std::string::size_type prev_pos = 0, pos = 0;
        while((pos = subDict.find(seperator, pos)) != std::string::npos)
        {
            std::string substring(subDict.substr(prev_pos, pos-prev_pos) );
            subDictList.push_back(substring);
            prev_pos = ++pos;
        }

        subDictList.push_back(subDict.substr(prev_pos, pos-prev_pos)); // Last word
    }
    return subDictList;
}

// -------------------------------------------------------------------------- //
const void* CDictReaderOF::getDict(std::vector<std::string> subDictList)
{
    IOdictionary* pIODict = static_cast<IOdictionary *>(m_pIODict);
    const dictionary* pDict = &pIODict->topDict();
    for (size_t i=0; i < subDictList.size(); i++)
        pDict = &pDict->optionalSubDict(subDictList[i]);
    return pDict;
}

// -------------------------------------------------------------------------- //
std::vector<std::string> CDictReaderOF::getKeyOrDictList(const bool bKey, std::vector<std::string> subDictList)
{
    std::vector<std::string> retList;

    const dictionary* pDict = static_cast<const dictionary *>(getDict(subDictList));

    // toc list
    wordList tocList = pDict->toc();
    const int nTocCount = tocList.size();
    for (int nTocIndex=0; nTocIndex < nTocCount; nTocIndex++)
    {
        const bool bIsDict = pDict->isDict(tocList[nTocIndex]);
        if (bIsDict != bKey) 
            retList.push_back(tocList[nTocIndex]);
    }
    
    return retList;
}

// -------------------------------------------------------------------------- //
template<class retT, class reqT> retT CDictReaderOF::getData(std::string paramName, std::vector<std::string> subDictList)
{
    const dictionary* pDict = static_cast<const dictionary *>(getDict(subDictList));
    retT ret = pDict->get<reqT>(paramName);
    return ret;
}

// -------------------------------------------------------------------------- //
template<class retT, class reqT> retT CDictReaderOF::getDataOrDefault(std::string paramName, std::vector<std::string> subDictList, retT defaultValue)
{
    const dictionary* pDict = static_cast<const dictionary *>(getDict(subDictList));
    retT ret = pDict->getOrDefault<reqT>(paramName, defaultValue);
    return ret;
}


// ************************************************************************* //
