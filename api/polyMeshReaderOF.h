#ifndef CPOLYMESHREADER_H
#define CPOLYMESHREADER_H

#include <ostream>
#include <vector>
#include <map>
#include "meshReader.h"

class CPolyMeshReaderOF : public IMeshReader
{
public:
    // Constructor
    CPolyMeshReaderOF(void* pRunTime);
    virtual ~CPolyMeshReaderOF();

    // Points
    virtual int getPointCount() const;
    virtual void getPoint(const int nPointIndex, double *pXYZ) const;
    virtual int getPointProcAddressing(const int nPointIndex) const;

    // Faces
    virtual int getFaceCount() const;
    virtual void getFaceAreaNormal(const int nFaceIndex, double* pN) const;
    virtual void getFaceCenter(const int nFaceIndex, double* pXYZ) const;
    virtual int getFaceOwner(const int nFaceIndex) const;
    virtual int getFaceNeighbour(const int nFaceIndex) const;
    virtual int getFaceBoundary(const int nFaceIndex) const;    
    virtual int getFaceIndexInsideBoundary(const int nFaceIndex) const;
    virtual int getFacePointCount(const int nFaceIndex) const;
    virtual int getFaceProcAddressing(const int nFaceIndex) const;
    virtual int getFacePointIndex(const int nFaceIndex, const int nIndex) const;
    virtual std::vector<int> getFacePointIndexList(const int nFaceIndex) const;
    virtual bool getFaceIsValid(const int nFaceIndex) const {return m_bValidFaceFlagList[nFaceIndex];}
    virtual int getCyclicFaceIndex(const int nFaceIndex) const;
    virtual int getFaceId(const int nFaceIndex) const;
    
    // Cells
    virtual int getCellCount() const;
    virtual double getCellVolume(const int nCellIndex) const;
    virtual void getCellCenter(const int nCellIndex, double* pXYZ) const;
    virtual int getCellFaceCount(const int nCellIndex) const;
    virtual int getCellFaceIndex(const int nCellIndex, const int nIndex) const;
    virtual bool getCellFaceOwner(const int nCellIndex, const int nIndex) const;
    virtual int getCellProcAddressing(const int nCellIndex) const;
    virtual std::vector<int> getCellPointList(const int nCellIndex) const;
    virtual int getCellValidFaceCount(const int nCellIndex) const;
    virtual int getCellValidFaceIndex(const int nCellIndex, const int nIndex) const;
    virtual double getCellScalar(const int nCellIndex, const int nFieldIndex) const;
    virtual void getCellVector(const int nCellIndex, const int nFieldIndex, double* pVec) const;    

    // Boundary
    virtual int getBoundaryCount() const;
    virtual std::string getBoundaryName(const int nBoundaryIndex) const;
    virtual std::string getBoundaryType(const int nBoundaryIndex) const;
    virtual int getBoundaryFaceStart(const int nBoundaryIndex) const;
    virtual int getBoundaryFaceEnd(const int nBoundaryIndex) const;
    virtual std::vector<int> getBoundaryFaceCellList(const int nBoundaryIndex) const;
    virtual int getBoundaryProcAddressing(const int nBoundaryIndex) const;
    virtual int getBoundaryProcessorRank(const int nBoundaryIndex) const;
    virtual int getBoundaryTag(const int nBoundaryIndex) const;
    virtual int getBoundaryCyclicPairIndex(const int nBoundaryIndex) const;
    virtual std::vector<int> getNeighbourCellList(const int nBoundaryIndex) const;

    // Submesh
    virtual int getMask(const int nIndex) const {return 1 << nIndex;}
    virtual void initializeSubmesh();
    virtual int getSubmeshCount() const {return m_nSubmeshCount;}
    virtual int getCellSubmeshIndex(const int nCellIndex) const {return m_nCellSubmeshList[nCellIndex];}
    virtual std::vector<int> getSubmeshPointIndexList(const int nSubmeshIndex) const {return getMaskedList(m_nPointMaskList, getMask(nSubmeshIndex));}
    virtual std::vector<int> getSubmeshFaceIndexList(const int nSubmeshIndex) const {return getMaskedList(m_nFaceMaskList, getMask(nSubmeshIndex));}
    virtual std::vector<int> getSubmeshCellIndexList(const int nSubmeshIndex) const {return getMaskedList(m_nCellMaskList, getMask(nSubmeshIndex));}
    virtual int getSubmeshPointIndex(const int nSubmeshIndex, const int nPointIndex)  {return m_nSubmeshPointMappingList[nSubmeshIndex][nPointIndex];}
    virtual int getSubmeshFaceIndex(const int nSubmeshIndex, const int nFaceIndex) {return m_nSubmeshFaceMappingList[nSubmeshIndex][nFaceIndex];}
    virtual int getSubmeshCellIndex(const int nSubmeshIndex, const int nCellIndex) {return m_nSubmeshCellMappingList[nSubmeshIndex][nCellIndex];}
    virtual std::vector<double> getCellDistanceFromBoundary(const int nSubmeshIndex, const std::vector<int> nBoundaryIndexList) const;

    // Locality Test
    virtual double getLocalityScore(const int nMethod = 0) const;

    // Scalar & Vector Fields
    virtual int readScalarField(const std::string sFieldName);
    virtual int readVectorField(const std::string sFieldName);
    virtual int createScalarField(const std::string sFieldName);

    // Update Fields
    virtual void updateScalarField(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex = -1, const double dFactor = 1.0);
    virtual void updateScalarField(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex = -1, const float fFactor = 1.0);
    virtual void updateScalarFieldSqrt(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex = -1, const double dFactor = 1.0);
    virtual void updateScalarFieldSqrt(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex = -1, const float fFactor = 1.0);
    virtual void updateScalarField(const int nFieldIndex, const double dValue, const int nCellIndex);
    virtual void updateScalarField(const int nFieldIndex, const float fValue, const int nCellIndex);
    virtual void updateVectorField(const int nFieldIndex, const std::vector<double>& dVectorList, const int nSubmeshIndex = -1);
    virtual void updateVectorField(const int nFieldIndex, const std::vector<float>& fVectorList, const int nSubmeshIndex = -1);
    virtual void updateVectorField(const int nFieldIndex, const double* pVec, const int nCellIndex);
    virtual void updateVectorField(const int nFieldIndex, const float* pVec, const int nCellIndex);

    // Get Field Id
    virtual int getScalarFieldIndex(const std::string& sFieldName) const;
    virtual int getVectorFieldIndex(const std::string& sFieldName) const;

protected:
    std::vector<int> getMaskedList(const std::vector<int>& checkMaskList, const int nMask) const
    {
        // All Points
        std::vector<int> retList;
    
        // Calculate required size
        int nCount = 0;
        for (unsigned int nIndex = 0; nIndex < checkMaskList.size(); nIndex++)
            if ((nMask & checkMaskList[nIndex]) != 0)
                nCount++;
    
        // Resize
        retList.resize(nCount);

        // Assign values
        nCount = 0;
        for (unsigned int nIndex = 0; nIndex < checkMaskList.size(); nIndex++)
            if ((nMask & checkMaskList[nIndex]) != 0)
                retList[nCount++] = nIndex;
    
        return retList;
    }

protected:
    // Mesh 
    void* m_pRunTime;
    void* m_pMesh;
    
    int m_nTimeIndex;

    // Submesh
    int m_nSubmeshCount;
    std::vector<int> m_nCellSubmeshList;

    // Point / Face / Cell Proc Addressing
    std::vector<int> m_nPointProcAddressingList;
    std::vector<int> m_nFaceProcAddressingList;
    std::vector<int> m_nCellProcAddressingList;
    std::vector<int> m_nBoundaryProcAddressingList;

    // Point / Face / Cell BC Mask    
    std::vector<int> m_nPointMaskList;
    std::vector<int> m_nFaceMaskList;
    std::vector<int> m_nCellMaskList;

    // Point / Face / Cell Submesh Mapping
    std::vector<std::map<int,int>> m_nSubmeshPointMappingList;
    std::vector<std::map<int,int>> m_nSubmeshFaceMappingList;
    std::vector<std::map<int,int>> m_nSubmeshCellMappingList;

    // Valid Faces
    std::vector<bool> m_bValidFaceFlagList;

    // Scalar / Vector Fields
    std::vector<void *> m_pVolScalarFieldList;
    std::vector<std::string> m_sScalarFieldList;

    std::vector<void *> m_pVolVectorFieldList;
    std::vector<std::string> m_sVectorFieldList;
};

#endif
