#ifndef IMESHREADER_H
#define IMESHREADER_H

#include <ostream>
#include <vector>

class IMeshReader
{
public:
    // Destructor
    virtual ~IMeshReader() {}

    // Points
    virtual int getPointCount() const = 0;
    virtual void getPoint(const int nPointIndex, double *pXYZ) const = 0;
    virtual int getPointProcAddressing(const int nPointIndex) const = 0;

    // Faces
    virtual int getFaceCount() const = 0;
    virtual void getFaceAreaNormal(const int nFaceIndex, double* pN) const = 0;
    virtual void getFaceCenter(const int nFaceIndex, double* pXYZ) const = 0;
    virtual int getFaceOwner(const int nFaceIndex) const = 0;
    virtual int getFaceNeighbour(const int nFaceIndex) const = 0;
    virtual int getFaceBoundary(const int nFaceIndex) const = 0;    
    virtual int getFaceIndexInsideBoundary(const int nFaceIndex) const = 0;
    virtual int getFacePointCount(const int nFaceIndex) const = 0;
    virtual int getFaceProcAddressing(const int nFaceIndex) const = 0;
    virtual int getFacePointIndex(const int nFaceIndex, const int nIndex) const = 0;
    virtual std::vector<int> getFacePointIndexList(const int nFaceIndex) const = 0;
    virtual bool getFaceIsValid(const int nFaceIndex) const = 0;
    virtual int getCyclicFaceIndex(const int nFaceIndex) const = 0;
    virtual int getFaceId(const int nFaceIndex) const = 0;

    // Cells
    virtual int getCellCount() const = 0;
    virtual double getCellVolume(const int nCellIndex) const = 0;
    virtual void getCellCenter(const int nCellIndex, double* pXYZ) const = 0;
    virtual int getCellFaceCount(const int nCellIndex) const = 0;
    virtual int getCellFaceIndex(const int nCellIndex, const int nIndex) const = 0;
    virtual bool getCellFaceOwner(const int nCellIndex, const int nIndex) const = 0;
    virtual int getCellProcAddressing(const int nCellIndex) const = 0;
    virtual std::vector<int> getCellPointList(const int nCellIndex) const = 0;
    virtual int getCellValidFaceCount(const int nCellIndex) const = 0;
    virtual int getCellValidFaceIndex(const int nCellIndex, const int nIndex) const = 0;
    virtual double getCellScalar(const int nCellIndex, const int nFieldIndex) const = 0;
    virtual void getCellVector(const int nCellIndex, const int nFieldIndex, double *pVec) const = 0;

    // Boundary
    virtual int getBoundaryCount() const = 0;
    virtual std::string getBoundaryName(const int nBoundaryIndex) const = 0;
    virtual std::string getBoundaryType(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryFaceStart(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryFaceEnd(const int nBoundaryIndex) const = 0;
    virtual std::vector<int> getBoundaryFaceCellList(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryProcAddressing(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryProcessorRank(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryTag(const int nBoundaryIndex) const = 0;
    virtual int getBoundaryCyclicPairIndex(const int nBoundaryIndex) const = 0;
    virtual std::vector<int> getNeighbourCellList(const int nBoundaryIndex) const = 0;

    // Submesh
    virtual int getMask(const int nIndex) const = 0;
    virtual void initializeSubmesh() = 0;
    virtual int getSubmeshCount() const = 0;
    virtual int getCellSubmeshIndex(const int nCellIndex) const = 0;
    virtual std::vector<int> getSubmeshPointIndexList(const int nSubmeshIndex) const = 0;
    virtual std::vector<int> getSubmeshFaceIndexList(const int nSubmeshIndex) const = 0;
    virtual std::vector<int> getSubmeshCellIndexList(const int nSubmeshIndex) const = 0;
    virtual int getSubmeshPointIndex(const int nSubmeshIndex, const int nPointIndex) = 0;
    virtual int getSubmeshFaceIndex(const int nSubmeshIndex, const int nFaceIndex) = 0;
    virtual int getSubmeshCellIndex(const int nSubmeshIndex, const int nCellIndex) = 0;
    virtual std::vector<double> getCellDistanceFromBoundary(const int nSubmeshIndex, const std::vector<int> nBoundaryIndexList) const = 0;

    // ----------------------------- //
    // Scalar & Vector volume fields //
    // ----------------------------- //
    // Scalar & Vector Fields
    virtual int readScalarField(const std::string sFieldName) = 0;
    virtual int readVectorField(const std::string sFieldName) = 0;
    virtual int createScalarField(const std::string sFieldName) = 0;

    // Update Fields
    virtual void updateScalarField(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex = -1, const double dFactor = 1.0) = 0;
    virtual void updateScalarField(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex = -1, const float fFactor = 1.0) = 0;
    virtual void updateScalarFieldSqrt(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex = -1, const double dFactor = 1.0) = 0;
    virtual void updateScalarFieldSqrt(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex = -1, const float fFactor = 1.0) = 0;
    virtual void updateScalarField(const int nFieldIndex, const double dValue, const int nCellIndex) = 0;
    virtual void updateScalarField(const int nFieldIndex, const float fValue, const int nCellIndex) = 0;
    virtual void updateVectorField(const int nFieldIndex, const std::vector<double>& dVectorList, const int nSubmeshIndex = -1) = 0;
    virtual void updateVectorField(const int nFieldIndex, const std::vector<float>& fVectorList, const int nSubmeshIndex = -1) = 0;
    virtual void updateVectorField(const int nFieldIndex, const double* pVec, const int nCellIndex) = 0;
    virtual void updateVectorField(const int nFieldIndex, const float* pVec, const int nCellIndex) = 0;

    // Get Field Id
    virtual int getScalarFieldIndex(const std::string& sFieldName) const = 0;
    virtual int getVectorFieldIndex(const std::string& sFieldName) const = 0;

    // Locality Test
    virtual double getLocalityScore(const int nMethod = 0) const = 0;
};

#endif
