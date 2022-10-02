#include "polyMeshReaderOF.h"
#include "fvCFD.H"

// -------------------------------------------------------------------------- //
CPolyMeshReaderOF::CPolyMeshReaderOF(void* pRunTime)
{
    // RunTime
    m_pRunTime = pRunTime;
    Foam::Time& runTime = *(static_cast<Foam::Time *>(m_pRunTime));

    // Region Name
    Foam::word regionName(Foam::polyMesh::defaultRegion);

    // Read Mesh
    Foam::fvMesh* pPolyMesh = new Foam::fvMesh(Foam::IOobject(regionName, runTime.timeName(), runTime, Foam::IOobject::MUST_READ), false);
    Foam::fvMesh& mesh = *pPolyMesh;
    pPolyMesh->init(true);
    m_pMesh = pPolyMesh;

    // Read Parallel Information
    labelList emptyList;
    labelIOList pointProcAddressing = labelIOList(IOobject("pointProcAddressing", mesh.facesInstance(), mesh.meshSubDir, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), emptyList);
    m_nPointProcAddressingList.resize(pointProcAddressing.size());
    for (int i=0; i < pointProcAddressing.size(); i++)
        m_nPointProcAddressingList[i] = pointProcAddressing[i];

    labelIOList faceProcAddressing = labelIOList(IOobject("faceProcAddressing", mesh.facesInstance(), mesh.meshSubDir, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), emptyList);
    m_nFaceProcAddressingList.resize(faceProcAddressing.size());
    for (int i=0; i < faceProcAddressing.size(); i++)
        m_nFaceProcAddressingList[i] = faceProcAddressing[i];

    labelIOList cellProcAddressing = labelIOList(IOobject("cellProcAddressing", mesh.facesInstance(), mesh.meshSubDir, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), emptyList);
    m_nCellProcAddressingList.resize(cellProcAddressing.size());
    for (int i=0; i < cellProcAddressing.size(); i++)
        m_nCellProcAddressingList[i] = cellProcAddressing[i];

    labelIOList cellSubmeshList = labelIOList(IOobject("cellSubmesh", mesh.facesInstance(), mesh.meshSubDir, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), emptyList);
    m_nCellSubmeshList.resize(cellSubmeshList.size());
    for (int i=0; i < cellSubmeshList.size(); i++)
        m_nCellSubmeshList[i] = cellSubmeshList[i];

    labelIOList boundaryProcAddressing = labelIOList(IOobject("boundaryProcAddressing", mesh.facesInstance(), mesh.meshSubDir, mesh, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), emptyList);
    m_nBoundaryProcAddressingList.resize(boundaryProcAddressing.size());
    for (int i=0; i < boundaryProcAddressing.size(); i++)
        m_nBoundaryProcAddressingList[i] = boundaryProcAddressing[i];
    
    // Initialize faces to be valid (those with boundary different than empty)
    const int nFaceCount = pPolyMesh->faces().size();
    m_bValidFaceFlagList.resize(nFaceCount);
    for (int nFaceIndex = 0; nFaceIndex < nFaceCount; nFaceIndex++)
        m_bValidFaceFlagList[nFaceIndex] = true;

    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    for (int nBoundaryIndex = 0; nBoundaryIndex < bnd.size(); nBoundaryIndex++)
    {
        // Empty
        if (bnd.types()[nBoundaryIndex].compare("empty") == 0)
        {
            const int nMin = bnd.patchRanges()[nBoundaryIndex].min();
            const int nMax = bnd.patchRanges()[nBoundaryIndex].max();
            for (int nFaceIndex = nMin; nFaceIndex <= nMax; nFaceIndex++)
                m_bValidFaceFlagList[nFaceIndex] = false;
        }
    }
    
    // Initialize Submesh
    initializeSubmesh();
}

// -------------------------------------------------------------------------- //
CPolyMeshReaderOF::~CPolyMeshReaderOF()
{
    // Delete Fields
    for (unsigned int nIndex = 0; nIndex < m_pVolScalarFieldList.size(); nIndex++)
    {
        volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nIndex]);
        delete pField;
    }
    for (unsigned int nIndex = 0; nIndex < m_pVolVectorFieldList.size(); nIndex++)
    {
        volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nIndex]);
        delete pField;
    }

    // Delete fvMesh
    Foam::fvMesh* pFvMesh = static_cast<Foam::fvMesh *>(m_pMesh);
    if (pFvMesh)
        delete pFvMesh;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getPointCount() const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->points().size();
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::getPoint(const int nPointIndex, double *pXYZ) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const point& p = pPolyMesh->points()[nPointIndex];
    pXYZ[0] = p[0];
    pXYZ[1] = p[1];
    pXYZ[2] = p[2];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getPointProcAddressing(const int nPointIndex) const
{
    return (nPointIndex < static_cast<int>(m_nPointProcAddressingList.size())) ? m_nPointProcAddressingList[nPointIndex] : nPointIndex;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceCount() const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->faces().size();
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::getFaceAreaNormal(const int nFaceIndex, double* pN) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    Foam::vector n = pPolyMesh->faceAreas()[nFaceIndex];
    pN[0] = n[0];
    pN[1] = n[1];
    pN[2] = n[2];
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::getFaceCenter(const int nFaceIndex, double* pXYZ) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const point& p = pPolyMesh->faceCentres()[nFaceIndex];
    pXYZ[0] = p[0];
    pXYZ[1] = p[1];
    pXYZ[2] = p[2];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceOwner(const int nFaceIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->faceOwner()[nFaceIndex];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceNeighbour(const int nFaceIndex) const
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const labelList& neighbourList = pPolyMesh->faceNeighbour();

    // Interior Face
    if (nFaceIndex < neighbourList.size()) 
        return neighbourList[nFaceIndex];

    // Cyclic Face
    const int nCyclicFaceIndex = getCyclicFaceIndex(nFaceIndex);
    return (nCyclicFaceIndex == nFaceIndex) ? -1 : pPolyMesh->faceOwner()[nCyclicFaceIndex];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceBoundary(const int nFaceIndex) const
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);

    // Boundary Patches
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    return bnd.whichPatch(nFaceIndex);    
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceIndexInsideBoundary(const int nFaceIndex) const
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);

    // Boundary Patches
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    const int nBoundaryIndex = bnd.whichPatch(nFaceIndex);

    // Get Boundary Offset
    const int nOffset = (nBoundaryIndex == -1) ? 0 : bnd.patchStarts()[nBoundaryIndex];

    // Return with offset
    return nFaceIndex - nOffset;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFacePointCount(const int nFaceIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->faces()[nFaceIndex].size();
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceProcAddressing(const int nFaceIndex) const
{
    return (nFaceIndex < static_cast<int>(m_nFaceProcAddressingList.size())) ? m_nFaceProcAddressingList[nFaceIndex] : nFaceIndex;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFacePointIndex(const int nFaceIndex, const int nIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->faces()[nFaceIndex][nIndex];
}

// -------------------------------------------------------------------------- //
std::vector<int> CPolyMeshReaderOF::getFacePointIndexList(const int nFaceIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const face& checkFace = pPolyMesh->faces()[nFaceIndex];
    const int nPointCount = checkFace.size();
    std::vector<int> pointList;
    pointList.resize(nPointCount);
    for (int nIndex=0; nIndex < nPointCount; nIndex++)
        pointList[nIndex] = checkFace[nIndex];
    return pointList;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCyclicFaceIndex(const int nFaceIndex) const
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);

    // Boundary Patches
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    const int nBoundaryIndex = bnd.whichPatch(nFaceIndex);

    // Not Coupled
    if (nBoundaryIndex == -1 || !bnd[nBoundaryIndex].coupled())
        return nFaceIndex;

    // Neighbor Boundary Index
    const int nNeighborBoundaryIndex = bnd[nBoundaryIndex].neighbPolyPatchID();

    // Invalid Index
    if (nNeighborBoundaryIndex == -1) 
        return nFaceIndex;
        
    // Return cyclic face index
    return nFaceIndex - bnd.patchStarts()[nBoundaryIndex] + bnd.patchStarts()[nNeighborBoundaryIndex];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getFaceId(const int nFaceIndex) const
{
    // Boundary Tag
    const int nTag = getBoundaryTag(getFaceBoundary(nFaceIndex));

    // No Tag / 'Simple' MPI boundary (i.e tag = -1/1)
    if (nTag >= -1 && nTag <= 1)
        return copysign(getFaceProcAddressing(nFaceIndex), 1);
    // Otherwise - return face index within boundary
    return getFaceIndexInsideBoundary(nFaceIndex);   
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellCount() const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->cells().size();
}

// -------------------------------------------------------------------------- //
double CPolyMeshReaderOF::getCellVolume(const int nCellIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->cellVolumes()[nCellIndex];    
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::getCellCenter(const int nCellIndex, double* pXYZ) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const point& p = pPolyMesh->cellCentres()[nCellIndex];
    pXYZ[0] = p[0];
    pXYZ[1] = p[1];
    pXYZ[2] = p[2];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellFaceCount(const int nCellIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->cells()[nCellIndex].nFaces();
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellFaceIndex(const int nCellIndex, const int nIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->cells()[nCellIndex][nIndex];
}

// -------------------------------------------------------------------------- //
bool CPolyMeshReaderOF::getCellFaceOwner(const int nCellIndex, const int nIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    return pPolyMesh->faceOwner()[pPolyMesh->cells()[nCellIndex][nIndex]] == nCellIndex;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellProcAddressing(const int nCellIndex) const
{
    return (nCellIndex < static_cast<int>(m_nCellProcAddressingList.size())) ? m_nCellProcAddressingList[nCellIndex] : nCellIndex;
}

// -------------------------------------------------------------------------- //
std::vector<int> CPolyMeshReaderOF::getCellPointList(const int nCellIndex) const
{
    std::vector<int> retList;
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const cell& c = pPolyMesh->cells()[nCellIndex];
    const int nFaceCount = c.size();
    for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
    {
        const face& f = pPolyMesh->faces()[c[nFIndex]];
        const int nPointCount = f.size();
        for (int nPIndex = 0; nPIndex < nPointCount; nPIndex++)
        {
            int bNewPoint = true;
            for (unsigned int nIndex = 0; nIndex < retList.size() && bNewPoint; nIndex++)
                bNewPoint = (retList[nIndex] != f[nPIndex]);
            if (bNewPoint)
                retList.push_back(f[nPIndex]);
        }
    }
    return retList;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellValidFaceCount(const int nCellIndex) const
{
    int nValidFaceCount = 0;
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const cell& c = pPolyMesh->cells()[nCellIndex];
    for (int nFIndex = 0; nFIndex < c.size(); nFIndex++)
        if (m_bValidFaceFlagList[c[nFIndex]])
            nValidFaceCount++;
    return nValidFaceCount;    
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getCellValidFaceIndex(const int nCellIndex, const int nIndex) const
{
    int nValidFaceIndex = 0;
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const cell& c = pPolyMesh->cells()[nCellIndex];
    for (int nFIndex = 0; nFIndex < c.size(); nFIndex++)
        if (m_bValidFaceFlagList[c[nFIndex]])
        {
            // Desired face
            if (nValidFaceIndex == nIndex) 
                return c[nFIndex];
            nValidFaceIndex++;
        }
    return -1;
}

// -------------------------------------------------------------------------- //
double CPolyMeshReaderOF::getCellScalar(const int nCellIndex, const int nFieldIndex) const
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);
    return field[nCellIndex];
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::getCellVector(const int nCellIndex, const int nFieldIndex, double* pVec) const
{
    volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nFieldIndex]);
    volVectorField& field = static_cast<volVectorField&>(*pField);
    for (int nD = 0; nD <= 2; nD++)
        pVec[nD] = field[nCellIndex][nD];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryCount() const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    return bnd.patchStarts().size();
}

// -------------------------------------------------------------------------- //
std::string CPolyMeshReaderOF::getBoundaryName(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    return bnd.names()[nBoundaryIndex];
}

// -------------------------------------------------------------------------- //
std::string CPolyMeshReaderOF::getBoundaryType(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    return bnd.types()[nBoundaryIndex];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryFaceStart(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    return bnd.patchStarts()[nBoundaryIndex];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryFaceEnd(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    
    return bnd.patchRanges()[nBoundaryIndex].max();
}

// -------------------------------------------------------------------------- //
std::vector<int> CPolyMeshReaderOF::getBoundaryFaceCellList(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    const labelUList faceCells = bnd[nBoundaryIndex].faceCells();
    std::vector<int> nFaceCellIndexList;
    const int nFaceCellCount = static_cast<int>(faceCells.size());
    nFaceCellIndexList.resize(nFaceCellCount);
    for (int i = 0; i < nFaceCellCount; i++)
        nFaceCellIndexList[i] = static_cast<int>(faceCells[i]);
    return nFaceCellIndexList;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryProcAddressing(const int nBoundaryIndex) const
{
    return (nBoundaryIndex < static_cast<int>(m_nBoundaryProcAddressingList.size())) ? m_nBoundaryProcAddressingList[nBoundaryIndex] : nBoundaryIndex;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryProcessorRank(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();

    // Non Process Boundary
    if (nBoundaryIndex < bnd.nNonProcessor()) return -1;

    // Cast
    const polyPatch* pPolyPatch = &bnd[nBoundaryIndex];
    const processorPolyPatch* pProcPolyPatch = static_cast<const processorPolyPatch*>(pPolyPatch);
    return pProcPolyPatch->neighbProcNo();
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryTag(const int nBoundaryIndex) const
{
    // Boundary
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();

    // Not Processor
    if (nBoundaryIndex < bnd.nNonProcessor())
        return 0;

    // Cast to processorPolyPatch
    const polyPatch* pPolyPatch = &bnd[nBoundaryIndex];
    const processorPolyPatch* pProcPolyPatch = static_cast<const processorPolyPatch*>(pPolyPatch);
    const int nSign = 2 * pProcPolyPatch->owner() - 1;
    return pProcPolyPatch->tag() * nSign;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getBoundaryCyclicPairIndex(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    
    return bnd[nBoundaryIndex].coupled() ? -1 : bnd[nBoundaryIndex].neighbPolyPatchID();
}


// -------------------------------------------------------------------------- //
std::vector<int> CPolyMeshReaderOF::getNeighbourCellList(const int nBoundaryIndex) const
{
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const polyBoundaryMesh& bnd = pPolyMesh->boundaryMesh();
    const labelUList& nbrCells = bnd[nBoundaryIndex].nbrCells();
    std::vector<int> nNeighbourCellList;
    nNeighbourCellList.resize(nbrCells.size());
    for (int nI = 0; nI < nbrCells.size(); nI++)
        nNeighbourCellList[nI] = nbrCells[nI];
    return nNeighbourCellList;
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::initializeSubmesh()
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const faceList& faces = pPolyMesh->faces();
    const cellList& cells = pPolyMesh->cells();

    // Get Geometry Sizes
    const unsigned int nPointCount = getPointCount();
    const unsigned int nFaceCount = faces.size();
    const unsigned int nCellCount = cells.size();
    const unsigned int nBoundaryCount = getBoundaryCount();

    // Submesh not read from dictionary
    if (m_nCellSubmeshList.size() != nCellCount)
    {
        // Points Interior Flag
        std::vector<bool> bPointIsInterior(pPolyMesh->points().size(), true);

        // Resize
        m_nCellSubmeshList.resize(cells.size(), 1);

        // Get all the bc points
        for (unsigned int nBoundaryIndex = 0; nBoundaryIndex < nBoundaryCount; nBoundaryIndex++)
        {
            // Ignore 'empty' boundary type
            if (getBoundaryType(nBoundaryIndex).compare("empty") == 0)
                continue;
            
            // Get Faces
            const int nFaceStart = getBoundaryFaceStart(nBoundaryIndex);
            const int nFaceEnd   = getBoundaryFaceEnd(nBoundaryIndex);

            // Mark their points as BC
            for (int nFaceIndex = nFaceStart; nFaceIndex <= nFaceEnd; nFaceIndex++)
                for (int nPIndex = 0; nPIndex < faces[nFaceIndex].size(); nPIndex++)
                    bPointIsInterior[faces[nFaceIndex][nPIndex]] = false;
        }

        // Check all the cells
        for (unsigned int nCellIndex = 0; nCellIndex < nCellCount; nCellIndex++)
        {
            bool bIsInteriorCell = true;
            for (int nFIndex = 0; nFIndex < cells[nCellIndex].size() && bIsInteriorCell; nFIndex++)
            {
                const int nFaceIndex = cells[nCellIndex][nFIndex];
                for (int nPIndex = 0; nPIndex < faces[nFaceIndex].size() && bIsInteriorCell; nPIndex++)
                    bIsInteriorCell = bPointIsInterior[faces[nFaceIndex][nPIndex]];
            }
            m_nCellSubmeshList[nCellIndex] = bIsInteriorCell ? 1 : 0;
        }
    }

    // Initialize Mask Arrays
    m_nPointMaskList.resize(nPointCount, 0);
    m_nFaceMaskList.resize(nFaceCount, 0);
    m_nCellMaskList.resize(nCellCount, 0);

    // Define Cells / Faces / Points Masks
    m_nSubmeshCount = 0;
    for (unsigned int nCellIndex = 0; nCellIndex < nCellCount; nCellIndex++)
    {
        // Cell Submesh
        const int nCellSubmesh = m_nCellSubmeshList[nCellIndex];
        if (nCellSubmesh >= m_nSubmeshCount)
            m_nSubmeshCount = nCellSubmesh + 1;

        // Cell Mask
        const int nMask = getMask(nCellSubmesh);
        m_nCellMaskList[nCellIndex] |= nMask;

        // Cell
        const cell& c = cells[nCellIndex];
        

        // Faces
        for (int nFIndex = 0; nFIndex < c.size(); nFIndex++)
        {
            // Face
            const int nFaceIndex = c[nFIndex];
            const face& f = faces[nFaceIndex];

            // Face Mask
            m_nFaceMaskList[nFaceIndex] |= nMask;

            //Points Mask
            for (int nPIndex = 0; nPIndex < f.size(); nPIndex++)
                m_nPointMaskList[f[nPIndex]] |= nMask;
        }
    }

    // Initialize Submesh Mapping
    m_nSubmeshPointMappingList.resize(m_nSubmeshCount);
    m_nSubmeshFaceMappingList.resize(m_nSubmeshCount);
    m_nSubmeshCellMappingList.resize(m_nSubmeshCount);
    for (int nSubmeshIndex = 0; nSubmeshIndex < m_nSubmeshCount; nSubmeshIndex++)
    {
        const std::vector<int> nPointIndexList = getSubmeshPointIndexList(nSubmeshIndex);
        for (unsigned int nPIndex = 0; nPIndex < nPointIndexList.size(); nPIndex++)
            m_nSubmeshPointMappingList[nSubmeshIndex].insert(std::pair<int,int>(nPointIndexList[nPIndex], nPIndex));

        const std::vector<int> nFaceIndexList = getSubmeshFaceIndexList(nSubmeshIndex);
        for (unsigned int nFIndex = 0; nFIndex < nFaceIndexList.size(); nFIndex++)
            m_nSubmeshFaceMappingList[nSubmeshIndex].insert(std::pair<int,int>(nFaceIndexList[nFIndex], nFIndex));

        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            m_nSubmeshCellMappingList[nSubmeshIndex].insert(std::pair<int,int>(nCellIndexList[nCIndex], nCIndex));
    }
}

// -------------------------------------------------------------------------- //
std::vector<double> CPolyMeshReaderOF::getCellDistanceFromBoundary(const int nSubmeshIndex, const std::vector<int> nBoundaryIndexList) const
{
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);

    // Get Cell Index List
    std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
    const int nCellCount = nCellIndexList.size();

    // Cell Distance List
    std::vector<double> dCellDistanceList;
    dCellDistanceList.resize(nCellCount, INT64_MAX);

    double dDistanceVec[3];
    double dDistance2;
    for (unsigned int nBIndex = 0; nBIndex < nBoundaryIndexList.size(); nBIndex++)
    {
        const int nBoundaryIndex = nBoundaryIndexList[nBIndex];

        // Get Boundary Face Range
        const int nFaceStart = getBoundaryFaceStart(nBoundaryIndex);
        const int nFaceEnd = getBoundaryFaceEnd(nBoundaryIndex);

        // Get Face Normal
        const int nFaceCount = nFaceEnd - nFaceStart + 1;
        std::vector<double> dFaceNormalList;
        dFaceNormalList.resize(nFaceCount * 3);
        for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
        {
            const int nFaceIndex = nFIndex + nFaceStart;
            Foam::vector n = pPolyMesh->faceAreas()[nFaceIndex];
            const double dMag = std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            for (int nD = 0; nD <= 2; nD++)
                dFaceNormalList[nFIndex * 3 + nD] = n[nD] / dMag;
        }

        // Calculate distance of each cell        
        for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
        {
            // Cell Index
            const int nCellIndex = nCellIndexList[nCIndex];

            // Cell Center
            const point &dCellCenter = pPolyMesh->cellCentres()[nCellIndex];

            // Check Each Face
            for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
            {
                // Face Index
                const int nFaceIndex = nFIndex + nFaceStart;

                // Face Center
                const point &dFaceCenter = pPolyMesh->faceCentres()[nFaceIndex];

                // Distance Vec                
                for (int nD = 0; nD <= 2; nD++)
                    dDistanceVec[nD] = (dFaceNormalList[nFIndex * 3 + nD] * (dFaceCenter[nD] - dCellCenter[nD]));
                dDistance2 = dDistanceVec[0] * dDistanceVec[0] + dDistanceVec[1] * dDistanceVec[1] + dDistanceVec[2] * dDistanceVec[2];                

                // Smaller Distance
                if (dDistance2 < dCellDistanceList[nCIndex])
                    dCellDistanceList[nCIndex] = dDistance2;
            }

            // Sqrt Distance
            dCellDistanceList[nCIndex] = std::sqrt(dCellDistanceList[nCIndex]);
        }
    }

    // Return the distance of each cell
    return dCellDistanceList;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::readScalarField(const std::string sFieldName)
{   
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;

    // Field Already Loaded
    int nFieldIndex = getScalarFieldIndex(sFieldName);
    if (nFieldIndex != -1) return nFieldIndex;
            
    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const fvMesh& mesh = static_cast<fvMesh&>(*pPolyMesh);

    // Read Field
    volScalarField* pField = new volScalarField(
        IOobject
        (
            sFieldName.c_str(),
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ), mesh
    );

    // Append to lists
    nFieldIndex = m_pVolScalarFieldList.size();          
    m_sScalarFieldList.push_back(sFieldName);
    m_pVolScalarFieldList.push_back(pField);

    // Return Field Index
    return nFieldIndex;    
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::readVectorField(const std::string sFieldName)
{
    Foam::Time* pRunTime = static_cast<Foam::Time *>(m_pRunTime);
    Foam::Time& runTime = *pRunTime;

    // Field Already Loaded
    int nFieldIndex = getVectorFieldIndex(sFieldName);
    if (nFieldIndex != -1) return nFieldIndex;

    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const fvMesh& mesh = static_cast<fvMesh&>(*pPolyMesh);

    // Data
    volVectorField* pField = new volVectorField(
        IOobject
        (
            sFieldName,
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ), mesh
    );

    // Append to lists
    nFieldIndex = m_pVolVectorFieldList.size();
    m_sVectorFieldList.push_back(sFieldName);
    m_pVolVectorFieldList.push_back(pField);

    // Return Field Index
    return nFieldIndex;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::createScalarField(const std::string sFieldName)
{
    // Field Already Loaded
    int nFieldIndex = getScalarFieldIndex(sFieldName);
    if (nFieldIndex != -1) return nFieldIndex;

    // PolyMesh
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const fvMesh& mesh = static_cast<fvMesh&>(*pPolyMesh);

    // Create Scalar Field    
    volScalarField* pField = new volScalarField(
             IOobject
             (
                 sFieldName,
                 mesh.time().timeName(),
                 mesh,
                 IOobject::NO_READ,
                 IOobject::AUTO_WRITE                 
             ),
             mesh,
             dimensionedScalar(dimless, Zero),
             zeroGradientFvPatchScalarField::typeName
         );

    // Append to lists
    nFieldIndex = m_pVolScalarFieldList.size();          
    m_sScalarFieldList.push_back(sFieldName);
    m_pVolScalarFieldList.push_back(pField);

    // Return Field Index
    return nFieldIndex;    
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarField(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex, const double dFactor)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);

    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            field[nCellIndex] = dScalarList[nCellIndex] * dFactor;
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            field[nCellIndexList[nCIndex]] = dScalarList[nCIndex] * dFactor;
    }
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarField(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex, const float fFactor)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);

    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            field[nCellIndex] = fScalarList[nCellIndex] * fFactor;
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            field[nCellIndexList[nCIndex]] = fScalarList[nCIndex] * fFactor;
    }    
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarFieldSqrt(const int nFieldIndex, const std::vector<double>& dScalarList, const int nSubmeshIndex, const double dFactor)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);

    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            field[nCellIndex] = std::sqrt(dScalarList[nCellIndex] * dFactor);
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            field[nCellIndexList[nCIndex]] = std::sqrt(dScalarList[nCIndex] * dFactor);
    }
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarFieldSqrt(const int nFieldIndex, const std::vector<float>& fScalarList, const int nSubmeshIndex, const float fFactor)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);

    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            field[nCellIndex] = std::sqrt(fScalarList[nCellIndex] * fFactor);
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            field[nCellIndexList[nCIndex]] = std::sqrt(fScalarList[nCIndex] * fFactor);
    }  
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarField(const int nFieldIndex, const double dValue, const int nCellIndex)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);
    field[nCellIndex] = dValue;
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateScalarField(const int nFieldIndex, const float fValue, const int nCellIndex)
{
    volScalarField* pField = static_cast<volScalarField *>(m_pVolScalarFieldList[nFieldIndex]);
    volScalarField& field = static_cast<volScalarField&>(*pField);
    field[nCellIndex] = fValue;
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateVectorField(const int nFieldIndex, const std::vector<double>& dVectorList, const int nSubmeshIndex)
{
    volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nFieldIndex]);
    volVectorField& field = static_cast<volVectorField&>(*pField);

    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            for (int nD = 0; nD <= 2; nD++)
                field[nCellIndex][nD] = dVectorList[nCellIndex * 3 + nD];
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            for (int nD = 0; nD <= 2; nD++)
                field[nCellIndexList[nCIndex]][nD] = dVectorList[nCIndex * 3 + nD];
    }
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateVectorField(const int nFieldIndex, const std::vector<float>& fVectorList, const int nSubmeshIndex)
{
    volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nFieldIndex]);
    volVectorField& field = static_cast<volVectorField&>(*pField);
    
    // Update all cells
    if (nSubmeshIndex == -1)
        for (int nCellIndex = 0; nCellIndex < field.size(); nCellIndex++)
            for (int nD = 0; nD <= 2; nD++)
                field[nCellIndex][nD] = fVectorList[nCellIndex * 3 + nD];
    else
    {
        const std::vector<int> nCellIndexList = getSubmeshCellIndexList(nSubmeshIndex);
        for (unsigned int nCIndex = 0; nCIndex < nCellIndexList.size(); nCIndex++)
            for (int nD = 0; nD <= 2; nD++)
                field[nCellIndexList[nCIndex]][nD] = fVectorList[nCIndex * 3 + nD];
    }
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateVectorField(const int nFieldIndex, const double* pVec, const int nCellIndex)
{
    volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nFieldIndex]);
    volVectorField& field = static_cast<volVectorField&>(*pField);
    for (int nD = 0; nD <= 2; nD++)
        field[nCellIndex][nD] = pVec[nD];
}

// -------------------------------------------------------------------------- //
void CPolyMeshReaderOF::updateVectorField(const int nFieldIndex, const float* pVec, const int nCellIndex)
{
    volVectorField* pField = static_cast<volVectorField *>(m_pVolVectorFieldList[nFieldIndex]);
    volVectorField& field = static_cast<volVectorField&>(*pField);
    for (int nD = 0; nD <= 2; nD++)
        field[nCellIndex][nD] = pVec[nD];
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getScalarFieldIndex(const std::string& sFieldName) const
{
    for (unsigned int nIndex = 0; nIndex < m_sScalarFieldList.size(); nIndex++)
        if (sFieldName.compare(m_sScalarFieldList[nIndex]) == 0)
            return nIndex;
    return -1;
}

// -------------------------------------------------------------------------- //
int CPolyMeshReaderOF::getVectorFieldIndex(const std::string& sFieldName) const
{
    // Field Already Loaded
    for (unsigned int nIndex = 0; nIndex < m_sVectorFieldList.size(); nIndex++)
        if (sFieldName.compare(m_sVectorFieldList[nIndex]) == 0)
            return nIndex;
    return -1;
}


// -------------------------------------------------------------------------- //
double CPolyMeshReaderOF::getLocalityScore(const int nMethod) const
{
    // Get average distance between each cell and its neighbours
    double dScore = 0;
    
    polyMesh* pPolyMesh = static_cast<polyMesh *>(m_pMesh);
    const cellList& cells = pPolyMesh->cells();
    const int nCellCount = cells.size();

    // Neighbour Cells
    const labelList& neighbourList = pPolyMesh->faceNeighbour();
    
    for (int nCellIndex = 0; nCellIndex < nCellCount; nCellIndex++)
    {
        // Cell Face Count
        const int nCellFaceCount = cells[nCellIndex].nFaces();

        double dCellScore = 0;
        int nNeighbourCount = 0;
        for (int nIndex = 0; nIndex < nCellFaceCount; nIndex++)
        {
            const int nFaceIndex = cells[nCellIndex][nIndex];
            if (nFaceIndex < neighbourList.size())
            {
                const int nNeighbourCellIndex = neighbourList[nFaceIndex];
                const double dDist = 1.0 * copysign(nNeighbourCellIndex - nCellIndex, 1) / nCellCount;
                dCellScore += (nMethod == 0) ? dDist : dDist * dDist;
                nNeighbourCount++;
            }
        }
        dScore += dCellScore / nNeighbourCount;
    }

    return dScore / nCellCount;

}
