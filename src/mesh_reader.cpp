#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <sstream>

#include "api/fastmesh.h"
#include "api/mpi_env.h"
#include "api/polyMeshReaderOF.h"

using namespace std;
using namespace fastmesh;

//***************************************************************************************************
//									PUBLIC FUNCTION
//***************************************************************************************************

//***************************************************************************************************
// Initialize simulation:
//		- MPI framework
//		- read LFM mesh files
//		- setup mesh parameters
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: initialize()
{	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// - Get global rank/processCount
	// - Get local rank/processCount (and create communicator)
	m_MPI_Env.initialize( true );

	// Read mesh files (one per MPI rank)
	if( m_MPI_Env.is_master()) cout << "read mesh files..." << endl;
	readMeshFiles();
	
	// Initialize Mesh Parameters
	for (unsigned int nSubmeshIndex = 0; nSubmeshIndex < m_SubmeshList.size(); nSubmeshIndex++)
		setParams(nSubmeshIndex);

	// Calculate vector to neighboring cells
	for (unsigned int nSubmeshIndex = 0; nSubmeshIndex < m_SubmeshList.size(); nSubmeshIndex++)
		calculateVNeighbour(nSubmeshIndex);

	MPI_Barrier(MPI_COMM_WORLD);
}

//***************************************************************************************************
//***************************************************************************************************
//									PRIVATE FUNCTIONS
//***************************************************************************************************
//***************************************************************************************************

//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: readMeshFiles()
{
	// Read polyMesh
	m_pMeshReader = new CPolyMeshReaderOF(m_pRunTime->getRunPointer());

	// Submesh Count
	const int nSubmeshCount = m_pMeshReader->getSubmeshCount();

	// Nodes / Faces / Cells
	const int nNodeCount = m_pMeshReader->getPointCount();
	const int nFaceCount = m_pMeshReader->getFaceCount();
	const int nCellCount = m_pMeshReader->getCellCount();

	// Resize gmsh_parts
	m_SubmeshList.resize(nSubmeshCount);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Read submesh
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	for(int nSubmeshIndex=0; nSubmeshIndex < nSubmeshCount; nSubmeshIndex++)
	{
		// Id
		m_SubmeshList[nSubmeshIndex].id = nSubmeshIndex;

		// Extract Nodes, Faces, Cells
		readMeshNodes(nSubmeshIndex);
		readMeshFaces(nSubmeshIndex);
		readMeshCells(nSubmeshIndex);
	}

}

// -------------------------------------------------------------------------- //
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT>::readMeshNodes(const int nSubmeshIndex)
{
	// Submesh
	gmsh_mesh *pSubMesh = &m_SubmeshList[nSubmeshIndex];

	// Nodes
	const std::vector<int> nNodeIndexList = m_pMeshReader->getSubmeshPointIndexList(nSubmeshIndex);
	const int nNodeCount = nNodeIndexList.size();

	// Allocate Vectors
	pSubMesh->nodes.resize(nNodeCount);

	// Extract Info
	for (size_t nIndex = 0; nIndex < nNodeCount; nIndex++)
	{
		const int nPointIndex = nNodeIndexList[nIndex];
		pSubMesh->nodes[nIndex].id = nPointIndex;
		m_pMeshReader->getPoint(nPointIndex, pSubMesh->nodes[nIndex].xn);
	}
}

// -------------------------------------------------------------------------- //
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT>::readMeshFaces(const int nSubmeshIndex)
{
	// Submesh
	gmsh_mesh *pSubMesh = &m_SubmeshList[nSubmeshIndex];

	// Faces
	const std::vector<int> nFaceIndexList = m_pMeshReader->getSubmeshFaceIndexList(nSubmeshIndex);
	const int nFaceCount = nFaceIndexList.size();

	// Max Node Per Face
	int nMaxNodePerFace = 0;
	for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
	{
		const int nFaceNodeCount = m_pMeshReader->getFacePointCount(nFaceIndexList[nFIndex]);
		if (nFaceNodeCount > nMaxNodePerFace)
			nMaxNodePerFace = nFaceNodeCount;
	}

	// Allocate Vectors
	pSubMesh->faces.resize(nFaceCount);
	pSubMesh->face2cell.resize(nFaceCount, std::vector<int>(2, -1));
	pSubMesh->face2node.resize(nFaceCount, std::vector<int>(nMaxNodePerFace, -1));

	// Extract Info
	for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
	{
		// Face Index
		const int nFaceIndex = nFaceIndexList[nFIndex];

		// Remember Id to meshReader
		pSubMesh->faces[nFIndex].id = nFaceIndex;

		// Face Center
		m_pMeshReader->getFaceCenter(nFaceIndex, pSubMesh->faces[nFIndex].xf);

		// Get Adjacent Cells
		pSubMesh->face2cell[nFIndex][0] = m_pMeshReader->getFaceOwner(nFaceIndex);
		pSubMesh->face2cell[nFIndex][1] = m_pMeshReader->getFaceNeighbour(nFaceIndex);

		// Get Nodes
		const std::vector<int> nFaceNodeIndexList = m_pMeshReader->getFacePointIndexList(nFaceIndex);
		const int nFaceNodeCount = nFaceNodeIndexList.size();
		for (int nNIndex = 0; nNIndex < nFaceNodeCount && nNIndex < nMaxNodePerFace; nNIndex++)
			pSubMesh->face2node[nFIndex][nNIndex] = nFaceNodeIndexList[nNIndex];
	}
}

// -------------------------------------------------------------------------- //
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT>::readMeshCells(const int nSubmeshIndex)
{
	// Submesh
	gmsh_mesh *pSubMesh = &m_SubmeshList[nSubmeshIndex];

	// Cells
	const std::vector<int> nCellIndexList = m_pMeshReader->getSubmeshCellIndexList(nSubmeshIndex);
	const int nCellCount = nCellIndexList.size();

	int nMaxCellFaceCount = 0;
	int nMaxCellNodeCount = 0;
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
	{
		// Cell Index
		const int nCellIndex = nCellIndexList[nCIndex];

		// Get Valid Face Count (without 'empty' bc)
		const int nCellFaceCount = m_pMeshReader->getCellValidFaceCount(nCellIndex);
		if (nCellFaceCount > nMaxCellFaceCount)
			nMaxCellFaceCount = nCellFaceCount;

		// Get Node Count
		const std::vector<int> nCellNodeIndexList = m_pMeshReader->getCellPointList(nCellIndex);
		const int nCellNodeCount = nCellNodeIndexList.size();
		if (nCellNodeCount > nMaxCellNodeCount)
			nMaxCellNodeCount = nCellNodeCount;
	}
	pSubMesh->faces_per_cell = nMaxCellFaceCount;
	pSubMesh->nodes_per_cell = nMaxCellNodeCount;

	// Allocate Vectors
	pSubMesh->cells     .resize( nCellCount );
	pSubMesh->cells_vol .resize( nCellCount );
	pSubMesh->cell2node .resize( nCellCount, std::vector<int>            ( nMaxCellNodeCount, -1 ) );
	pSubMesh->cell2face .resize( nCellCount, std::vector<int>            ( nMaxCellFaceCount, -1 ) );
	pSubMesh->cell2neigh.resize( nCellCount, std::vector<t_gmsh_neighbor>( nMaxCellFaceCount) );

	// Extract Info
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
	{
		// Cell Index
		const int nCellIndex = nCellIndexList[nCIndex];
		pSubMesh->cells[nCIndex].id = nCellIndex;

		// Cell Volume
		pSubMesh->cells_vol[nCIndex] = m_pMeshReader->getCellVolume(nCellIndex);

		// Cell Center
		m_pMeshReader->getCellCenter(nCellIndex, pSubMesh->cells[nCIndex].xc);

		// Cell Nodes
		const std::vector<int> nCellNodeIndexList = m_pMeshReader->getCellPointList(nCellIndex);
		const int nCellNodeCount = nCellNodeIndexList.size();
		pSubMesh->cells[nCIndex].nodeCount = nCellNodeCount;
		for (int nNIndex = 0; nNIndex < nCellNodeCount; nNIndex++)
			pSubMesh->cell2node[nCIndex][nNIndex] = m_pMeshReader->getSubmeshPointIndex(nSubmeshIndex, nCellNodeIndexList[nNIndex]);

		// Cell Faces (ignore faces lying on 'empty' boundary condition)
		const int nCellFaceCount = m_pMeshReader->getCellValidFaceCount(nCellIndex);
		pSubMesh->cells[nCIndex].faceCount = nCellFaceCount;
		for (int nFIndex = 0; nFIndex < nCellFaceCount; nFIndex++)
		{
			// Face Index
			const int nFaceIndex = m_pMeshReader->getCellValidFaceIndex(nCellIndex, nFIndex);			
			pSubMesh->cell2face[nCIndex][nFIndex] = m_pMeshReader->getSubmeshFaceIndex(nSubmeshIndex, nFaceIndex);

			const int nFaceCellOwner = m_pMeshReader->getFaceOwner(nFaceIndex);
			const int nFaceCellNeighbour = m_pMeshReader->getFaceNeighbour(nFaceIndex);

			// Neighbour Cell
			const int nNeighbourCell = (nFaceCellOwner == nCellIndex) ? nFaceCellNeighbour : nFaceCellOwner;			

			// Neighbour submesh
			const int nNeighbourSubmesh = (nNeighbourCell == -1) ? -1 : m_pMeshReader->getCellSubmeshIndex(nNeighbourCell);
			const int nNeighbourCellId = (nNeighbourCell == -1) ? -1 : m_pMeshReader->getSubmeshCellIndex(nNeighbourSubmesh, nNeighbourCell);

			// Boundary Index
			const int nBoundaryIndex = m_pMeshReader->getFaceBoundary(nFaceIndex);

			// Boundary Processor Rank
			const int nBoundaryRankId = (nBoundaryIndex == -1) ? -1 : m_pMeshReader->getBoundaryProcessorRank(nBoundaryIndex);

			// Neighbour Type
			neigh_type eNeighbourType = (nNeighbourCell != -1) ? (nNeighbourSubmesh == nSubmeshIndex) ? CELL_REGULAR : CELL_SUBMESH :
										(nBoundaryRankId == -1) ? CELL_NONE : CELL_MPI;

			// Note: neighbour.sm contains: sm (if local) / remote rank
			const int nStoredSM = (nBoundaryRankId == -1) ? nNeighbourSubmesh : nBoundaryRankId;

			// Set Neighbour Info
			pSubMesh->cell2neigh[nCIndex][nFIndex].type = eNeighbourType;	
			pSubMesh->cell2neigh[nCIndex][nFIndex].sm = nStoredSM;
			pSubMesh->cell2neigh[nCIndex][nFIndex].id = nNeighbourCellId;
			pSubMesh->cell2neigh[nCIndex][nFIndex].bound = nBoundaryIndex;
		}
	}
		
}

// -------------------------------------------------------------------------- //
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT>::setParams(const int nSubmeshIndex)
{
	// Get SubMesh
	gmsh_mesh *pSubMesh = &m_SubmeshList[nSubmeshIndex];

	// Evaluate Face Area & Center & Normal
	const int nFaceCount = pSubMesh->faces.size();
	double dP0[3], dP1[3];
	double dAreaNormal[3] = {0, 0, 0};
	double dTangVect[3] = {0, 0, 0};
	pSubMesh->faces_area.resize(nFaceCount);
	pSubMesh->v_norm.resize(nFaceCount, std::vector<double>(DIM_CNT));
	pSubMesh->v_tang.resize(nFaceCount, std::vector<double>(DIM_CNT));
	for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
	{
			// Face Index
			const int nFaceIndex = pSubMesh->faces[nFIndex].id;

			// Face Area Normal
			m_pMeshReader->getFaceAreaNormal(nFaceIndex, dAreaNormal);
		
			// Face Area
			const double dFaceArea = sqrt(pow(dAreaNormal[0], 2) + 
										  pow(dAreaNormal[1], 2) + 
										  pow(dAreaNormal[2], 2));

			pSubMesh->faces_area[nFIndex] = dFaceArea;


			// Calculate Tangential Direction
			// Try cross operator with [0, 0, 1] vector.
			// If its 0 (same direction as normal), multiply by [0, 1, 0]
			if (dAreaNormal[0] != 0 || dAreaNormal[1] != 0)
			{
				dTangVect[0] = -dAreaNormal[1];
				dTangVect[1] = dAreaNormal[0];
				dTangVect[2] = 0;
			}
			else
			{
				dTangVect[0] = dAreaNormal[2];
				dTangVect[1] = 0;
				dTangVect[2] = -dAreaNormal[0];
			}

			// Store Normal & Tangentail Normal
			for (int iDim = 0; iDim < DIM_CNT; iDim++)
			{
				pSubMesh->v_norm[nFIndex][iDim] = dAreaNormal[iDim];
				pSubMesh->v_tang[nFIndex][iDim] = dTangVect[iDim];
			}				
	}

	// Calculate the face normal sign for each face
	const int nCellCount = pSubMesh->cells.size();
	pSubMesh->v_sign.resize(nCellCount, std::vector<double>(pSubMesh->faces_per_cell));
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
	{
		// Cell Index
		const int nCellIndex = pSubMesh->cells[nCIndex].id;

		// Faces
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
		{
			// Face Index
			const int nFIndex = pSubMesh->cell2face[nCIndex][nCellFIndex];
			const int nFaceIndex = pSubMesh->faces[nFIndex].id;
			
			//Owner
			const int nOwner = (nFaceIndex == -1) ? nCellIndex : m_pMeshReader->getFaceOwner(nFaceIndex);

			// face contains the correct normal for its owner cell
			pSubMesh->v_sign[nCIndex][nCellFIndex] = (nOwner == nCellIndex) ? 1.0 : -1.0;
		}
	}

	// Get MPI Neighbours
	std::vector<bool> unregistered_mpi_neigh( m_MPI_Env.size(), true );
	const int nBoundaryCount = m_pMeshReader->getBoundaryCount();
	for (int nBoundaryIndex = 0; nBoundaryIndex < nBoundaryCount; nBoundaryIndex++)
	{
		// Boundary Rank
		const int nBoundaryRank = m_pMeshReader->getBoundaryProcessorRank(nBoundaryIndex);
		if (nBoundaryRank != -1 && unregistered_mpi_neigh[nBoundaryRank])
		{
			// Append to submesh flag
			bool bAppendToMesh = false;
			
			// Get Boundary Faces
			const int nFaceStartIndex = m_pMeshReader->getBoundaryFaceStart(nBoundaryIndex);
			const int nFaceEndIndex = m_pMeshReader->getBoundaryFaceEnd(nBoundaryIndex);
			for (int nFaceIndex = nFaceStartIndex; nFaceIndex <= nFaceEndIndex && !bAppendToMesh; nFaceIndex++)
			{
				// Face Cell Owner
				const int nCellOwner = m_pMeshReader->getFaceOwner(nFaceIndex);

				// Cell SM
				bAppendToMesh = (m_pMeshReader->getCellSubmeshIndex(nCellOwner) == nSubmeshIndex);
			}

			if (bAppendToMesh)
			{
				unregistered_mpi_neigh[nBoundaryRank] = false;
				pSubMesh->mpi_neighbors.push_back(nBoundaryRank);
			}			
		}
	}
}

// -------------------------------------------------------------------------- //
template <typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT>::calculateVNeighbour(const int nSubmeshIndex)
{
	// Submesh
	gmsh_mesh* pSubMesh = &m_SubmeshList[nSubmeshIndex];	

	// Assign Vector
	const int nFaceCount = pSubMesh->faces.size();
	pSubMesh->v_neigh.resize(nFaceCount, std::vector<double>(DIM_CNT));

	// CellCenter
	double dCellCenter[3];
	double dNeigbhourCenter[3];
	double dFaceCenter[3];
	double dCyclicFaceCenter[3];

	// MPI Communication
	const int nMPINeigbhourCount = pSubMesh->mpi_neighbors.size();

	// Boundary Count
	const int nBoundaryCount = m_pMeshReader->getBoundaryCount();
	
	// Map from global MPI_Rank to local MPI_Index
	std::map<int,int> rank_global2local;
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)
		rank_global2local.insert(pair<int,int>(pSubMesh->mpi_neighbors[nLocalIndex], nLocalIndex));

	// Handle neighbours on local process
	// Prepare all the faces with MPI neigbhour
	std::vector<std::vector<int>> nMPIFIndexList(nMPINeigbhourCount);
	std::vector<std::map<std::pair<int,int>,int>> faceTagMapping(nMPINeigbhourCount);

	// Tag Mapping
	std::vector<std::map<int,int>> tagMappingList(nMPINeigbhourCount);

	for (int nFIndex = 0; nFIndex < nFaceCount; nFIndex++)
	{
		// Face Index
		const int nFaceIndex = pSubMesh->faces[nFIndex].id;

		// Owner Cell Index
		const int nOwnerCellIndex = m_pMeshReader->getFaceOwner(nFaceIndex);

		// Neighbour Cell Index
		const int nNeigbhourCellIndex = m_pMeshReader->getFaceNeighbour(nFaceIndex);

		// Interior Face / Cyclic Face
		if (nNeigbhourCellIndex != -1)
		{
			// Cell Center
			m_pMeshReader->getCellCenter(nOwnerCellIndex, dCellCenter);
			m_pMeshReader->getCellCenter(nNeigbhourCellIndex, dNeigbhourCenter);

			for (int nD = 0; nD < DIM_CNT; nD++)
				pSubMesh->v_neigh[nFIndex][nD] = dNeigbhourCenter[nD] - dCellCenter[nD];

			// Cyclic Face
			const int nCyclicFaceIndex = m_pMeshReader->getCyclicFaceIndex(nFaceIndex);

			// Translate Cyclic BC - Assume faces must have same center
			if (nCyclicFaceIndex != nFaceIndex)
			{
				// Face Center
				m_pMeshReader->getFaceCenter(nFaceIndex, dFaceCenter);
				m_pMeshReader->getFaceCenter(nCyclicFaceIndex, dCyclicFaceCenter);
				for (int nD = 0; nD < DIM_CNT; nD++)				
					pSubMesh->v_neigh[nFIndex][nD] += dFaceCenter[nD] - dCyclicFaceCenter[nD];					
			}		
			continue;
		}

		// Boundary Face (must exist otherwise we are interior face)
		const int nBoundaryIndex = m_pMeshReader->getFaceBoundary(nFaceIndex);
		const int nBoundaryRank = m_pMeshReader->getBoundaryProcessorRank(nBoundaryIndex);		

		// Non-MPI boundary: Mirror ghost center
		if (nBoundaryRank == -1)
		{
			m_pMeshReader->getCellCenter(nOwnerCellIndex, dCellCenter);
			m_pMeshReader->getFaceCenter(nFaceIndex, dNeigbhourCenter);
			for (int nD = 0; nD < DIM_CNT; nD++)
				pSubMesh->v_neigh[nFIndex][nD] = 2 * (dNeigbhourCenter[nD] - dCellCenter[nD]);
			continue;
		}
		else
		{
			// Boundary Tag (use reverse since thats what the other processor will send)
			const int nTag = m_pMeshReader->getBoundaryTag(nBoundaryIndex);

			// Local MPI Index
			const int nLocalIndex = rank_global2local[nBoundaryRank];

			// MPI Faces List
			nMPIFIndexList[nLocalIndex].push_back(nFIndex);

			// FaceId
			const int nFaceId = m_pMeshReader->getFaceId(nFaceIndex);

			// Add Mapping
			faceTagMapping[nLocalIndex].insert(pair<pair<int,int>,int>(std::make_pair(nFaceId,-nTag), nFIndex));
		}
	}

	// No MPI Communication
	if (nMPINeigbhourCount == 0)
		return;
		
	// Get number of data required to send to each MPI neigbhour
	std::vector<int> num_mpi_neighs(nMPINeigbhourCount, 0);
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)
		num_mpi_neighs[nLocalIndex] = nMPIFIndexList[nLocalIndex].size();
		
		
	// Data Structure to send / recieve
	// Note: The Id we will send here will be the faceId and tag associated with this boundary!
	// That way, the other side can also get the matching data we have send to it
	struct t_cell_center
	{
		int id[3];
		double x[DIM_CNT];
	};

	// Define DataStructure in MPI	
	const int center_nfields = 2;
	MPI_Aint  center_disps[center_nfields];
	int       center_blocklens[] = { 3, DIM_CNT };
	MPI_Datatype center_types[] = { MPI_INT, MPI_DOUBLE };

	// Displacement in structure
	struct t_cell_center tmp_center;
	center_disps[0] = (char*) &tmp_center.id - (char*) &tmp_center;
	center_disps[1] = (char*)  tmp_center.x  - (char*) &tmp_center;

	MPI_Datatype type_center;
	MPI_Type_create_struct (center_nfields, center_blocklens, center_disps, center_types, &type_center);
	MPI_Type_commit(&type_center);
	
	// Assigned send/recieve vectors to/from each mpi rank
	std::vector<std::vector<struct t_cell_center>> send_cells_coord(nMPINeigbhourCount);
	std::vector<std::vector<struct t_cell_center>> recv_cells_coord(nMPINeigbhourCount);
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)
	{
		send_cells_coord[nLocalIndex].resize(num_mpi_neighs[nLocalIndex]);
		recv_cells_coord[nLocalIndex].resize(num_mpi_neighs[nLocalIndex]);
	}

	// Fill Send arrays
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)
	{
		// Face Index List
		const std::vector<int>& nFIndexList = nMPIFIndexList[nLocalIndex];
		for (int nIndex = 0; nIndex < nFIndexList.size(); nIndex++)
		{
			// FIndex
			const int nFIndex = nFIndexList[nIndex];

			// Face Index
			const int nFaceIndex = pSubMesh->faces[nFIndex].id;

			// FaceId
			const int nFaceId = m_pMeshReader->getFaceId(nFaceIndex);

			// Boundary Tag
			const int nBoundaryIndex = m_pMeshReader->getFaceBoundary(nFaceIndex);
			const int nTag = m_pMeshReader->getBoundaryTag(nBoundaryIndex);

			// Fill Send Data
			send_cells_coord[nLocalIndex][nIndex].id[0] = nFaceId;
			send_cells_coord[nLocalIndex][nIndex].id[1] = nTag;

			// Store Owner Cell Submesh Index
			// Note: This will be used to update our cell neighbour id
			const int nCellOwner = m_pMeshReader->getFaceOwner(nFaceIndex);
			send_cells_coord[nLocalIndex][nIndex].id[2] = m_pMeshReader->getSubmeshCellIndex(nSubmeshIndex, nCellOwner);

			// Store Displacement From FaceCenter
			m_pMeshReader->getCellCenter(nCellOwner, dCellCenter);
			m_pMeshReader->getFaceCenter(nFaceIndex, dFaceCenter);
			for (int nD = 0; nD < DIM_CNT; nD++)
				send_cells_coord[nLocalIndex][nIndex].x[nD] = dCellCenter[nD] - dFaceCenter[nD];
		}
	}

	// Send cell centers to MPI neighbors
	std::vector<MPI_Request> comm_reqs(2 * nMPINeigbhourCount);	
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)
	{
		// Neigbhour Rank
		const int nNeighbourRank = pSubMesh->mpi_neighbors[nLocalIndex];

		// Send / Recv Size (Same)
		const int nSendSize = send_cells_coord[nLocalIndex].size();
		const int nRecvSize = recv_cells_coord[nLocalIndex].size();

		// Perform Non-Blocking Send/Recieve
		MPI_Isend( &send_cells_coord[nLocalIndex][0], nSendSize, type_center, nNeighbourRank, 0, MPI_COMM_WORLD, &comm_reqs[nLocalIndex] );
		MPI_Irecv( &recv_cells_coord[nLocalIndex][0], nRecvSize, type_center, nNeighbourRank, 0, MPI_COMM_WORLD, &comm_reqs[nLocalIndex + nMPINeigbhourCount]);
	}

	// Wait for communication to end
	MPI_Waitall(comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE);

	// Extract cell centers from received data
	for (int nLocalIndex = 0; nLocalIndex < nMPINeigbhourCount; nLocalIndex++)	
		for (int nRecvIndex = 0; nRecvIndex < num_mpi_neighs[nLocalIndex]; nRecvIndex++)
		{
			// Face Id
			const int nFaceId = recv_cells_coord[nLocalIndex][nRecvIndex].id[0];

			// Tag
			const int nTag = recv_cells_coord[nLocalIndex][nRecvIndex].id[1];

			// FIndex
			const int nFIndex = faceTagMapping[nLocalIndex][std::make_pair(nFaceId, nTag)];			
			
			// Face Index
			const int nFaceIndex = pSubMesh->faces[nFIndex].id;

			// Cell Owner
			const int nCellOwner = m_pMeshReader->getFaceOwner(nFaceIndex);
			const int nCIndex = m_pMeshReader->getSubmeshCellIndex(nSubmeshIndex, nCellOwner);

			// Cell Center
			m_pMeshReader->getCellCenter(nCellOwner, dCellCenter);

			// Face Center
			m_pMeshReader->getFaceCenter(nFaceIndex, dFaceCenter);

			// Store Neigbhor Direction
			for (int nD = 0; nD < DIM_CNT; nD++)
				pSubMesh->v_neigh[nFIndex][nD] = recv_cells_coord[nLocalIndex][nRecvIndex].x[nD] + dFaceCenter[nD] - dCellCenter[nD];

			// Update Corresponding inner face neighbour
			for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
				if (pSubMesh->cell2face[nCIndex][nCellFIndex] == nFIndex)
				{
					pSubMesh->cell2neigh[nCIndex][nCellFIndex].id = recv_cells_coord[nLocalIndex][nRecvIndex].id[2];
					break;
				}
		}
	MPI_Barrier(MPI_COMM_WORLD);
}




