#include <cstddef>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unistd.h>
#include <signal.h>
#include <execinfo.h>
#include <string.h>
#include <fenv.h>
#include "api/inputReader.h"
#include "api/fastmesh.h"
#include "api/runTimeManagerOF.h"
#include "api/polyMeshReaderOF.h"
#include "omp.h"

//***************************************************************************************************
int main( int argc, char *argv[] )
{
	// Close Profiling
	EXTRAE_SHUTDOWN

	// Initialize MPI
    MPI_Init(&argc, &argv);

	// Get rank & process count
	int nRank, nRankCount;
	MPI_Comm_rank(MPI_COMM_WORLD, &nRank);
	MPI_Comm_size( MPI_COMM_WORLD, &nRankCount );

	// Master
	const bool bIsMaster = (nRank == 0);

	// Invalid process count
	if(nRankCount > MAX_MPI_RANKS)
	{
		if (bIsMaster)
			std::cerr << "Exceeded maximum allowed number of ranks: " << nRankCount << " > " << MAX_MPI_RANKS << std::endl;

		MPI_Finalize();
		return 0;
	}

	// Read Input
	CInputReader inputReader;
	if (!inputReader.initialize(argc, argv))
	{
		MPI_Finalize();
		return 0;
	}

	// Invalid Solver
	if (inputReader.m_nSolver < 0 || inputReader.m_nSolver >= fastmesh::SOLVER_COUNT)
	{
		if (bIsMaster)
			std::cerr << "unknown solver type: " << inputReader.m_nSolver << std::endl;
		MPI_Finalize();
		return 0;
	}

	// Invalid Dimension
	if (inputReader.m_nDimension != 2 && inputReader.m_nDimension != 3)
	{
		if (bIsMaster)
			std::cerr << "invalid dimension: " << inputReader.m_nDimension << std::endl;
		MPI_Finalize();
		return 0;
	}

	// Invalid Comm Type
	if (inputReader.m_nCommType < 0 || inputReader.m_nCommType >= MPI_EXCHANGE_COUNT)
	{
		if (bIsMaster)
			std::cerr << "invalid communication type: " << inputReader.m_nCommType << std::endl;
		MPI_Finalize();
		return 0;
	}
	
	// Invalid Halo Comm Type
	if (inputReader.m_nHaloCommType < 0 || inputReader.m_nHaloCommType >= MPI_TYPE_COUNT)
	{
		if (bIsMaster)
			std::cerr << "invalid halo communication type: " << inputReader.m_nHaloCommType << std::endl;
		MPI_Finalize();
		return 0;
	}

	#ifndef SINGLE_PRECISION
	if (!inputReader.m_bIsDoublePrecision)
	{
		if (bIsMaster)
			std::cerr << "Single precision is not supported unless 'SINGLE_PRECISION' is defined" << std::endl;
		MPI_Finalize();
		return 0;
	}
	#endif

	// Use required mesh
	fastmesh::IMesh* pMesh = nullptr;
	if (inputReader.m_bIsDoublePrecision)
	{
		if (inputReader.m_nDimension == 2) 
			pMesh = new fastmesh::Mesh<double, 2>;
		else 
			pMesh = new fastmesh::Mesh<double, 3>;
	}
	else
	{
#ifdef SINGLE_PRECISION
		if (inputReader.m_nDimension == 2) 
			pMesh = new fastmesh::Mesh<float, 2>;
		else 
			pMesh = new fastmesh::Mesh<float, 3>;
#endif
	}

	// RunTime Manager
	const int nParallelRank = inputReader.m_bParallel ? nRank : -1;
	CRunTimeManagerOF runTime(nParallelRank);
	pMesh->setRunTimeManager(&runTime);

	// Store Input Reader Pointer
	pMesh->setInputReader(&inputReader);

	// Read PolyMesh
	pMesh->initialize();

	// Initialize mesh solver	
	if(pMesh->initializeSolver() != 0)
	{
		MPI_Finalize();
		return 0;
	}
	
	// Solve CFD
	pMesh->solve();

	// Finalize
	pMesh->finalize();

	// Delete
	delete pMesh;

	return 0;
}