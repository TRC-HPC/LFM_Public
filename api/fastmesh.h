/**
 * Copyright (C) Huawei Technologies Co., Ltd. 2020. ALL RIGHTS RESERVED.
 *
 * See file LICENSE for terms.
 *
 *
 * (lib)FastMesh
 * =============
 *
 * This is a prototype for a new memory layout for representing a mesh, intended
 * to facilitate faster iterative solvers. This is achieved by enhanced memory
 * proximity and a layout which is friendly to both SIMD instructions and
 * network offloading APIs (i.e. RDMA).
 */

#ifndef NEIGHBOR_H_
#define NEIGHBOR_H_

#include <math.h>
#include <vector>
#include <memory>
#include <ostream>
#include <iostream>
#include <functional>

#include "elements.h"
#include "mpi_env.h"
#include "cfdv0_solver.h"
#include "api/meshReader.h"
#include "api/inputReader.h"
#include "api/runTimeManager.h"

namespace fastmesh {

//***************************************************************************************************
// List of solvers
//***************************************************************************************************
typedef enum fastmesh_solvers {
    SOLVER_CAAFOAM = 0, SOLVER_M2 = 1, SOLVER_M2PDISS = 2, 
	SOLVER_COUNT,
} fastmesh_solver_t;

class IMesh
{
public:
	// Initialize
	virtual void initialize() = 0;

	// Initialize Solver
	virtual int initializeSolver() = 0;

	// Solver
	virtual void solve() = 0;

	// Finalize
	virtual void finalize() = 0;

	// Input Reader
	virtual void setInputReader(CInputReader *pInputReader) = 0;

	// RunTime Manager
	virtual void setRunTimeManager(IRunTimeManager *pRunTimeManager) = 0;

};

//***************************************************************************************************
// The Mesh class focuses on data exchange between Submesh instances (while the Submesh focuses on
// solvers). This means the mesh is responsible for moving the contents of boundary cells between
// Submeshes - after every iteration. This exchange can be either within the same host or over the
// network.
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT>
class Mesh : public IMesh {

public:
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    /**< Create an empty mesh */
    Mesh() : m_pMeshReader(nullptr) {};
	virtual ~Mesh() {if (m_pMeshReader) delete m_pMeshReader;}

    /**< Clone an existing mesh */
    Mesh( Mesh<PRECISION, DIM_CNT> *original) {};

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Geometrical and Solver objects
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::string m_sFileBase;
	IMeshReader* m_pMeshReader;
    std::vector<gmsh_mesh> m_SubmeshList;
    std::vector<ISolver *> m_pSolverList;
	std::vector<ISolver *> m_pInteriorSolverList;

	// Input Reader
	CInputReader* m_pInput;

	// RunTime Manager
	IRunTimeManager *m_pRunTime;
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // MPI environment
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MPI_env m_MPI_Env;
    std::vector<int> m_nMPINeighbourList;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Functions
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Input Reader
	virtual void setInputReader(CInputReader *pInputReader) {m_pInput = pInputReader;}

	// Run Time Manager
	virtual void setRunTimeManager(IRunTimeManager *pRunTimeManager) {m_pRunTime = pRunTimeManager;}

    virtual void initialize();
	virtual int initializeSolver();

	// New Output Methods
	void updateAverages(const int nTimeStep);
	void updateSolutionResidual();
	void updateSolution();
	void updateSolutionBlendFactor();
	
	void output_Forces  ( bool firstRun, size_t numberBodies);

	virtual void solve();

	// Wrap-up simulation, including MPI environmnet
	virtual void finalize()
	{
		for (ISolver* pSolver : m_pSolverList)
		{
			pSolver->deallocate();
			delete pSolver;
		}

		m_MPI_Env.finalize();
	}

	void readMeshFiles();
	void readMeshNodes(const int nSubmeshIndex);
	void readMeshFaces(const int nSubmeshIndex);
	void readMeshCells(const int nSubmeshIndex);
	void setParams(const int nSubmeshIndex);
	void calculateVNeighbour(const int nSubmeshIndex);
};

#ifdef SINGLE_PRECISION
template class Mesh<float ,2>;
template class Mesh<float ,3>;
#endif

template class Mesh<double ,2>;
template class Mesh<double ,3>;

} // namespace fastmesh

#endif // NEIGHBOR_H_








































