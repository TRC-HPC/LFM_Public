#include <iostream>
#include <iomanip>
#include <fstream>
#include <ctime>
#include <sstream>
#include <cstring>
#include <array>
#include <limits>
#include <map>

#include <sys/stat.h>

#include "api/fastmesh.h"
#include "api/mpi_env.h"
#include "api/dictReaderOF.h"

using namespace std;
using namespace fastmesh;

std::string upperCase(std::string input) {
  for (std::string::iterator it = input.begin(); it != input.end(); ++ it)
    *it = toupper(*it);
  return input;
}

//***************************************************************************************************
// Initialize solver
// 	- for now, only jacobi solver
// 	- next: extend for multiple solvers w/ function pointers
//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT >
int Mesh<PRECISION, DIM_CNT> :: initializeSolver()
{	
	// Submesh Count
	const int nSubmeshCount = m_SubmeshList.size();
	
	// Pointers to:
	//	cells_ptr: Cells variables (vector for each submesh)
	//  ghost_ptr: Ghost MPI boundary cells variables (vector for each MPI neighbour)
	//  ghost_bnd_ptr: Ghost local boundary cells variables (vector for each boundary type)
	vector<vector<void*>> cells_ptr(nSubmeshCount);
	vector<vector<void*>> ghost_ptr;
	vector<vector<void*>> ghost_bnd_ptr;

	MPI_Barrier(MPI_COMM_WORLD);

	// Check Each Submesh	
	for (int nSubmeshIndex = 0; nSubmeshIndex < nSubmeshCount; nSubmeshIndex++)
	{
		// SubMesh
		gmsh_mesh* pSubMesh = &m_SubmeshList[nSubmeshIndex];

		// Max Faces Per cell
		const int nMaxFacesPerCell = pSubMesh->faces_per_cell;

		switch (nMaxFacesPerCell)
		{
			case 3:
				m_pSolverList.push_back(new CFDv0_solver<PRECISION, DIM_CNT, 3>);
				break;
			case 4:
				m_pSolverList.push_back(new CFDv0_solver<PRECISION, DIM_CNT, 4>);
				break;
			case 5:
				m_pSolverList.push_back(new CFDv0_solver<PRECISION, 3, 5>);
				break;
			case 6:
				m_pSolverList.push_back(new CFDv0_solver<PRECISION, 3, 6>);
				break;
			default:
				std::cerr << "Invalid face cell count: " << nMaxFacesPerCell << std::endl;
				MPI_Abort( MPI_COMM_WORLD, 100 );			
				break;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Interior Solvers	
	for (int nSubmeshIndex = INDEX_BND_SUBMESH + 1; nSubmeshIndex < nSubmeshCount; nSubmeshIndex++)
		m_pInteriorSolverList.push_back(m_pSolverList[nSubmeshIndex]);

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Allocate cells and ghost cells
	for (int nSubmeshIndex = 0; nSubmeshIndex < nSubmeshCount; nSubmeshIndex++)
	{
		// SubMesh
		gmsh_mesh* pSubMesh = &m_SubmeshList[nSubmeshIndex];

		// Solver
		ISolver* pSolver = m_pSolverList[nSubmeshIndex];

		pSolver->m_pMeshReader = m_pMeshReader;
		pSolver->m_nSubmeshIndex  = nSubmeshIndex;
		pSolver->m_pInput = m_pInput;
		pSolver->allocate(pSubMesh, cells_ptr[nSubmeshIndex]);
		pSolver->allocate_ghost_cells(m_MPI_Env, pSubMesh, ghost_ptr, ghost_bnd_ptr);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	
	// Set pointers to neighbouring cells
	for (int nSubmeshIndex = 0; nSubmeshIndex < nSubmeshCount; nSubmeshIndex++)
	{
		// SubMesh
		gmsh_mesh* pSubMesh = &m_SubmeshList[nSubmeshIndex];

		// Solver
		ISolver* pSolver = m_pSolverList[nSubmeshIndex];

		pSolver->assign_pointers(pSubMesh, m_MPI_Env, cells_ptr, ghost_ptr, ghost_bnd_ptr);
		pSolver->initialize(m_MPI_Env, pSubMesh, m_pInput->m_nSolver);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Reorder Faces
	for (ISolver* pSolver : m_pSolverList)
		pSolver->reorder_faces();

	MPI_Barrier(MPI_COMM_WORLD);

	return 0;
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: updateAverages(const int nTimeStep)
{
	// Ask all submeshes to update averages scalar field
	if(nTimeStep > 0){
		for (ISolver* pSolver : m_pSolverList)
			pSolver->updateAverageField(nTimeStep);
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: updateSolutionResidual()
{
	// Ask all submeshes to update residual scalar field
	for (ISolver* pSolver : m_pSolverList)
		pSolver->updateSolutionResidual();
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: updateSolution()
{
	// Ask all submeshes to update residual scalar field
	for (ISolver* pSolver : m_pSolverList)
		pSolver->updateSolutionPrimitives();	
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: updateSolutionBlendFactor()
{
	// Ask all submeshes to update residual scalar field
	for (ISolver* pSolver : m_pSolverList)
		pSolver->updateSolutionBlendFactor();
}


//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
void Mesh<PRECISION, DIM_CNT> :: solve()
{
	if( m_MPI_Env.is_master() )
	{
		cout << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "CFD solver: caafoam" << endl;
		cout << "    - Flow type    : Navier-Stokes equation" << endl;
		cout << "    - MPI framework: ";
		switch(m_MPI_Env.m_nCommType){
			case MPI_EXCHANGE_FULL_BND:
				cout << "full bnd Hpath exchange with all neighbors" << endl;
				break;
			case MPI_EXCHANGE_PACKED:
				cout << "packed exchange to each neighbor" << endl;
				break;
			case MPI_EXCHANGE_SPLIT:
				cout << "split solvars/viscous communication steps" << endl;
				break;
		}
		cout << "    - MPI protocol : ";
		switch( m_MPI_Env.get_halo_comm_type() ){
			case MPI_ONESIDED_NONB:
				cout << "one-sided, non-blocking" << endl;
				break;
			case MPI_ONESIDED_BLCK:
				cout << "one-sided, blocking" << endl;
				break;
			case MPI_TWOSIDED_BLCK:
				cout << "two-sided, blocking" << endl;
				break;
			case MPI_TWOSIDED_NONB:
				cout << "two-sided, non-blocking" << endl;
				break;
			case MPI_TWOSIDED_PERS:
				cout << "two-sided, persistent" << endl;
				break;
			case MPI_NEIGHCOLL_BLCK:
				cout << "neighborhood collective, blocking" << endl;
				break;
			case MPI_NEIGHCOLL_NONB:
				cout << "neighborhood collective, non-blocking" << endl;
				break;
			case MPI_NEIGHCOLL_PERS:
				cout << "neighborhood collective, persistent" << endl;
				break;
			default:
				cout << "None: serial simulation" << endl;
				break;
		}
		cout << "    - MPI cores    :" << setw(10) << m_MPI_Env.size() << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		if(m_pInput->m_nSolver == fastmesh::SOLVER_CAAFOAM){
			cout << "Running with M1 discretization\n";
		}
		else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2){
			cout << "Running with M2 discretization\n";
		}
		else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2PDISS){
			cout << "Running with M2 AUSM discretization\n";
		}
		if(m_pInput->m_bisLaminar){
			cout << "Laminar simulation\n";
		}else{
			cout << "Turbulent simulation with Smagorinsky model\n";
		}
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl << endl;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialization phase
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	// Residual
	PRECISION *RES = NULL;
	
	PRECISION dt = m_pRunTime->getDeltaTime();
	PRECISION t_end = m_pRunTime->getEndTime();
	PRECISION time = m_pRunTime->getCurrentTime();
	PRECISION cflMax = m_pInput->m_dCFLMax;

	PRECISION tStartAverage = m_pInput->m_tStartAverage;
	int indexStartAverage = (tStartAverage-time)/dt;

	MPI_Barrier(MPI_COMM_WORLD);
	
	if( m_MPI_Env.is_master() ){
		cout << "Statistics will start being computed at t = " << tStartAverage << ", at iteration " << indexStartAverage << endl;
	}
	int time_step     = 0;

	vector<probe> x_probes;
	vector<probe> x_samples;

	// Create Average / RMS Fields
	if (m_pInput->m_bHaveAverage)
	{
		m_pMeshReader->createScalarField("pAvg");
		m_pMeshReader->createScalarField("pRMS");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Create Primitive Fields
	{
		m_pMeshReader->createScalarField("rho");		
		m_pMeshReader->createScalarField("E");
		m_pMeshReader->createScalarField("p");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Create Residual Fields
	if (m_pInput->m_bHaveResidual && m_pInput->m_bSaveResiduals)
	{
		string uvwFields[3] = {"uRES", "vRES", "wRES"};
		m_pMeshReader->createScalarField("rhoRES");
		m_pMeshReader->createScalarField("eRES");
		for (int nD = 0; nD < DIM_CNT; nD++)
			m_pMeshReader->createScalarField(uvwFields[nD]);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Create Blend Output (sponge_layer)
	if (m_pInput->m_bSaveBlendFactor)
		m_pMeshReader->createScalarField("sigma0");

	// Save Rank
	if (m_pInput->m_bSaveRank)
	{
		const PRECISION rank = m_MPI_Env.rank();
		const int nRankField = m_pMeshReader->createScalarField("rank");
		for (int nCellIndex = m_pMeshReader->getCellCount() - 1; nCellIndex >= 0; nCellIndex--)
			m_pMeshReader->updateScalarField(nRankField, rank, nCellIndex);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Save Cell Order
	if (m_pInput->m_bSaveCellOrder)
	{
		const int nCellOrderField = m_pMeshReader->createScalarField("cellOrder");
		for (int nCellIndex = m_pMeshReader->getCellCount() - 1; nCellIndex >= 0; nCellIndex--)
			m_pMeshReader->updateScalarField(nCellOrderField, static_cast<float>(nCellIndex), nCellIndex);
	}
		
	MPI_Barrier(MPI_COMM_WORLD);

	// Allocate arrays for averaging, RMS and forces	
	for (int nSubmeshIndex = 0; nSubmeshIndex < m_pSolverList.size(); nSubmeshIndex++)
		m_pSolverList[nSubmeshIndex]->prePostProc( m_MPI_Env, m_pInput->m_bHaveProbes, m_pInput->m_bHaveSampling, 
		m_pInput->m_bHaveAverage, m_pInput->m_bHaveForces, x_probes, x_samples, nSubmeshIndex);

	// Allocate residual array
	if(m_pInput->m_bHaveResidual)
	{
		RES = new PRECISION[DIM_CNT+2];
		for (int i = 0; i < DIM_CNT + 2; i++) RES[i] = 0;
	}
	
	struct stat info;
	if( stat( "./output/", &info ) != 0 )
	{
		mkdir("./output", 0777);
		mkdir("./output/Times", 0777);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Wall BC Count
	ISolver* pBCSolver = m_pSolverList[INDEX_BND_SUBMESH];
	const int numberBodies = pBCSolver->m_nWallBCList.size();
	if ( m_pInput->m_bHaveForces ) {
		if( m_MPI_Env.is_master() ) 
		{
			std::stringstream tmp_filename;
			tmp_filename << "./output/postProcess_Forces_1";

			char filename[tmp_filename.str().size()+1];
			strcpy(filename, tmp_filename.str().c_str());

			// Write forces file header (if not found)
			if( stat( filename, &info ) != 0 )
				output_Forces(true, numberBodies);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// CPU timing
	uint64_t t0, t1, t2, t3, t5, t6, t10, t11, t12, t13;
	double *t_submesh = new double [m_pSolverList.size()];
	double t_init     = 0.;
	double t_solver   = 0.;
	double t_step     = 0.;
	double t_compute_step  = 0.;
	double t_postPro  = 0.;
	double t_output   = 0.;
	double t_visc     = 0.;
	double t_ghost    = 0.;
	double t_mpi      = 0.;
	double t_mpi_step = 0.;
	double t_grads    = 0.;
	double t_wait0    = 0.;
	double t_wait1    = 0.;
	double t_wait2    = 0.;
	double t_wait3    = 0.;
	double t_int0     = 0.;
	double t_int1     = 0.;
	double t_bnd      = 0.;
	double t_int      = 0.;

	for(int nSubmeshIndex=0; nSubmeshIndex < m_pSolverList.size(); nSubmeshIndex++ )
		t_submesh[nSubmeshIndex] = 0.;

	MPI_Datatype mpi_precision;
	if( sizeof(PRECISION) == sizeof(double) )
		mpi_precision = MPI_DOUBLE;
	else
		mpi_precision = MPI_FLOAT;


	MPI_Barrier(MPI_COMM_WORLD);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// CFD solver algorithm. For each time-step:
	// 		- prepapre solver for time-step
	// 		- solve for boundary submesh
	// 		- start non-blocking MPI exchange of ghost cells
	// 		- solve for interior submeshes
	// 		- wait for MPI comm to be done
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Boundary Solver
	MPI_Barrier( MPI_COMM_WORLD );

	// Actual MPI communication
	pBCSolver->mpi_communication(m_MPI_Env, 0);

	// Set Boundary Condition
	pBCSolver->set_boundary_conditions();

	// Wait for ghost cells sent through MPI to receive
	pBCSolver->mpi_wait(m_MPI_Env, 0);

	// Send viscous data (if sent is split to two parts)
	pBCSolver->mpi_communication(m_MPI_Env, 1);
	pBCSolver->mpi_wait(m_MPI_Env, 1);

	// Calculate viscous 
	if(m_pInput->m_bisLaminar){
		pBCSolver->calc_VIS(m_MPI_Env);
	} else {
		pBCSolver->calc_VIS_Smagorinsky(m_MPI_Env);
	}
	if ( m_pInput->m_bHaveForces )
	 	pBCSolver->postProcForces( m_MPI_Env, m_pRunTime->getCurrentTime(), mpi_precision );

	// Update blend factor (sponge_sigma)
	if (m_pInput->m_bSaveBlendFactor)
		updateSolutionBlendFactor();

	MPI_Barrier( MPI_COMM_WORLD );

	// Enable MPI Profiling
 	MPI_Pcontrol(1);

	// ----------- //
	// Timing Info //
	// ----------- //
	double dCompleteTime = 0;

	// Prepare Time
	double dPrepareTime = 0;

	// RK Loop Time
	double dRKLoopTime = 0;
	double dRkPrepareTime = 0;
	double dRkMPIWaitTime = 0;
	double dRkSetBCTime = 0;
	double dRkMinmodTime = 0;
	double dRkVisTime = 0;
	double dRkCommTime = 0;
	double dRkStepTime = 0;

	// Outside Rk Loop
	double dOutsideMPIWaitTime = 0;
	double dAdjustCflOrDTTime = 0;
	double dAdvanceTime = 0;
	double dOutputTime = 0;

	struct timespec nStartTime, nEndTime, nStartRkLoopTime, nEndRkLoopTime, nTs, nTe;

	// Start Measure Time	
	GETTIME(nStartTime);

	t0 = clock();

	EXTRAE_RESTART
	EXTRAE_EVENT(1110000, 1)

	int nLoopIndex = 1;
	while (m_pRunTime->isRunning())
	{	
		// Start Loop Event
		EXTRAE_EVENT(1110001, nLoopIndex++);

		// Timing
		GETTIME(nTs);

		// Get TimeStep (if changed)
		dt = m_pRunTime->getDeltaTime();

		// Set DQ and RES to 0 for all cell_cfd
		EXTRAE_EVENT(1110002, 1)
		for (ISolver* pSolver : m_pSolverList)
			pSolver->prepare_for_timestep();
		EXTRAE_EVENT(1110002, 0)

		// Timing
		GETTIME(nTe);
		dPrepareTime += DIFFTIME(nTe, nTs);
		
		// For each RK-step...
		t_mpi_step = 0.0;
		t5 = clock();
		GETTIME(nStartRkLoopTime);
		EXTRAE_EVENT(1110003, 1)
		for( int rk_step=0; rk_step < m_pInput->m_nRkStepCount; rk_step++ ){
			// Prepare for next iteration
			GETTIME(nTs);
			t2 = clock();
			EXTRAE_EVENT(1110004, 1)
			for (ISolver* pSolver : m_pSolverList)
				pSolver->prepare_for_RKstep( rk_step, m_pInput->m_nSolver );
			EXTRAE_EVENT(1110004, 0)			
			GETTIME(nTe);
			dRkPrepareTime += DIFFTIME(nTe, nTs);
			t3 = clock(); t_init += (double) (t3-t2) / CLOCKS_PER_SEC;

			// Wait for MPI comm of last iteration to finish
			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110005, 1)
			pBCSolver->mpi_wait(m_MPI_Env, 0);
			EXTRAE_EVENT(1110005, 0)
			GETTIME(nTe);
			dRkMPIWaitTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_wait0 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_mpi_step += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Boundary submesh: gradients & viscous first step
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// physical boundary conditions must be set before calc_VIS and calc_gradients are called
			GETTIME(nTs);
			EXTRAE_EVENT(1110006, 1)
			pBCSolver->set_boundary_conditions( );
			EXTRAE_EVENT(1110006, 0)
			GETTIME(nTe);
			dRkSetBCTime += DIFFTIME(nTe, nTs);

			// Note: mpi_env is not used in calc_gradients and calc_VIS
			GETTIME(nTs);
			EXTRAE_EVENT(1110007, 1)
			if(m_pInput->m_bMinmodInterpolation)
			{
				t10 = clock();
				
				if (m_pInput->m_nSolver == fastmesh::SOLVER_M2PDISS) 
					pBCSolver->calc_gradients_M2AUSM(m_MPI_Env);
				else
					pBCSolver->calc_gradients( m_MPI_Env );

				t11 = clock(); t_grads += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_bnd += (double) (t11-t10) / CLOCKS_PER_SEC;
			}
			EXTRAE_EVENT(1110007, 0)
			GETTIME(nTe);
			dRkMinmodTime += DIFFTIME(nTe, nTs);
			t10 = clock();

			GETTIME(nTs);
			EXTRAE_EVENT(1110008, 1)
			if(m_pInput->m_bisLaminar){
				pBCSolver->calc_VIS(m_MPI_Env);
			} else {
				pBCSolver->calc_VIS_Smagorinsky(m_MPI_Env);
			}
			EXTRAE_EVENT(1110008, 0)
			GETTIME(nTe);
			dRkVisTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_visc += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_bnd += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110009, 2)
			pBCSolver->mpi_communication(m_MPI_Env, 1);
			EXTRAE_EVENT(1110009, 0)
			GETTIME(nTe);
			dRkCommTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_ghost += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_mpi_step += (double)(t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Interior submesh: gradients & viscous first step
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			GETTIME(nTs);
			EXTRAE_EVENT(1110007, 1)
			if( m_pInput->m_bMinmodInterpolation ){
				t10 = clock();
				for (ISolver* pInteriorSolver : m_pInteriorSolverList){					
					if (m_pInput->m_nSolver == fastmesh::SOLVER_M2PDISS) 
						pInteriorSolver->calc_gradients_M2AUSM(m_MPI_Env);
					else
						pInteriorSolver->calc_gradients( m_MPI_Env);
				}
				
				t11 = clock(); t_grads += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_int0 += (double) (t11-t10) / CLOCKS_PER_SEC;
				t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;
			}
			EXTRAE_EVENT(1110007, 0)
			GETTIME(nTe);
			dRkMinmodTime += DIFFTIME(nTe, nTs);

			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110008, 1)
			for (ISolver* pInteriorSolver : m_pInteriorSolverList){
				if(m_pInput->m_bisLaminar){
					pInteriorSolver->calc_VIS(m_MPI_Env);
				} else {
					pInteriorSolver->calc_VIS_Smagorinsky(m_MPI_Env);
				}
			}
			EXTRAE_EVENT(1110008, 0)
			GETTIME(nTe);
			dRkVisTime += DIFFTIME(nTe, nTs);
			
			t11 = clock(); t_visc += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int0 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110005, 2)
			pBCSolver->mpi_wait(m_MPI_Env, 1);
			EXTRAE_EVENT(1110005, 0)
			GETTIME(nTe);
			dRkMPIWaitTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_wait1 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_mpi_step += (double) (t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Boundary submesh: advance time-step (convection + viscous, step 2)
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110010, 1)
			if(m_pInput->m_nSolver == fastmesh::SOLVER_CAAFOAM){
				pBCSolver->one_rk_step_M1( rk_step, dt, m_MPI_Env, RES );
			}
			else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2){
				pBCSolver->one_rk_step_M2( rk_step, dt, m_MPI_Env, RES);
			}
			else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2PDISS){
				pBCSolver->one_rk_step_M2AUSM( rk_step, dt, m_MPI_Env, RES);
			}
			EXTRAE_EVENT(1110010, 0)
			GETTIME(nTe);
			dRkStepTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_solver += (double) (t11-t10) / CLOCKS_PER_SEC;

			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110009, 1)
			pBCSolver->mpi_communication(m_MPI_Env, 0);
			EXTRAE_EVENT(111009, 0)
			GETTIME(nTe);
			dRkCommTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_ghost += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_mpi_step += (double)(t11-t10) / CLOCKS_PER_SEC;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Interior submesh: advance time-step (convection + viscous, step 2)
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			t10 = clock();
			GETTIME(nTs);
			EXTRAE_EVENT(1110010, 1)
			for (ISolver* pInteriorSolver : m_pInteriorSolverList)
				if(m_pInput->m_nSolver == fastmesh::SOLVER_CAAFOAM){
					pInteriorSolver->one_rk_step_M1( rk_step, dt, m_MPI_Env, RES );
				}
				else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2){
					pInteriorSolver->one_rk_step_M2( rk_step, dt, m_MPI_Env, RES);
				}
				else if(m_pInput->m_nSolver == fastmesh::SOLVER_M2PDISS){
					pInteriorSolver->one_rk_step_M2AUSM( rk_step, dt, m_MPI_Env, RES);
				}
			EXTRAE_EVENT(1110010, 0)
			GETTIME(nTe);
			dRkStepTime += DIFFTIME(nTe, nTs);
			t11 = clock(); t_solver += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int1 += (double) (t11-t10) / CLOCKS_PER_SEC;
			t_int  += (double) (t11-t10) / CLOCKS_PER_SEC;
		}

		// Finish Rk Event
		EXTRAE_EVENT(1110003, 0)
		GETTIME(nEndRkLoopTime);
		dRKLoopTime += DIFFTIME(nEndRkLoopTime, nStartRkLoopTime);


		// Wait for halo exchange after the last RK step
		t10 = clock();
		GETTIME(nTs);
		EXTRAE_EVENT(1110005, 1)
		pBCSolver->mpi_wait(m_MPI_Env, 0);
		EXTRAE_EVENT(1110005, 0)
		GETTIME(nTe);
		dOutsideMPIWaitTime += DIFFTIME(nTe, nTs);
		t11 = clock(); t_wait2 += (double) (t11-t10) / CLOCKS_PER_SEC;
		t_mpi_step += (double) (t11-t10) / CLOCKS_PER_SEC;

		t6 = clock(); t_step = (double) (t6-t5) / CLOCKS_PER_SEC;

		t_compute_step = t_step - t_mpi_step;

		// CFL
		GETTIME(nTs);
		EXTRAE_EVENT(1110011, 1)
		if (m_pRunTime->isAdjustDeltaTime()) 
		{
			PRECISION dt_local = numeric_limits<PRECISION>::max();
			for (ISolver* pSolver : m_pSolverList)
			{
				PRECISION tmp_dt = pSolver->compute_dt( cflMax );
				if( dt_local > tmp_dt )
					dt_local = tmp_dt;
			}

			MPI_Allreduce( &dt_local, &dt, 1, mpi_precision, MPI_MIN, MPI_COMM_WORLD );
			m_pRunTime->setDeltaTime(dt);
		} 
		else 
		{
			PRECISION cfl_local = 0.0;
			for (ISolver* pSolver : m_pSolverList){
				PRECISION tmp_cfl = pSolver->compute_cfl( dt );
				if( cfl_local < tmp_cfl )
					cfl_local = tmp_cfl;
			}
			cflMax = cfl_local;
		}
		EXTRAE_EVENT(1110011, 0)
		GETTIME(nTe);
		dAdjustCflOrDTTime += DIFFTIME(nTe, nTs);

		// Advance Time
		GETTIME(nTs);
		EXTRAE_EVENT(1110012, 1)
		m_pRunTime->advanceTime();
		GETTIME(nTe);		
		time += dt;
		time_step++;
		dAdvanceTime += DIFFTIME(nTe, nTs);		

		// Post-Process functions execution
		GETTIME(nTs);		
		if ( m_pInput->m_bHaveAverage ) {
			t12 = clock();
			for (ISolver* pSolver : m_pSolverList)
				pSolver->postProcAverage(time_step-indexStartAverage);
			
			t13 = clock(); t_postPro += (double) (t13-t12) / CLOCKS_PER_SEC;
		}
		
		if ( m_pInput->m_bHaveForces && time_step % m_pInput->m_nSaveForcesSteps == 0 ) {
			t12 = clock();
			pBCSolver->postProcForces( m_MPI_Env, time, mpi_precision );
			t13 = clock(); t_postPro += (double) (t13-t12) / CLOCKS_PER_SEC;
		}
		EXTRAE_EVENT(1110012, 0)

		// Sim info
		EXTRAE_EVENT(1110013, 1)
		if( time_step % m_pInput->m_nPrintInfoFrequency == 0 || time >= t_end ){
			vector<struct t_domain_info<PRECISION, DIM_CNT>> domain_info( m_pSolverList.size() );
			PRECISION RES_global[DIM_CNT+2];
			if(m_pInput->m_bHaveResidual)
				MPI_Reduce( RES, RES_global, DIM_CNT+2, mpi_precision, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD );

			double t_compute_step_global, t_mpi_step_global, t_step_global;
			MPI_Reduce( &t_compute_step, &t_compute_step_global, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_mpi_step, &t_mpi_step_global, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_step, &t_step_global, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD );

			double t_compute_step_min, t_mpi_step_min, t_step_min;
			MPI_Reduce( &t_compute_step, &t_compute_step_min, 1, MPI_DOUBLE, MPI_MIN, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_mpi_step, &t_mpi_step_min, 1, MPI_DOUBLE, MPI_MIN, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_step, &t_step_min, 1, MPI_DOUBLE, MPI_MIN, m_MPI_Env.get_master(), MPI_COMM_WORLD );

			double t_compute_step_max, t_mpi_step_max, t_step_max;
			MPI_Reduce( &t_compute_step, &t_compute_step_max, 1, MPI_DOUBLE, MPI_MAX, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_mpi_step, &t_mpi_step_max, 1, MPI_DOUBLE, MPI_MAX, m_MPI_Env.get_master(), MPI_COMM_WORLD );
			MPI_Reduce( &t_step, &t_step_max, 1, MPI_DOUBLE, MPI_MAX, m_MPI_Env.get_master(), MPI_COMM_WORLD );

			if( m_MPI_Env.is_master() ) {
				cout << setprecision(10) << defaultfloat <<
						setw(10) << "Time step: " << time_step <<
						setw(10) << "Time: " << time;
						if (m_pRunTime->isAdjustDeltaTime())
							cout << setw(10) << "dt: " << dt;
						else
							cout << setw(10) << "CFL: " << cflMax;
						
						if(m_pInput->m_bHaveResidual) {
							cout << setw(15) << "Residuals: ";
							for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
								cout << setw(25) << setprecision(15) << scientific << sqrt(RES_global[idim]);
						}
						cout << endl << endl;
						cout << setprecision(4) << defaultfloat << fixed
							 << setw(15) << "Compute time: ("
							 << setw(8) << t_compute_step_min * m_pInput->m_nPrintInfoFrequency << ", "
							 << setw(8) << t_compute_step_global * m_pInput->m_nPrintInfoFrequency / m_MPI_Env.size() << ", "
							 << setw(8) << t_compute_step_max * m_pInput->m_nPrintInfoFrequency << ")"
							 << setw(15) << "MPI time: ("
							 << setw(8) << t_mpi_step_min * m_pInput->m_nPrintInfoFrequency << ", "
							 << setw(8) << t_mpi_step_global * m_pInput->m_nPrintInfoFrequency / m_MPI_Env.size() << ", "
							 << setw(8) << t_mpi_step_max * m_pInput->m_nPrintInfoFrequency << ")"
							 << setw(20) << "Execution time: ("
							 << setw(8) << t_step_min * m_pInput->m_nPrintInfoFrequency << ", "
							 << setw(8) << t_step_global * m_pInput->m_nPrintInfoFrequency / m_MPI_Env.size() << ", "
							 << setw(8) << t_step_max * m_pInput->m_nPrintInfoFrequency << ")" << endl << endl << endl << endl;
			}
		}

		EXTRAE_EVENT(1110013, 0)

		EXTRAE_EVENT(1110012, 1)
		if(m_pInput->m_bHaveResidual)
			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
				RES[idim] = 0.0 ;
			}		

		// Output solution
		if (m_pRunTime->isWriteTime()){
			t10 = clock();

			// Output Solution Residual
			if (m_pInput->m_bHaveResidual && m_pInput->m_bSaveResiduals)
				updateSolutionResidual();

			// Output Solution
			updateSolution();			

			if ( m_pInput->m_bHaveAverage ) {
				updateAverages(time_step-indexStartAverage);
			}

			if ( m_pInput->m_bHaveForces ) {
				if( m_MPI_Env.is_master() )
					output_Forces( false, numberBodies );
			}

			m_MPI_Env.barrier();
			t11 = clock(); t_output += (double) (t11-t10) / CLOCKS_PER_SEC;

			// Write Results
			m_pRunTime->writeResults();
		}
		EXTRAE_EVENT(1110012, 0)

		// Close Loop Event
		EXTRAE_EVENT(1110001, 0)

		GETTIME(nTe);
		dOutputTime += DIFFTIME(nTe, nTs);
	}

	GETTIME(nEndTime);
	dCompleteTime = DIFFTIME(nEndTime, nStartTime);
	t1 = clock();

	EXTRAE_EVENT(1110000, 0)
	EXTRAE_FINI
	EXTRAE_SHUTDOWN

	{
		// Meassure Total Times By All Processors
		double dCompleteTimeG = 0, dCompleteTimeMaxG;
		MPI_Reduce(&dCompleteTime, &dCompleteTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dCompleteTime, &dCompleteTimeMaxG, 1, MPI_DOUBLE, MPI_MAX, m_MPI_Env.get_master(), MPI_COMM_WORLD);

		// Prepare Time
		double dPrepareTimeG = 0;
		MPI_Reduce(&dPrepareTime, &dPrepareTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);

		// RK Loop Time
		double dRKLoopTimeG = 0;
		double dRkPrepareTimeG = 0;
		double dRkMPIWaitTimeG = 0;
		double dRkSetBCTimeG = 0;
		double dRkMinmodTimeG = 0;
		double dRkVisTimeG = 0;
		double dRkCommTimeG = 0;
		double dRkStepTimeG = 0;
		MPI_Reduce(&dRKLoopTime, &dRKLoopTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkPrepareTime, &dRkPrepareTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkMPIWaitTime, &dRkMPIWaitTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkSetBCTime, &dRkSetBCTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkMinmodTime, &dRkMinmodTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkVisTime, &dRkVisTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkCommTime, &dRkCommTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dRkStepTime, &dRkStepTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);

		// Outside Rk Loop
		double dOutsideMPIWaitTimeG = 0;
		double dAdjustCflOrDTTimeG = 0;
		double dAdvanceTimeG = 0;
		double dOutputTimeG = 0;
		MPI_Reduce(&dOutsideMPIWaitTime, &dOutsideMPIWaitTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dAdjustCflOrDTTime, &dAdjustCflOrDTTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dAdvanceTime, &dAdvanceTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);
		MPI_Reduce(&dOutputTime, &dOutputTimeG, 1, MPI_DOUBLE, MPI_SUM, m_MPI_Env.get_master(), MPI_COMM_WORLD);

		if( m_MPI_Env.is_master() )
		{
			const int nRankCount = m_MPI_Env.size();
			const double dCommTimeG = dRkMPIWaitTimeG + dRkCommTimeG + dOutsideMPIWaitTimeG;
			const double dCompTimeG = dCompleteTimeG - dCommTimeG;
			fprintf(stdout, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nTotal Time [%d Ranks]\n", nRankCount);
			fprintf(stdout, "[    ][%.5lf]: Complete Simulation [Maximum: %.5lf]\n", dCompleteTimeG / nRankCount, dCompleteTimeMaxG);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Prepare\n", dPrepareTimeG / dCompleteTimeG * 100, dPrepareTimeG / nRankCount);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Rk Loop\n", dRKLoopTimeG / dCompleteTimeG * 100, dRKLoopTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk Prepare\n", dRkPrepareTimeG / dRKLoopTimeG * 100, dRkPrepareTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk Wait\n", dRkMPIWaitTimeG / dRKLoopTimeG * 100, dRkMPIWaitTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk Comm\n", dRkCommTimeG / dRKLoopTimeG * 100, dRkCommTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk SetBC\n", dRkSetBCTimeG / dRKLoopTimeG * 100, dRkSetBCTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk Vis\n", dRkVisTimeG / dRKLoopTimeG * 100, dRkVisTimeG / nRankCount);
			fprintf(stdout, "\t\t[%.2lf %%] / [%.5lf]: Rk OneStep\n", dRkStepTimeG / dRKLoopTimeG * 100, dRkStepTimeG / nRankCount);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Outside Wait\n", dOutsideMPIWaitTimeG / dCompleteTimeG * 100, dOutsideMPIWaitTimeG / nRankCount);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Adjust CFL-DT\n", dAdjustCflOrDTTimeG / dCompleteTimeG * 100, dAdjustCflOrDTTimeG / nRankCount);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Advance Time\n", dAdvanceTimeG / dCompleteTimeG * 100, dAdvanceTimeG / nRankCount);
			fprintf(stdout, "\t[%.2lf %%] / [%.5lf]: Output Time\n\n", dOutputTimeG / dCompleteTimeG * 100, dOutputTimeG / nRankCount);

			fprintf(stdout, "\n\n");
			fprintf(stdout, "Total: [%.5lf]\n", dCompleteTimeG / nRankCount);
			fprintf(stdout, "MPI  : [%.2lf %%] / [%.5lf]\n", dCommTimeG / dCompleteTimeG * 100, dCommTimeG / nRankCount);
			fprintf(stdout, "Comp : [%.2lf %%] / [%.5lf]\n", dCompTimeG / dCompleteTimeG * 100, dCompTimeG / nRankCount);

		}
	}

	// Disable MPI Profiling
 	MPI_Pcontrol(0);

	double t_total = (double) (t1 - t0) / CLOCKS_PER_SEC;
	double tmp_compute;
	double *t_compute = new double[m_MPI_Env.size()];
	double global_cpu = -100;
	double global_ghost, global_postPro, global_output, global_visc, global_init;
	double global_grads, global_solver, global_compute;
	double global_wait0, global_wait1, global_wait2, global_wait3, global_int0, global_int1, global_int, global_bnd;
	for( size_t ism=0; ism < m_pSolverList.size(); ism++ ){
		t_solver += t_submesh[ism];
	}


	double min_mpi, max_mpi, avg_mpi, t_mpi_local;
	t_mpi_local = t_ghost + t_wait0 + t_wait1 + t_wait2 + t_wait3;

	int myrank = m_MPI_Env.rank();
	std::stringstream filename_local;
	filename_local << "./output/Times/performances_rank" << myrank;
	ifstream fin;
	fin.open( filename_local.str() );
	vector<string> tpVec;
	string tp;
	for ( int i=0; i<2; i++ ){
		getline(fin, tp);
		tpVec.push_back( tp );
	}
	fin.close();
	ofstream fout;
	fout.open( filename_local.str() );
	for ( int i=0; i<tpVec.size(); i++)
		fout << tpVec[i] << endl;
	fout << "Compute =" << setw(15) << std::setprecision(5) << t_solver + t_init + t_visc + t_grads << endl;
	fout << "    - solver:  " << setw(15) << std::setprecision(5) << t_solver << endl;
	fout << "    - init:    " << setw(15) << std::setprecision(5) << t_init << endl;
	fout << "    - visc:    " << setw(15) << std::setprecision(5) << t_visc << endl;
	fout << "    - grad:    " << setw(15) << std::setprecision(5) << t_grads << endl;
	fout << "    - interior:" << setw(15) << std::setprecision(5) << t_int << endl;
	fout << "    - boundary:" << setw(15) << std::setprecision(5) << t_bnd << endl;
    fout << "MPI comm =" << setw(15) << std::setprecision(5) << t_mpi_local << endl;
    fout << "    - ghost   :" << setw(15) << std::setprecision(5) << t_ghost << endl;
    fout << "    - wait0   :" << setw(15) << std::setprecision(5) << t_wait0 << endl;
	fout << "    - wait1   :" << setw(15) << std::setprecision(5) << t_wait1 << endl;
	fout << "    - wait2   :" << setw(15) << std::setprecision(5) << t_wait2 << endl;
	fout << "    - wait3   :" << setw(15) << std::setprecision(5) << t_wait3 << endl;
	fout << "Post-Process =" << setw(15) << std::setprecision(5) << t_postPro << endl;
	fout << "File output  =" << setw(15) << std::setprecision(5) << t_output << endl;
	fout << endl;
	fout.close();


	MPI_Allreduce( &t_total  , &global_cpu    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_ghost  , &global_ghost  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_postPro, &global_postPro, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_output , &global_output , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_visc   , &global_visc   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_init   , &global_init   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_grads  , &global_grads  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_solver , &global_solver , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	MPI_Allreduce( &t_wait0  , &global_wait0  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_wait1  , &global_wait1  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_wait2  , &global_wait2  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_wait3  , &global_wait3  , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int0   , &global_int0   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int1   , &global_int1   , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_int    , &global_int    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
	MPI_Allreduce( &t_bnd    , &global_bnd    , 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	MPI_Allreduce( &t_mpi_local, &min_mpi, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &t_mpi_local, &max_mpi, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &t_mpi_local, &avg_mpi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

	tmp_compute = t_init + t_grads + t_visc + t_grads;
	global_compute = global_solver + global_init + global_visc + global_grads;
	t_mpi = global_ghost + global_wait0 + global_wait1 + global_wait2 + global_wait3;

	MPI_Gather( &tmp_compute, 1, MPI_DOUBLE, &t_compute[0], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	if( m_MPI_Env.is_master() )
	{
		cout << endl;
		cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
		cout << "Compute =" << setw(15) << std::setprecision(5) << global_compute    / m_MPI_Env.size() << endl;
		cout << "    - solver:  " << setw(15) << std::setprecision(5) << global_solver / m_MPI_Env.size() << endl;
		cout << "    - init:    " << setw(15) << std::setprecision(5) << global_init   / m_MPI_Env.size() << endl;
		cout << "    - visc:    " << setw(15) << std::setprecision(5) << global_visc   / m_MPI_Env.size() << endl;
		cout << "    - grad:    " << setw(15) << std::setprecision(5) << global_grads  / m_MPI_Env.size() << endl;
		cout << "    - interior:" << setw(15) << std::setprecision(5) << global_int    / m_MPI_Env.size() << endl;
		cout << "    - boundary:" << setw(15) << std::setprecision(5) << global_bnd    / m_MPI_Env.size() << endl;
		cout << "MPI comm =" << setw(15) << std::setprecision(5) << t_mpi / m_MPI_Env.size() << endl;
		cout << "    - ghost   :" << setw(15) << std::setprecision(5) << global_ghost  / m_MPI_Env.size() << endl;
		cout << "    - wait0   :" << setw(15) << std::setprecision(5) << global_wait0  / m_MPI_Env.size() << endl;
		cout << "    - wait1   :" << setw(15) << std::setprecision(5) << global_wait1  / m_MPI_Env.size() << endl;
		cout << "    - wait2   :" << setw(15) << std::setprecision(5) << global_wait2  / m_MPI_Env.size() << endl;
		cout << "    - wait3   :" << setw(15) << std::setprecision(5) << global_wait3  / m_MPI_Env.size() << endl;
		cout << "    - avg mpi:" << setw(15) << std::setprecision(5) << avg_mpi / m_MPI_Env.size() << endl;
		cout << "    - min mpi:" << setw(15) << std::setprecision(5) << min_mpi << endl;
		cout << "    - max mpi:" << setw(15) << std::setprecision(5) << max_mpi << endl;
		cout << "Post-Process =" << setw(15) << std::setprecision(5) << global_postPro / m_MPI_Env.size() << endl;
		cout << "File output  =" << setw(15) << std::setprecision(5) << global_output  / m_MPI_Env.size() << endl;
		cout << endl;
		cout << "Total CPU-time   =" << setw(15) << std::setprecision(5) << global_cpu / m_MPI_Env.size() << endl;
		cout << endl;
		cout << "Simulation finished successfully." << endl;
		cout << endl;

		std::stringstream filename_global;
		filename_global << "./output/Times/performances_global";
		ifstream fin;
		fin.open( filename_global.str() );
		tpVec.clear();
		for ( int i=0; i<4; i++ ){
			getline(fin, tp);
			tpVec.push_back( tp );
		}
		fin.close();
		fout.open( filename_global.str() );
		for ( int i=0; i<tpVec.size(); i++)
			fout << tpVec[i] << endl;
		fout << "Compute =" << setw(15) << std::setprecision(5) << global_compute    / m_MPI_Env.size() << endl;
		fout << "    - solver:  " << setw(15) << std::setprecision(5) << global_solver / m_MPI_Env.size() << endl;
		fout << "    - init:    " << setw(15) << std::setprecision(5) << global_init   / m_MPI_Env.size() << endl;
		fout << "    - visc:    " << setw(15) << std::setprecision(5) << global_visc   / m_MPI_Env.size() << endl;
		fout << "    - grad:    " << setw(15) << std::setprecision(5) << global_grads  / m_MPI_Env.size() << endl;
		fout << "    - interior:" << setw(15) << std::setprecision(5) << global_int    / m_MPI_Env.size() << endl;
		fout << "    - boundary:" << setw(15) << std::setprecision(5) << global_bnd    / m_MPI_Env.size() << endl;
		fout << "MPI comm =" << setw(15) << std::setprecision(5) << t_mpi / m_MPI_Env.size() << endl;
		fout << "    - ghost   :" << setw(15) << std::setprecision(5) << global_ghost  / m_MPI_Env.size() << endl;
		fout << "    - wait0   :" << setw(15) << std::setprecision(5) << global_wait0  / m_MPI_Env.size() << endl;
		fout << "    - wait1   :" << setw(15) << std::setprecision(5) << global_wait1  / m_MPI_Env.size() << endl;
		fout << "    - wait2   :" << setw(15) << std::setprecision(5) << global_wait2  / m_MPI_Env.size() << endl;
		fout << "    - wait3   :" << setw(15) << std::setprecision(5) << global_wait3  / m_MPI_Env.size() << endl;
		fout << "    - avg mpi:" << setw(15) << std::setprecision(5) << avg_mpi / m_MPI_Env.size() << endl;
		fout << "    - min mpi:" << setw(15) << std::setprecision(5) << min_mpi << endl;
		fout << "    - max mpi:" << setw(15) << std::setprecision(5) << max_mpi << endl;
		fout << "Post-Process =" << setw(15) << std::setprecision(5) << global_postPro / m_MPI_Env.size() << endl;
		fout << "File output  =" << setw(15) << std::setprecision(5) << global_output  / m_MPI_Env.size() << endl;
		fout << endl;
		fout << "Total CPU-time   =" << setw(15) << std::setprecision(5) << global_cpu / m_MPI_Env.size() << endl;
		fout.close();
	}

	m_MPI_Env.barrier();
	delete[] RES;
}


//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT>
void Mesh<PRECISION, DIM_CNT> :: output_Forces( bool firstRun, size_t numberBodies ){

	int myrank = m_MPI_Env.rank();

	vector<std::stringstream> tmp_filename;
	tmp_filename.resize( numberBodies );
	for ( size_t itype = 0; itype < numberBodies; itype++) {
		tmp_filename[itype] << "./output/postProcess_Forces_" << itype+1;
	}

	vector<std::string> filename;
	filename.resize( numberBodies );
	for ( size_t itype = 0; itype < numberBodies; itype++) {
		filename[itype] = tmp_filename[itype].str();
	}
	ofstream fout;

	// Write forces header
	char sXYZ[3] = {'x', 'y', 'z'};
	ISolver* pBCSolver = m_pSolverList[INDEX_BND_SUBMESH];
	if ( firstRun ) {
		for ( size_t itype = 0; itype < numberBodies; itype++) 
		{
			fout.open( filename[itype] );
			fout << setw(25) << "time";
			for (int nD = 0; nD < DIM_CNT; nD++)
				fout << setw(25) << "Fp," << sXYZ[nD];
			for (int nD = 0; nD < DIM_CNT; nD++)
				fout << setw(25) << "Fv," << sXYZ[nD];
			for (int nD = 0; nD < DIM_CNT; nD++)
				fout << setw(25) << "Fp," << sXYZ[nD] << " + Fv," << sXYZ[nD];
			fout << endl;
			fout.close();
		}
	} else {
		// Write actual forces on each wall bc
		for ( size_t itype = 0; itype < numberBodies; itype++) {
			fout.open(filename[itype], ios_base::app);
			pBCSolver->postProcForcesOutput( m_MPI_Env, fout, itype );
			fout.close();
		}
	}
}
