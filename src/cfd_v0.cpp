#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>
#include <initializer_list>

#include "math.h"
#include "api/cfdv0_solver.h"
#include "api/dictReaderOF.h"
#include "omp.h"


using namespace std;

//***************************************************************************************************
//											Templates
//***************************************************************************************************

//***************************************************************************************************
//***************************************************************************************************
//									Public functions
//***************************************************************************************************
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: allocate( gmsh_mesh *pSubMesh, std::vector<void*> &cells_addr )
{
	// Allocate vectors
	const int nCellCount = pSubMesh->cells.size();
	cells_cfd .resize(nCellCount);
	cells_addr.resize(nCellCount);
	dummy_cell.resize(1);

	// Initialize cfd cells & record address

	for(int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
	{
		// Cell Identification
		cells_cfd[nCIndex].vars.id    = nCIndex;
		cells_cfd[nCIndex].vars.sm_id = this->m_nSubmeshIndex;
		cells_cfd[nCIndex].vars.is_ghost = false;
		cells_cfd[nCIndex].m_nCellIndex = pSubMesh->cells[nCIndex].id;

		// Cell solution address
		cells_addr[nCIndex] = &(cells_cfd[nCIndex].vars);
	}

	// Dummy cell, used for NEIGHBOR_NONE:
	// 		- instead of having a special case for out-of-bounds neighbors, have a "dummy cell" with
	// 			all variables initialized to zero
	// 		- saves a condition statement or extra loop depending on solution chosen
	dummy_var.q_old[0] = 123456.;

	// MPI Rank
	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: deallocate( ){

	// Free
	cells_cfd .clear( );
	dummy_cell.clear( );

	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
		ghost_mpi[ineigh].clear( );
	}
	ghost_mpi.clear ( );

	for( size_t ibound=0; ibound < ghost_bnd.size(); ibound++ ){
		ghost_bnd[ibound].clear( );
	}
	ghost_bnd.clear ( );
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize the CFD solver: parameters, mpi type, initial/boundary conditions, etc...
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize( MPI_env &mpi_env, gmsh_mesh *pSubMesh, int nSolver)
{
	switch (m_pInput->m_nRkStepCount)
	{
		case 0:
		case 1:
			// RK0-1 (Euler)
			 m_dAk.push_back(  0.0 );
			 m_dBk.push_back(  1.0 );
			break;
		case 2:
			// RK2
			m_dAk.push_back(  0.0 );
			m_dAk.push_back( -1.0 );
			m_dBk.push_back(  1.0 );
			m_dBk.push_back(  0.5 );
			break;
		case 3:
		case 4:
			// RK3-4
			m_dAk.push_back(  0.0 );
			m_dAk.push_back( -756391.0 / 934407.0 );
			m_dAk.push_back( -36441873.0 / 15625000.0 );
			m_dAk.push_back( -1953125.0 / 1085297.0 );
			m_dBk.push_back(  8.0 / 141.0 );
			m_dBk.push_back(  6627.0 / 2000.0 );
			m_dBk.push_back(  609375.0 / 1085297.0 );
			m_dBk.push_back(  198961.0 / 526383.0 );
			break;
		case 5:
			// RK5
			m_dAk.push_back(  0.0 );
			m_dAk.push_back( -0.4178904745 );
			m_dAk.push_back( -1.192151694643 );
			m_dAk.push_back( -1.697784692471 );
			m_dAk.push_back( -1.514183444257 );
			m_dBk.push_back(  0.1496590219993 );
			m_dBk.push_back(  0.3792103129999 );
			m_dBk.push_back(  0.8229550293869 );
			m_dBk.push_back(  0.6994504559488 );
			m_dBk.push_back(  0.1530572479681 );
			break;
		default:
			break;
	}

	// Initialize all parameters, mostly geometrical variables for the CFD solver
	init_params(pSubMesh);

	// Initialize MPI environment: only for boundary submesh
	if( pSubMesh->id == INDEX_BND_SUBMESH ){
		initialize_mpi_env( mpi_env, pSubMesh);
	}

	// Setup boundary conditions, not always necessary
	if( pSubMesh->id == INDEX_BND_SUBMESH ){
		init_boundary_conditions( );
	}

	// Setup initial conditions		
	set_init_conditions( );

	// Initialize vars to zero
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		init_array( cells_cfd[icell].vars.rho_grad  );
		init_array( cells_cfd[icell].vars.rhoU_grad );
		init_array( cells_cfd[icell].vars.rhoE_grad );
		init_array( cells_cfd[icell].vars.Rpsi_grad );
		init_array( cells_cfd[icell].vars.c_grad    );
		init_array( cells_cfd[icell].vars.dudx      );
		init_array( cells_cfd[icell].vars.dTdx      );
		init_array( cells_cfd[icell].vars.tauMC     );
	}
	if( nSolver == 2){
		for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
			// M2 gradients AUSM
			init_array( cells_cfd[icell].vars.p_grad    );
			init_array( cells_cfd[icell].vars.U_grad    );
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *pSubMesh)
{
	// Get MPI neighbors
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	// Cell Count
	const int nCellCount = pSubMesh->cells.size();

	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)	
		{
			// Isolate MPI/periodic neighbors
			auto *neigh = &( pSubMesh->cell2neigh[nCIndex][nCellFIndex]);
			switch( neigh->type ){
				case CELL_REGULAR:
				case CELL_SUBMESH:
				case CELL_NONE:
					continue;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					break;
			}

			// Filter out MPI neighbors we already registered
			if( registered_rank[neigh->sm] ){
				continue;
			}

			// Register neighbor
			mpi_neighbors.push_back( neigh->sm );
			registered_rank[neigh->sm] = true;
		}

	mpi_env.set_mpi_neighbors( mpi_neighbors );

	// Should only happen for a serial case without periodic boundary conditions
	if( mpi_neighbors.size() == 0 )
		return;

	const t_mpi_halo_comm nHaloCommType = static_cast<t_mpi_halo_comm>(m_pInput->m_nHaloCommType);
	init_mpi_types( mpi_env );
	mpi_env.set_halo_comm_reqs(2 * mpi_env.mpi_neighbors.size() );	
	mpi_env.m_nCommType = static_cast<t_mpi_comm_type>(m_pInput->m_nCommType);

	// Use different method for one sided packed / split
	if (m_pInput->m_nHaloCommType == MPI_ONESIDED_BLCK || m_pInput->m_nHaloCommType == MPI_ONESIDED_NONB)
	{
		std::vector<void*> pBufferList;
		std::vector<size_t> nBufferElementCountList;
		std::vector<size_t> nBufferElementSizeList;
		switch (m_pInput->m_nCommType)
		{
			case MPI_EXCHANGE_FULL_BND:
				pBufferList.push_back(&cells_cfd[0]);
				nBufferElementCountList.push_back(cells_cfd.size());
				nBufferElementSizeList.push_back(sizeof(cells_cfd[0]));				
				break;
			case MPI_EXCHANGE_PACKED:
				pBufferList.push_back(&m_SendBufList[0]);
				nBufferElementCountList.push_back(m_SendBufList.size());
				nBufferElementSizeList.push_back(sizeof(PRECISION));
				break;
			case MPI_EXCHANGE_SPLIT:
				pBufferList.push_back(&m_SendBufSolvarsList[0]);
				nBufferElementCountList.push_back(m_SendBufSolvarsList.size());
				nBufferElementSizeList.push_back(sizeof(PRECISION));

				pBufferList.push_back(&m_SendBufViscousList[0]);
				nBufferElementCountList.push_back(m_SendBufViscousList.size());
				nBufferElementSizeList.push_back(sizeof(PRECISION));
				break;
			default:
			break;
		}
		mpi_env.set_halo_comm_type_onesided(nHaloCommType, pBufferList, nBufferElementCountList, nBufferElementSizeList);
	}	
	else
		mpi_env.set_halo_comm_type(nHaloCommType, cells_cfd, ghost_mpi);

	// Initialize MPI Persistant send/recv
	if (m_pInput->m_nHaloCommType == MPI_TWOSIDED_PERS)
	{
		// Initialize communcation with each neighbor
		for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) 
		{
			switch (m_pInput->m_nCommType)
			{
				case MPI_EXCHANGE_FULL_BND:
					mpi_env.initializeCompletePersistantSendRecv<PRECISION,DIM_CNT,FACE_CNT>(ineigh, cells_cfd, ghost_mpi);
					break;
				case MPI_EXCHANGE_PACKED:					
					mpi_env.initializePersistantSendRecv<PRECISION>(ineigh, m_nSendBufCountList[ineigh], &m_SendBufList[m_nSendBufOffsetList[ineigh]],
																		    m_nRecvBufCountList[ineigh], &m_RecvBufList[m_nRecvBufOffsetList[ineigh]], 0);
					break;
				case MPI_EXCHANGE_SPLIT:
					mpi_env.initializePersistantSendRecv<PRECISION>(ineigh, m_nSendBufSolvarCountList[ineigh], &m_SendBufSolvarsList[m_nSendBufSolvarOffsetList[ineigh]],
														 					m_nRecvBufSolvarCountList[ineigh], &m_RecvBufSolvarsList[m_nRecvBufSolvarOffsetList[ineigh]], 0);
					mpi_env.initializePersistantSendRecv<PRECISION>(ineigh, m_nSendBufViscousCountList[ineigh], &m_SendBufViscousList[m_nSendBufViscousOffsetList[ineigh]],
														 					m_nRecvBufViscousCountList[ineigh], &m_RecvBufViscousList[m_nRecvBufViscousOffsetList[ineigh]], 1);
					break;				
			}
		}
	}
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::reorder_faces()
{	
	// Temporary Data	
	t_solution_vars<PRECISION, DIM_CNT> *tmpNeighs[FACE_CNT];
	PRECISION tmpS       [FACE_CNT][DIM_CNT];		// Face surface vector
	PRECISION tmpD       [FACE_CNT][DIM_CNT];		// Vector btw cell center to neighbor's center
	PRECISION tmp_weight_linear[FACE_CNT];			// Coefficients for interpolation

	// Reorder the faces of each cell to have those with neighbour required to calculate the flux first
	for (size_t icell=0; icell < cells_cfd.size(); icell++)
	{
		// Valid / Invalid Faces
		int nValidFaceIndex = 0;
		int nInvalidFaceIndex = FACE_CNT - 1;

		// Cell
		auto *cell  = &cells_cfd[icell].vars;
		auto *param = &cells_cfd[icell];

		// Place valid faces first
		for(unsigned iface=0; iface < FACE_CNT; iface++ )
		{			
			// Fetch neighbor solution vars
			auto *neigh = param->neighs[iface];

			// Is Valid
			bool bIsValid = true;

			// Not Valid
			if( neigh->sm_id < cell->sm_id )
				bIsValid = false;				
			else if( neigh->sm_id == cell->sm_id && neigh->id < cell->id )
				bIsValid = false;

			// Dest Index
			const int nDestIndex = bIsValid ? nValidFaceIndex : nInvalidFaceIndex;
			(bIsValid) ? nValidFaceIndex++ : nInvalidFaceIndex--;

			// Copy Data
			tmpNeighs[nDestIndex] = neigh;
			tmp_weight_linear[nDestIndex] = param->weight_linear[iface];
			for (int nD = 0; nD < DIM_CNT; nD++)
			{
				tmpS[nDestIndex][nD] = param->S[iface][nD];
				tmpD[nDestIndex][nD] = param->d[iface][nD];
			}
		}
		param->m_nFaceCount = nValidFaceIndex;

		// Update Cell
		for (int nFIndex = 0; nFIndex < FACE_CNT; nFIndex++)
		{
			param->neighs[nFIndex] = tmpNeighs[nFIndex];
			param->weight_linear[nFIndex] = tmp_weight_linear[nFIndex];
			for (int nD = 0; nD < DIM_CNT; nD++)
			{
				param->S[nFIndex][nD] = tmpS[nFIndex][nD];
				param->d[nFIndex][nD] = tmpD[nFIndex][nD];
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New method
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::
				allocate_ghost_cells( MPI_env &mpi_env,
									  gmsh_mesh *pSubMesh,
									  vector<vector<void*>> &ghost_ptr,
									  vector<vector<void*>> &ghost_ptr_bnd)
{

	if(pSubMesh->id != INDEX_BND_SUBMESH)
		return;

	// Cell Count
	const int nCellCount = pSubMesh->cells.size();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize boundary ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Get maximum index of local boundary condition (not MPI)
	// Note: bound value now starts from 0 (and not 1 as previous)
	int max_bound_local = 0;
	for(int nCIndex = 0; nCIndex < nCellCount; nCIndex++)	
		for(int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++ )
		{
			t_gmsh_neighbor *neigh = &(pSubMesh->cell2neigh[nCIndex][nFIndex]);

			if( neigh->type != CELL_NONE )
				continue;

			if (neigh->bound > max_bound_local)
				max_bound_local = neigh->bound;
		}

	int max_bound;
	MPI_Allreduce( &max_bound_local, &max_bound, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	// Add 1 since our bound starts at 0
	max_bound++;
	
	// Assign a vector for each boundary type, that contains [nCIndex, nFIndex] of the boundary
	// Overwrite the mesh neighbour id to contain the index within that vector
	boundaries.resize(max_bound);
	for(int nCIndex = 0; nCIndex < nCellCount; nCIndex++)	
		for(int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++ )
		{
			t_gmsh_neighbor *neigh = &(pSubMesh->cell2neigh[nCIndex][nFIndex]);

			if( neigh->type != CELL_NONE )
				continue;

			// Note: Changed boundary index to start from 0 and not 1
			int ibound = neigh->bound;// - 1;

			// Note: Changing the neighbour Id will allow us to get the pointer
			// of the corresponding bounding ghost cell during the assign_pointers
			boundaries[ibound].push_back( { (int) nCIndex, (int) nFIndex } );
			neigh->id = boundaries[ibound].size() - 1 ;
		}
	

	// Allocate and initialize ghost cells of each boundary type
	ghost_bnd    .resize( max_bound );
	ghost_ptr_bnd.resize( max_bound );
	for( size_t ibound=0; ibound < max_bound; ibound++ )
	{
		ghost_bnd    [ibound].resize( boundaries[ibound].size() );
		ghost_ptr_bnd[ibound].resize( boundaries[ibound].size() );
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ )
		{
			for( unsigned idim=0; idim < DIM_CNT+2; idim++ )
			{
				ghost_bnd[ibound][icell].vars.q_old[idim] = 0.;
			}
		}
	}

	// Set Identification information for boundary cfd_cells
	// Each local boundary ghost cell will have a unique Id which is bigger than any Id of interior cell
	int nBoundaryId = pSubMesh->cells.size();
	for( size_t ibound=0; ibound < max_bound; ibound++ )
	{
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ )
		{
			ghost_ptr_bnd[ibound][icell] = &( ghost_bnd[ibound][icell].vars );

			ghost_bnd[ibound][icell].vars.sm_id = MAX_MPI_RANKS;
			ghost_bnd[ibound][icell].vars.id = nBoundaryId++;
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Initialize MPI ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++)
		{
			// Isolate MPI/periodic neighbors
			auto *neigh = &( pSubMesh->cell2neigh[nCIndex][nFIndex]);
			switch( neigh->type )
			{
				case CELL_REGULAR:
				case CELL_SUBMESH:
				case CELL_NONE:
					continue;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					break;
			}

			// Filter out MPI neighbors we already registered
			if( registered_rank[neigh->sm] )
				continue;
			
			// Register neighbor
			mpi_neighbors.push_back(neigh->sm);
			registered_rank[neigh->sm] = true;
		}

	// Should only happen for a serial case without periodic boundary conditions
	// 		==> true for now but will never happen once we implement support for ALL BCs (wall, etc...)
	if( mpi_neighbors.size() == 0 )
		return;

	// Initialize mapping from global rank ID to local index
	rank2local.resize( mpi_env.size(), -1 );
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
		rank2local[mpi_neighbors[ineigh]] = ineigh;
	}

	// Count number of ghost cells from each neighbor
	// This tells each one of our neighbors that we will send them ALL of our cells data (and not just the onces with shared boundary!)
	// We will also receive it from each one of our neighbors, so that we may assign arrays to receive this data
	std::vector<MPI_Request> comm_reqs( 2*mpi_neighbors.size() );
	std::vector<int> ghost_size( mpi_neighbors.size(), -1 );
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int neigh_rank = mpi_neighbors[ineigh];

		// Hpath cellls
		int send_tag = 0;
		int recv_tag = 0;

		// Note: This may be very expensive if only part of our boundaries are in MPI boundary
		// Warning: The index within the mpi ghost (defined during assign_pointers) assume all the cells are sent!
		int local_size = (int) cells_cfd.size();

		MPI_Isend( &local_size, 1, MPI_INT, neigh_rank, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh] );
		MPI_Irecv( &ghost_size[ineigh], 1, MPI_INT, neigh_rank, recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh+mpi_neighbors.size()] );
	}
	MPI_Waitall( comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE );

	// Allocate ghost cells
	ghost_mpi.resize( mpi_neighbors.size() );
	ghost_ptr  .resize( mpi_neighbors.size() );
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
			ghost_mpi[ineigh].resize( ghost_size[ineigh] );
			ghost_ptr[ineigh].resize( ghost_size[ineigh] );
	}

	// Store ghost cells' memory addresses
	for( size_t ineigh=0; ineigh < ghost_mpi.size(); ineigh++ ){
		for( size_t icell=0; icell < ghost_mpi[ineigh].size(); icell++ ){
			ghost_ptr[ineigh][icell] = &( ghost_mpi[ineigh][icell].vars );
			ghost_mpi[ineigh][icell].vars.sm_id = MAX_MPI_RANKS;
			ghost_mpi[ineigh][icell].vars.id    = icell;

			// Initialize the conservatives of ghost mpi cells
			for (unsigned idim = 0; idim < DIM_CNT + 2; idim++)
				ghost_mpi[ineigh][icell].vars.q_old[0] = 0.0;			
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// New MPI protocol
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Send / Recv buffer
	int buffer_size = DIM_CNT+2 + DIM_CNT*DIM_CNT + DIM_CNT + DIM_CNT*DIM_CNT + DIM_CNT;

	// Cells for each MPI neighbor
	local_cells_to_send.resize(mpi_neighbors.size());
	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) {
		local_cells_to_send[ineigh].resize(ghost_mpi[ineigh].size(), -1);
	}

	// Cell mapping for each MPI neighbor
	neigh_cells_to_recv.resize(mpi_neighbors.size());
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New method
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: assign_pointers(
				gmsh_mesh *pSubMesh,
				MPI_env &mpi_env,
				std::vector<std::vector<void*>>  cells_ptr,
				std::vector<std::vector<void*>>  ghost_ptr,
				std::vector<std::vector<void*>>  ghost_ptr_bnd )
{

	// Assign pointers
	const int nCellCount = cells_cfd.size();
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++)
		{
			// Neighbour Info
			t_gmsh_neighbor *neigh = &(pSubMesh->cell2neigh[nCIndex][nFIndex]);

			switch (neigh->type)
			{
			// Same Submesh Neighbour
			case CELL_PERIODIC:
			case CELL_REGULAR:
				cells_cfd[nCIndex].neighs[nFIndex] = static_cast<t_solution_vars<PRECISION, DIM_CNT> *>(cells_ptr[this->m_nSubmeshIndex][neigh->id]);
				break;

			// Different Submesh, Same MPI Rank
			case CELL_SUBMESH:
				cells_cfd[nCIndex].neighs[nFIndex] = static_cast<t_solution_vars<PRECISION, DIM_CNT> *>(cells_ptr[neigh->sm][neigh->id]);
				break;

			// Local Boundary
			case CELL_NONE:
				cells_cfd[nCIndex].neighs[nFIndex] = static_cast<t_solution_vars<PRECISION, DIM_CNT> *>(ghost_ptr_bnd[neigh->bound][neigh->id]);
				break;

			// MPI Boundary
			case CELL_MPI:
			case CELL_PERIODIC_MPI:				
				cells_cfd[nCIndex].neighs[nFIndex] = static_cast<t_solution_vars<PRECISION, DIM_CNT> *>(ghost_ptr[rank2local[neigh->sm]][neigh->id]);
				break;
			default:
				break;
			}
		}

	// Get MPI neighbors
	std::vector<int > mpi_neighbors;
	std::vector<bool> registered_rank( mpi_env.size(), false );

	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++)
		{
			// Isolate MPI/periodic neighbors
			auto *neigh = &( pSubMesh->cell2neigh[nCIndex][nFIndex]);
			switch( neigh->type ){
				case CELL_REGULAR:
				case CELL_SUBMESH:
				case CELL_NONE:
					continue;
				case CELL_MPI:
				case CELL_PERIODIC:
				case CELL_PERIODIC_MPI:
					break;
			}

			// Filter out MPI neighbors we already registered
			if( registered_rank[neigh->sm] ){
				continue;
			}

			// Register neighbor
			mpi_neighbors.push_back( neigh->sm );
			registered_rank[neigh->sm] = true;
		}
	

	// Record neighbors cells for each MPI neighbor
	int myrank = mpi_env.rank();

	std::vector<int> ind_local(mpi_neighbors.size(), 0);

	std::vector<std::vector<bool>> recorded_cells(mpi_neighbors.size(), std::vector<bool>(cells_cfd.size(), false));
	
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nFIndex = 0; nFIndex < pSubMesh->cells[nCIndex].faceCount; nFIndex++)
		{
			auto *neigh = &( pSubMesh->cell2neigh[nCIndex][nFIndex]);

			if (neigh->type != CELL_MPI && neigh->type != CELL_PERIODIC_MPI)
				continue;

			int neigh_rank = rank2local[neigh->sm];

			if (recorded_cells[neigh_rank][nCIndex])
				continue;

			local_cells_to_send[neigh_rank][ind_local[neigh_rank]] = nCIndex;

			ind_local[neigh_rank]++;
			recorded_cells[neigh_rank][nCIndex] = true;
		}
	
	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) {
		local_cells_to_send[ineigh].resize(ind_local[ineigh]);
	}

	// Exchange mapping size
	int nSendCellOffset = 0;
	int send_info[2];
	std::vector<MPI_Request> comm_reqs(2*mpi_neighbors.size());
	std::vector<int> ghost_info(mpi_neighbors.size() * 2);
	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) {

		int neigh_rank = mpi_neighbors[ineigh];

		// Tags
		int send_tag = 0;
		int recv_tag = 0;
		
		send_info[0] = local_cells_to_send[ineigh].size();
		send_info[1] = nSendCellOffset;
		MPI_Isend( &send_info, 2, MPI_INT, neigh_rank, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh] );
		MPI_Irecv( &ghost_info[ineigh * 2], 2, MPI_INT, neigh_rank, recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh+mpi_neighbors.size()] );
		nSendCellOffset += local_cells_to_send[ineigh].size();
	}
	MPI_Waitall( comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE );

	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) {
		neigh_cells_to_recv[ineigh].resize(ghost_info[ineigh * 2], -1);
	}

	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) {
		int neigh_rank = mpi_neighbors[ineigh];

		int send_tag = 0;
		int recv_tag = 0;

		int send_size = (int) local_cells_to_send[ineigh].size();
		int recv_size = (int) neigh_cells_to_recv[ineigh].size();

		MPI_Isend(&local_cells_to_send[ineigh][0], send_size, MPI_INT, neigh_rank, send_tag, MPI_COMM_WORLD, &comm_reqs[ineigh]);
		MPI_Irecv(&neigh_cells_to_recv[ineigh][0], recv_size, MPI_INT, neigh_rank, recv_tag, MPI_COMM_WORLD, &comm_reqs[ineigh + mpi_neighbors.size()]);
	}
	MPI_Waitall( comm_reqs.size(), &comm_reqs[0], MPI_STATUS_IGNORE );

	// Resize send and recv buffers to the exact number of cells to exchange with each neighbor
	int buffer_size = DIM_CNT+2 + DIM_CNT*DIM_CNT + DIM_CNT + DIM_CNT*DIM_CNT + DIM_CNT;
	int nSendBufSize = 0;
	int nRecvBufSize = 0;
	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) 
	{
		// Send / Recv Offset
		m_nSendBufOffsetList.push_back(nSendBufSize);
		m_nRecvBufOffsetList.push_back(nRecvBufSize);

		// Send / Recv Count
		const int nSendCount = local_cells_to_send[ineigh].size() * buffer_size;
		const int nRecvCount = neigh_cells_to_recv[ineigh].size() * buffer_size;		
		m_nSendBufCountList.push_back(nSendCount);
		m_nRecvBufCountList.push_back(nRecvCount);
		nSendBufSize += nSendCount;
		nRecvBufSize += nRecvCount;

		// Window Offset
		m_nRecvWindowOffsetList.push_back(ghost_info[ineigh * 2 + 1] * buffer_size);
	}
	
	// Separate MPI comms
	int buffer_size_solvars = DIM_CNT+2;
	int buffer_size_viscous = DIM_CNT*DIM_CNT + DIM_CNT + DIM_CNT*DIM_CNT + DIM_CNT;
	int nSendSolvarsBufSize = 0;
	int nRecvSolvarsBufSize = 0;
	int nSendViscousBufSize = 0;
	int nRecvViscousBufSize = 0;
	for (size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++) 
	{
		// Send / Recv Offset
		m_nSendBufSolvarOffsetList.push_back(nSendSolvarsBufSize);
		m_nRecvBufSolvarOffsetList.push_back(nRecvSolvarsBufSize);
		m_nSendBufViscousOffsetList.push_back(nSendViscousBufSize);
		m_nRecvBufViscousOffsetList.push_back(nRecvViscousBufSize);

		// Send / Recv Count
		const int nSendSolvarsCount = local_cells_to_send[ineigh].size() * buffer_size_solvars;
		const int nRecvSolvarsCount = neigh_cells_to_recv[ineigh].size() * buffer_size_solvars;
		const int nSendViscousCount = local_cells_to_send[ineigh].size() * buffer_size_viscous;
		const int nRecvViscousCount = neigh_cells_to_recv[ineigh].size() * buffer_size_viscous;
		m_nSendBufSolvarCountList.push_back(nSendSolvarsCount);
		m_nRecvBufSolvarCountList.push_back(nRecvSolvarsCount);
		m_nSendBufViscousCountList.push_back(nSendViscousCount);
		m_nRecvBufViscousCountList.push_back(nRecvViscousCount);
		nSendSolvarsBufSize += nSendSolvarsCount;
		nRecvSolvarsBufSize += nRecvSolvarsCount;
		nSendViscousBufSize += nSendViscousCount;
		nRecvViscousBufSize += nRecvViscousCount;

		// Window Offset
		m_nRecvSolvarsWindowOffsetList.push_back(ghost_info[ineigh * 2 + 1] * buffer_size_solvars);
		m_nRecvViscousWindowOffsetList.push_back(ghost_info[ineigh * 2 + 1] * buffer_size_viscous);
	}

	// Allocate Send/Recive Buffers
	m_SendBufList.resize(nSendBufSize);
	m_RecvBufList.resize(nRecvBufSize);
	m_SendBufSolvarsList.resize(nSendSolvarsBufSize);
	m_RecvBufSolvarsList.resize(nRecvSolvarsBufSize);
	m_SendBufViscousList.resize(nSendViscousBufSize);
	m_RecvBufViscousList.resize(nRecvViscousBufSize);
}

//***************************************************************************************************
//***************************************************************************************************
//									PRIVATE FUNCTIONS
//***************************************************************************************************
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize the CFD solver
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_params( gmsh_mesh *pSubMesh )
{
	// Cell Count
	const int nCellCount = cells_cfd.size();

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Hpath cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	cells_cfd.set_start(0);
	cells_cfd.set_end  (nCellCount);

	// Cell's center
	xc.resize(nCellCount, std::vector<PRECISION>( DIM_CNT ));
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nD = 0; nD < DIM_CNT; nD++)
			xc[nCIndex][nD] = pSubMesh->cells[nCIndex].xc[nD];

	// Compute face surface vector
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
	{
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
		{	
			const int nFIndex = pSubMesh->cell2face[nCIndex][nCellFIndex];			
			for (int nD = 0; nD < DIM_CNT; nD++)
				cells_cfd[nCIndex].S[nCellFIndex][nD] = pSubMesh->v_sign[nCIndex][nCellFIndex] * pSubMesh->v_norm[nFIndex][nD];				
			
		}
	}
		
	// Vector from face center to cell's center
	std::vector<std::vector<std::vector<PRECISION>>> dP( cells_cfd.size(), std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
			dP[nCIndex][nCellFIndex].resize( DIM_CNT, 0.0 );

	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
		{
			const int nFIndex = pSubMesh->cell2face[nCIndex][nCellFIndex];			
			for (int nD = 0; nD < DIM_CNT; nD++)
				dP[nCIndex][nCellFIndex][nD] = pSubMesh->cells[nCIndex].xc[nD] - pSubMesh->faces[nFIndex].xf[nD];
		}

	// Vector from cell center to neighbor's center
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
		{
			const int nFIndex = pSubMesh->cell2face[nCIndex][nCellFIndex];			
			auto *neigh = &(pSubMesh->cell2neigh[nCIndex][nCellFIndex]);

			// Local Bounday
			if (neigh->type == CELL_NONE)
			{
				PRECISION S_mag_sqrt = 0.0;
				PRECISION dP_dot_S = 0.0;
				for (int nD = 0; nD < DIM_CNT; nD++)
					S_mag_sqrt += cells_cfd[nCIndex].S[nCellFIndex][nD] * cells_cfd[nCIndex].S[nCellFIndex][nD];

				for (int nD = 0; nD < DIM_CNT; nD++)
					dP_dot_S += dP[nCIndex][nCellFIndex][nD] * cells_cfd[nCIndex].S[nCellFIndex][nD];

				for (int nD = 0; nD < DIM_CNT; nD++)
					cells_cfd[nCIndex].d[nCellFIndex][nD] = 2.0 * abs(dP_dot_S) * cells_cfd[nCIndex].S[nCellFIndex][nD] / S_mag_sqrt;
			}
			else
			{
				for (int nD = 0; nD < DIM_CNT; nD++)
					cells_cfd[nCIndex].d[nCellFIndex][nD] = pSubMesh->v_sign[nCIndex][nCellFIndex] * pSubMesh->v_neigh[nFIndex][nD];
			}
		}

	// Vector from face center to neighbor's center
	std::vector<std::vector<std::vector<PRECISION>>> dN(nCellCount, std::vector<std::vector<PRECISION>> ( FACE_CNT ) );
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
			dN[nCIndex][nCellFIndex].resize(DIM_CNT, 0.0);
			
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
			for (int nD = 0; nD < DIM_CNT; nD++)
				dN[nCIndex][nCellFIndex][nD] = dP[nCIndex][nCellFIndex][nD] + cells_cfd[nCIndex].d[nCellFIndex][nD] ;

	// beta_pos and beta_neg calculations
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		for (int nCellFIndex = 0; nCellFIndex < pSubMesh->cells[nCIndex].faceCount; nCellFIndex++)
		{
			cells_cfd[nCIndex].weight_linear[nCellFIndex] = 0.0;
			PRECISION divider = 0.0;
			for (int nD = 0; nD < DIM_CNT; nD++)
			{
				cells_cfd[nCIndex].weight_linear[nCellFIndex] += cells_cfd[nCIndex].S[nCellFIndex][nD] * dN[nCIndex][nCellFIndex][nD] ;
				divider += ( - cells_cfd[nCIndex].S[nCellFIndex][nD] * dP[nCIndex][nCellFIndex][nD] + cells_cfd[nCIndex].S[nCellFIndex][nD] * dN[nCIndex][nCellFIndex][nD] );
			}

			cells_cfd[nCIndex].weight_linear[nCellFIndex] /= divider ;
		}

	// Cell's volume
	for (int nCIndex = 0; nCIndex < nCellCount; nCIndex++)
		cells_cfd[nCIndex].vars.vol_inv = ONE / pSubMesh->cells_vol[nCIndex];	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initial conditions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_init_conditions( )
{
	// Initialize Cells
	for(int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
	{
		for (int nD = 0; nD < DIM_CNT; nD++)
		{
			cells_cfd[nCIndex].vars.q_old[nD] = ZERO;
			cells_cfd[nCIndex].vars.delta_q[nD] = ZERO;
		}
		cells_cfd[nCIndex].sponge_sigma = 0.0;
	}

	m_dpInf = m_pInput->m_dPInf;
	m_dTInf = m_pInput->m_dTInf;
	for (int nD = 0; nD < DIM_CNT; nD++)
		m_dUInf[nD] = m_pInput->m_dUInf[nD];
	m_dMu0 = m_pInput->m_dMu0;
	m_dMolWeight = m_pInput->m_dMolWeight;
	m_dCp = m_pInput->m_dCp;	
	const PRECISION dLs = m_pInput->m_dLs;
	const PRECISION dPrandtl = m_pInput->m_dPrandtl;
	
	// Calculate values from input Data
	m_dRgas = m_dRuniversal / m_dMolWeight * 1000;
	const PRECISION dCv = m_dCp - m_dRgas;
	m_dGamma = m_dCp / dCv;
	m_drhoInf = m_dpInf /  (m_dRgas * m_dTInf);
	m_dGammaMinusOne = m_dGamma - ONE;
	m_dGamma_m_one_inv = ONE / m_dGammaMinusOne;
	PRECISION dUMagInf2 = 0;
	for (int nD = 0; nD < DIM_CNT; nD++)
		dUMagInf2 += m_dUInf[nD] * m_dUInf[nD];
	m_dEInf = m_dpInf / (m_drhoInf * m_dGammaMinusOne) + 0.5 * dUMagInf2;
	m_dRgas_inv = ONE / m_dRgas;
	m_dPr_inv = ONE / dPrandtl;

	// Cell Index List
	const std::vector<int> nCellIndexList = m_pMeshReader->getSubmeshCellIndexList(m_nSubmeshIndex);
	
	// Sponge Layer Initialization
	if (dLs > 0)
	{
		// Extract additional information from input
		const PRECISION dK = m_pInput->m_dK;
		const PRECISION dMach = m_pInput->m_dMach;

		// Sigma0
		const PRECISION dSigma0 = ( 3.0 * (1.0 - dMach * dMach) / ( dK * dLs ) );

		// Distance
		// read always from timestep = 0
		const int nAlphaField = m_pMeshReader->readScalarField("alpha");
		
		// Get Distance of each cell from the boundary outlet
		for (int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
		{
			const double dDist = m_pMeshReader->getCellScalar(nCellIndexList[nCIndex], nAlphaField);
			if (dDist < dLs)
				cells_cfd[nCIndex].sponge_sigma = dSigma0 * pow((dLs - dDist) / dLs, 2.0 );
		}
	}

	// Read Scalar Fields (for simulation starting time)
	const int nPField   = m_pMeshReader->readScalarField("p");
	const int nTField   = m_pMeshReader->readScalarField("T");
	const int nUVWField = m_pMeshReader->readVectorField("U");	

	// Initialize Conservatives
	double dUVW[3];	
	for (unsigned int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
	{
		// Cell Index
		const int nCellIndex = nCellIndexList[nCIndex];

		// Velocity
		m_pMeshReader->getCellVector(nCellIndex, nUVWField, dUVW);

		// Pressure, Temperature
		const PRECISION dPressure = m_pMeshReader->getCellScalar(nCellIndex, nPField);
		const PRECISION dTemperature = m_pMeshReader->getCellScalar(nCellIndex, nTField);

		// Velocity Magnitude
		PRECISION dUMag2 = 0;
		for (int nD = 0; nD < DIM_CNT; nD++)
			dUMag2 += dUVW[nD] * dUVW[nD];

		// rho = p/(R * T)
		const PRECISION dRho = dPressure /  (m_dRgas * dTemperature);

		// Energy
		const PRECISION dEnergy = dPressure / (dRho * m_dGammaMinusOne) + 0.5 * dUMag2;

		// Store To cells_cfd		
		cells_cfd[nCIndex].vars.q_old[0] = dRho;
		for (int nD = 0; nD < DIM_CNT; nD++)
			cells_cfd[nCIndex].vars.q_old[nD + 1] = dRho * dUVW[nD];
		cells_cfd[nCIndex].vars.q_old[DIM_CNT+1] = dRho * dEnergy;
	}

}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initial boundary conditions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_boundary_conditions( )
{
	const int nBoundaryCount = m_pMeshReader->getBoundaryCount();
	for (int nBCIndex = 0; nBCIndex < nBoundaryCount; nBCIndex++)
	{			
		// BC Type
		const std::string sBCType = m_pMeshReader->getBoundaryType(nBCIndex);
		const std::string sBCName = m_pMeshReader->getBoundaryName(nBCIndex);		

		// Wall
		if (sBCType.compare("wall") == 0)
			m_nWallBCList.push_back(nBCIndex);

		// Farfield
		if (sBCType.compare("farfield") == 0)
			m_nFarfieldBCList.push_back(nBCIndex);
		
		// Inlet / Outlet
		if (sBCType.compare("patch") == 0)
		{
			if (sBCName.compare("inlet") == 0 ||
			    sBCName.compare("Inlet") == 0 ||
			    sBCName.compare("Inflow") == 0 ||
			    sBCName.compare("inflow") == 0)
				m_nInletBCList.push_back(nBCIndex);
			if (sBCName.compare("outlet") == 0 ||
			    sBCName.compare("Outlet") == 0 ||
			    sBCName.compare("Outflow") == 0 ||
			    sBCName.compare("outflow") == 0)
				m_nOutletBCList.push_back(nBCIndex);
		}			
	}
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Set boundary conditions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: set_boundary_conditions( ){

	// WALL
	for ( size_t itype = 0; itype < m_nWallBCList.size(); itype++) {
		for ( size_t icell = 0; icell < ghost_bnd[m_nWallBCList[itype]].size(); icell++ ) {
			int bnd_id = boundaries[m_nWallBCList[itype]][icell][0];

			ghost_bnd[m_nWallBCList[itype]][icell].vars.q_old[0] = cells_cfd[bnd_id].vars.q_old[0];
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nWallBCList[itype]][icell].vars.q_old[idim+1] = - cells_cfd[bnd_id].vars.q_old[idim+1];
			}
			ghost_bnd[m_nWallBCList[itype]][icell].vars.q_old[DIM_CNT+1] = cells_cfd[bnd_id].vars.q_old[DIM_CNT+1];

			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nWallBCList[itype]][icell].vars.rho_grad[idim] = - cells_cfd[bnd_id].vars.rho_grad[idim];
				for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
					ghost_bnd[m_nWallBCList[itype]][icell].vars.rhoU_grad[idim][jdim] = cells_cfd[bnd_id].vars.rhoU_grad[idim][jdim];
					ghost_bnd[m_nWallBCList[itype]][icell].vars.dudx[idim][jdim]      = cells_cfd[bnd_id].vars.dudx[idim][jdim];
					ghost_bnd[m_nWallBCList[itype]][icell].vars.tauMC[idim][jdim]     = cells_cfd[bnd_id].vars.tauMC[idim][jdim];
				}
				ghost_bnd[m_nWallBCList[itype]][icell].vars.rhoE_grad[idim] =   cells_cfd[bnd_id].vars.rhoE_grad[idim];
				ghost_bnd[m_nWallBCList[itype]][icell].vars.Rpsi_grad[idim] = - cells_cfd[bnd_id].vars.Rpsi_grad[idim];
				ghost_bnd[m_nWallBCList[itype]][icell].vars.c_grad[idim]    = - cells_cfd[bnd_id].vars.c_grad[idim];

				ghost_bnd[m_nWallBCList[itype]][icell].vars.dTdx[idim]      = - cells_cfd[bnd_id].vars.dTdx[idim];
			}
		}
	}

	// INLET
	for ( size_t itype = 0; itype < m_nInletBCList.size(); itype++) {
		for ( size_t icell = 0; icell < ghost_bnd[m_nInletBCList[itype]].size(); icell++ ) {
			int bnd_id = boundaries[m_nInletBCList[itype]][icell][0];

			PRECISION rhoInt = cells_cfd[bnd_id].vars.q_old[0];
			PRECISION UvecInt[DIM_CNT];
			PRECISION Umag_sqrtInt = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				UvecInt[idim] = cells_cfd[bnd_id].vars.q_old[idim+1] / rhoInt;
				Umag_sqrtInt += UvecInt[idim] * UvecInt[idim];
			}
			PRECISION pInt = ( cells_cfd[bnd_id].vars.q_old[DIM_CNT+1] - 0.5 * Umag_sqrtInt * rhoInt ) * m_dGammaMinusOne;
			PRECISION TInt = pInt * m_dRgas_inv / rhoInt;

			PRECISION TExt = m_dTInf;
			PRECISION rhoExt = rhoInt * TInt / TExt;
			PRECISION UvecExt[DIM_CNT];
			PRECISION Umag_sqrtExt = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				UvecExt[idim] = m_dUInf[idim];
				Umag_sqrtExt += UvecExt[idim] * UvecExt[idim];
			}
			PRECISION pExt = pInt;
			PRECISION EExt = pExt / (rhoExt * m_dGammaMinusOne) + 0.5 * Umag_sqrtExt;

			ghost_bnd[m_nInletBCList[itype]][icell].vars.q_old[0] = rhoExt;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nInletBCList[itype]][icell].vars.q_old[idim+1] = rhoExt * UvecExt[idim];
			}
			ghost_bnd[m_nInletBCList[itype]][icell].vars.q_old[DIM_CNT+1] = rhoExt * EExt;

			// We want that all gradients on faces will be equal to zero (zero diffusion on farfield boundaries)
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nInletBCList[itype]][icell].vars.rho_grad[idim] = - cells_cfd[bnd_id].vars.rho_grad[idim];
				for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
					ghost_bnd[m_nInletBCList[itype]][icell].vars.rhoU_grad[idim][jdim] = - cells_cfd[bnd_id].vars.rhoU_grad[idim][jdim];
					ghost_bnd[m_nInletBCList[itype]][icell].vars.dudx[idim][jdim]      = - cells_cfd[bnd_id].vars.dudx[idim][jdim];
					ghost_bnd[m_nInletBCList[itype]][icell].vars.tauMC[idim][jdim]     = - cells_cfd[bnd_id].vars.tauMC[idim][jdim];
				}
				ghost_bnd[m_nInletBCList[itype]][icell].vars.rhoE_grad[idim] = - cells_cfd[bnd_id].vars.rhoE_grad[idim];
				ghost_bnd[m_nInletBCList[itype]][icell].vars.Rpsi_grad[idim] = - cells_cfd[bnd_id].vars.Rpsi_grad[idim];
				ghost_bnd[m_nInletBCList[itype]][icell].vars.c_grad[idim]    = - cells_cfd[bnd_id].vars.c_grad[idim];

				ghost_bnd[m_nInletBCList[itype]][icell].vars.dTdx[idim]      = - cells_cfd[bnd_id].vars.dTdx[idim];
			}
		}
	}

	// OUTLET
	for ( size_t itype = 0; itype < m_nOutletBCList.size(); itype++) {
		for ( size_t icell = 0; icell < ghost_bnd[m_nOutletBCList[itype]].size(); icell++ ) {
			int bnd_id = boundaries[m_nOutletBCList[itype]][icell][0];
		
			PRECISION rhoInt = cells_cfd[bnd_id].vars.q_old[0];
			PRECISION UvecInt[DIM_CNT];
			PRECISION Umag_sqrtInt = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				UvecInt[idim] = cells_cfd[bnd_id].vars.q_old[idim+1] / rhoInt;
				Umag_sqrtInt += UvecInt[idim] * UvecInt[idim];
			}
			PRECISION pInt = ( cells_cfd[bnd_id].vars.q_old[DIM_CNT+1] - 0.5 * Umag_sqrtInt * rhoInt ) * m_dGammaMinusOne;

			PRECISION pExt = m_dpInf;
			PRECISION rhoExt = rhoInt * pExt / pInt;
			PRECISION UvecExt[DIM_CNT];
			PRECISION Umag_sqrtExt = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				UvecExt[idim] = UvecInt[idim];
				Umag_sqrtExt += UvecExt[idim] * UvecExt[idim];
			}
			PRECISION EExt = pExt / (rhoExt * m_dGammaMinusOne) + 0.5 * Umag_sqrtExt;
			
			ghost_bnd[m_nOutletBCList[itype]][icell].vars.q_old[0] = rhoExt;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nOutletBCList[itype]][icell].vars.q_old[idim+1] = rhoExt * UvecExt[idim];
			}
			ghost_bnd[m_nOutletBCList[itype]][icell].vars.q_old[DIM_CNT+1] = rhoExt * EExt;

			// We want that all gradients on faces will be equal to zero (zero diffusion on farfield boundaries)
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nOutletBCList[itype]][icell].vars.rho_grad[idim] = - cells_cfd[bnd_id].vars.rho_grad[idim];
				for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
					ghost_bnd[m_nOutletBCList[itype]][icell].vars.rhoU_grad[idim][jdim] = - cells_cfd[bnd_id].vars.rhoU_grad[idim][jdim];
					ghost_bnd[m_nOutletBCList[itype]][icell].vars.dudx[idim][jdim]      = - cells_cfd[bnd_id].vars.dudx[idim][jdim];
					ghost_bnd[m_nOutletBCList[itype]][icell].vars.tauMC[idim][jdim]     = - cells_cfd[bnd_id].vars.tauMC[idim][jdim];
				}
				ghost_bnd[m_nOutletBCList[itype]][icell].vars.rhoE_grad[idim] = - cells_cfd[bnd_id].vars.rhoE_grad[idim];
				ghost_bnd[m_nOutletBCList[itype]][icell].vars.Rpsi_grad[idim] = - cells_cfd[bnd_id].vars.Rpsi_grad[idim];
				ghost_bnd[m_nOutletBCList[itype]][icell].vars.c_grad[idim]    = - cells_cfd[bnd_id].vars.c_grad[idim];

				ghost_bnd[m_nOutletBCList[itype]][icell].vars.dTdx[idim]      = - cells_cfd[bnd_id].vars.dTdx[idim];
			}
		}
	}

	// FAR-FIELD
	for ( size_t itype = 0; itype < m_nFarfieldBCList.size(); itype++) {
		for ( size_t icell = 0; icell < ghost_bnd[m_nFarfieldBCList[itype]].size(); icell++ ) {
			int bnd_id_cell = boundaries[m_nFarfieldBCList[itype]][icell][0];
			int bnd_id_face = boundaries[m_nFarfieldBCList[itype]][icell][1];

			PRECISION Udirection[DIM_CNT] = { (PRECISION) cos(m_dAoA * M_PI / 180.0), (PRECISION) sin(m_dAoA * M_PI / 180.0) };
			PRECISION T = 0.0;
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				T += Udirection[idim] * cells_cfd[bnd_id_cell].S[idim][bnd_id_face];
			}

			if ( T <= 0.0 ) {		// INLET

				PRECISION rhoInt = cells_cfd[bnd_id_cell].vars.q_old[0];
				PRECISION UvecInt[DIM_CNT];
				PRECISION Umag_sqrtInt = 0.0;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					UvecInt[idim] = cells_cfd[bnd_id_cell].vars.q_old[idim+1] / rhoInt;
					Umag_sqrtInt += UvecInt[idim] * UvecInt[idim];
				}
				PRECISION pInt = ( cells_cfd[bnd_id_cell].vars.q_old[DIM_CNT+1] - 0.5 * Umag_sqrtInt * rhoInt ) * m_dGammaMinusOne;
				PRECISION TInt = pInt * m_dRgas_inv / rhoInt;

				PRECISION TExt = m_dTInf;
				PRECISION rhoExt = rhoInt * TInt / TExt;
				PRECISION UvecExt[DIM_CNT];
				PRECISION Umag_sqrtExt = 0.0;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					UvecExt[idim] = m_dUInf[idim];
					Umag_sqrtExt += UvecExt[idim] * UvecExt[idim];
				}
				PRECISION pExt = pInt;
				PRECISION EExt = pExt / (rhoExt * m_dGammaMinusOne) + 0.5 * Umag_sqrtExt;

				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[0] = rhoExt;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[idim+1] = rhoExt * UvecExt[idim];
				}
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[DIM_CNT+1] = rhoExt * EExt;

			} else {				 // OUTLET

				PRECISION rhoInt = cells_cfd[bnd_id_cell].vars.q_old[0];
				PRECISION UvecInt[DIM_CNT];
				PRECISION Umag_sqrtInt = 0.0;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					UvecInt[idim] = cells_cfd[bnd_id_cell].vars.q_old[idim+1] / rhoInt;
					Umag_sqrtInt += UvecInt[idim] * UvecInt[idim];
				}
				PRECISION pInt = ( cells_cfd[bnd_id_cell].vars.q_old[DIM_CNT+1] - 0.5 * Umag_sqrtInt * rhoInt ) * m_dGammaMinusOne;

				PRECISION pExt = m_dpInf;
				PRECISION rhoExt = rhoInt * pExt / pInt;
				PRECISION UvecExt[DIM_CNT];
				PRECISION Umag_sqrtExt = 0.0;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					UvecExt[idim] = UvecInt[idim];
					Umag_sqrtExt += UvecExt[idim] * UvecExt[idim];
				}
				PRECISION EExt = pExt / (rhoExt * m_dGammaMinusOne) + 0.5 * Umag_sqrtExt;
				
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[0] = rhoExt;
				for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
					ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[idim+1] =  rhoExt * UvecExt[idim];
				}
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.q_old[DIM_CNT+1] = rhoExt * EExt;

			}

			// We want that all gradients on faces will be equal to zero (zero diffusion on farfield boundaries)
			for ( size_t idim = 0; idim < DIM_CNT; idim++ ) {
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.rho_grad[idim] = - cells_cfd[bnd_id_cell].vars.rho_grad[idim];
				for ( size_t jdim = 0; jdim < DIM_CNT; jdim++ ) {
					ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.rhoU_grad[idim][jdim] = - cells_cfd[bnd_id_cell].vars.rhoU_grad[idim][jdim];
					ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.dudx[idim][jdim]      = - cells_cfd[bnd_id_cell].vars.dudx[idim][jdim];
					ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.tauMC[idim][jdim]     = - cells_cfd[bnd_id_cell].vars.tauMC[idim][jdim];
				}
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.rhoE_grad[idim] = - cells_cfd[bnd_id_cell].vars.rhoE_grad[idim];
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.Rpsi_grad[idim] = - cells_cfd[bnd_id_cell].vars.Rpsi_grad[idim];
				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.c_grad[idim]    = - cells_cfd[bnd_id_cell].vars.c_grad[idim];

				ghost_bnd[m_nFarfieldBCList[itype]][icell].vars.dTdx[idim]      = - cells_cfd[bnd_id_cell].vars.dTdx[idim];
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Initialize MPI types
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: init_mpi_types( MPI_env &mpi_env ){

	if( mpi_env.size() == 1 ){
		return;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI Datatype for t_solution_vars
	// 		- 15 fields total
	// 		- without sm_id, id, and D2U
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	const int solvars_nfields = 16;
	MPI_Aint  solvars_disps[solvars_nfields];
	int       solvars_blocklens[] = { DIM_CNT+2, DIM_CNT+2, DIM_CNT, DIM_CNT*DIM_CNT,
									  DIM_CNT, DIM_CNT, DIM_CNT, DIM_CNT*DIM_CNT, DIM_CNT,
									  DIM_CNT*DIM_CNT, DIM_CNT, 1, DIM_CNT+2, 1, 1, 1 };

	MPI_Datatype mpi_precision;
	if( std::is_same< PRECISION, float >::value ){
		mpi_precision = MPI_FLOAT;
	}else if( std::is_same< PRECISION, double >::value ){
		mpi_precision = MPI_DOUBLE;
	}else{
		MPI_Abort( MPI_COMM_WORLD, 254 );
	}

	MPI_Datatype solvars_types[] = { mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, mpi_precision, mpi_precision,
									 mpi_precision, mpi_precision, mpi_precision, MPI_UNSIGNED, MPI_UNSIGNED_SHORT, MPI_LOGICAL };

	t_solution_vars<PRECISION, DIM_CNT> tmp_solvars;

	solvars_disps[ 0] = (char*) &tmp_solvars.q_old     - (char*) &tmp_solvars;	
	solvars_disps[ 1] = (char*) &tmp_solvars.delta_q   - (char*) &tmp_solvars;
	solvars_disps[ 2] = (char*) &tmp_solvars.rho_grad  - (char*) &tmp_solvars;
	solvars_disps[ 3] = (char*) &tmp_solvars.rhoU_grad - (char*) &tmp_solvars;
	solvars_disps[ 4] = (char*) &tmp_solvars.rhoE_grad - (char*) &tmp_solvars;
	solvars_disps[ 5] = (char*) &tmp_solvars.Rpsi_grad - (char*) &tmp_solvars;
	solvars_disps[ 6] = (char*) &tmp_solvars.c_grad    - (char*) &tmp_solvars;
	solvars_disps[ 7] = (char*) &tmp_solvars.dudx      - (char*) &tmp_solvars;
	solvars_disps[ 8] = (char*) &tmp_solvars.dTdx      - (char*) &tmp_solvars;
	solvars_disps[ 9] = (char*) &tmp_solvars.tauMC     - (char*) &tmp_solvars;
	solvars_disps[10] = (char*) &tmp_solvars.sigmaU    - (char*) &tmp_solvars;
	solvars_disps[11] = (char*) &tmp_solvars.vol_inv   - (char*) &tmp_solvars;
	solvars_disps[12] = (char*) &tmp_solvars.RES       - (char*) &tmp_solvars;
	solvars_disps[13] = (char*) &tmp_solvars.id        - (char*) &tmp_solvars;
	solvars_disps[14] = (char*) &tmp_solvars.sm_id     - (char*) &tmp_solvars;
	solvars_disps[15] = (char*) &tmp_solvars.is_ghost  - (char*) &tmp_solvars;

	MPI_Datatype type_solvars;
	MPI_Type_create_struct ( solvars_nfields, solvars_blocklens, solvars_disps, solvars_types, &type_solvars );
	MPI_Type_commit( &type_solvars );

	// MPI Datatype for the CFDv0 cell
	const int cfdcell_nfields = 6;
	MPI_Aint  cfdcell_disps[cfdcell_nfields];
	int       cfdcell_blocklens[] = { 1,				// solution vars struct
									  FACE_CNT,			// array of pointers to neighbors' sol vars
									  FACE_CNT*DIM_CNT,	// S
									  FACE_CNT*DIM_CNT,	// d
									  FACE_CNT,			// weight_linear
									  1,				// sponge_sigma
									 };

	MPI_Datatype cfdcell_types[] = { type_solvars,			// t_solution_vars
									 mpi_precision,			// neighbors
									 mpi_precision,			// S
									 mpi_precision,			// d
									 mpi_precision,			// weight_linear
									 mpi_precision,			// sponge_sigma
									};

	// Compute displacements in the CFDv0 cell
	CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT> tmp_cfdcell;

	cfdcell_disps[ 0] = (char*) &tmp_cfdcell.vars          - (char*) &tmp_cfdcell;
	cfdcell_disps[ 1] = (char*) &tmp_cfdcell.neighs        - (char*) &tmp_cfdcell;
	cfdcell_disps[ 2] = (char*) &tmp_cfdcell.S             - (char*) &tmp_cfdcell;
	cfdcell_disps[ 3] = (char*) &tmp_cfdcell.d             - (char*) &tmp_cfdcell;
	cfdcell_disps[ 4] = (char*) &tmp_cfdcell.weight_linear - (char*) &tmp_cfdcell;
	cfdcell_disps[ 5] = (char*) &tmp_cfdcell.sponge_sigma  - (char*) &tmp_cfdcell;

	MPI_Aint total_memory_disp = (char*) &cells_cfd[1] - (char*) &cells_cfd[0];

	MPI_Datatype tmp_type, type_cfdcell;
	MPI_Type_create_struct ( cfdcell_nfields, cfdcell_blocklens, cfdcell_disps, cfdcell_types, &tmp_type );
	MPI_Type_create_resized( tmp_type, 0, total_memory_disp, &type_cfdcell );
	MPI_Type_commit( &type_cfdcell );

	// Commit type to MPI env
	mpi_env.set_type_cfdcell( type_cfdcell );
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// cell->delta_q[i] = cell->Res[i] = 0.0
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: prepare_for_timestep(){

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		for( unsigned idim=0; idim < DIM_CNT + 2; idim++ ){
			cells_cfd[icell].vars.delta_q[idim] = cells_cfd[icell].vars.RES[idim] = 0.0;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Make preparations for the time-step integration
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: prepare_for_RKstep( int rk_step, int nSolver ){

	// Copy new ----> old
	const size_t nCellCount = cells_cfd.size();

	for( size_t icell=0; icell < nCellCount; icell++ ){
		for( unsigned idim=0; idim < DIM_CNT+2; idim++ ){
			cells_cfd[icell].vars.delta_q[idim] *= m_dAk[rk_step];
		}

		// Initialize vars to zero
		init_array( cells_cfd[icell].vars.rho_grad  );
		init_array( cells_cfd[icell].vars.rhoU_grad );
		init_array( cells_cfd[icell].vars.rhoE_grad );
		init_array( cells_cfd[icell].vars.Rpsi_grad );
		init_array( cells_cfd[icell].vars.c_grad    );
		init_array( cells_cfd[icell].vars.dudx      );
		init_array( cells_cfd[icell].vars.dTdx      );
		init_array( cells_cfd[icell].vars.tauMC     );
	}

	if( nSolver == 2) {
		for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
			// M2 gradients AUSM
			init_array( cells_cfd[icell].vars.p_grad    );
			init_array( cells_cfd[icell].vars.U_grad    );
		}
	}

	if( m_nSubmeshIndex != INDEX_BND_SUBMESH )
		return;

	// Boundary ghost cells
	for( size_t ibound=0; ibound < ghost_bnd.size(); ibound++ ){
		for( size_t icell=0; icell < ghost_bnd[ibound].size(); icell++ ){
			ghost_bnd[ibound][icell].vars.is_ghost = true;
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate needed gradients for cell to face interpolation schemes, only for M2-AUSM components
// 	- only linear interpolation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_gradients_M2AUSM( MPI_env &mpi_env)
{
	// M2 AUSM
	PRECISION UU[DIM_CNT]; // Compute U before Ugrad
	PRECISION cell_surf_over_vol[DIM_CNT];
	PRECISION adjc_surf_over_vol[DIM_CNT];
	PRECISION rhoU[DIM_CNT];

	// Run on all cells
	for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){
		// Get Cell & Vars
		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		PRECISION cell_Rpsi = compute_Rpsi( vars->q_old );
		PRECISION cell_c    = sqrt( m_dGamma * cell_Rpsi );

		// Check all faces
		for( unsigned iface=0; iface < cell->m_nFaceCount; iface++ )
		{
			// Get Neighbour
			auto *neigh = cells_cfd[icell].neighs[iface];

			// Neighbor primitive vars
			PRECISION adjc_Rpsi = compute_Rpsi( neigh->q_old );
			PRECISION adjc_c    = sqrt( m_dGamma * adjc_Rpsi );

			// Evaluate UU
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				UU[idim] = INTERP_LINEAR( cell->weight_linear[iface], vars ->q_old[idim+1]/vars->q_old[0],
																	neigh->q_old[idim+1]/neigh->q_old[0] );
			}

			// Fluxes
			PRECISION rho = INTERP_LINEAR( cell->weight_linear[iface],	vars ->q_old[0],
																		neigh->q_old[0] );
			

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoU[idim] = INTERP_LINEAR( cell->weight_linear[iface], vars ->q_old[idim+1],
																		neigh->q_old[idim+1] );
			}

			PRECISION rhoE = INTERP_LINEAR( cell->weight_linear[iface], vars ->q_old[DIM_CNT+1],
																		neigh->q_old[DIM_CNT+1] );

			PRECISION Rpsi = INTERP_LINEAR( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			PRECISION c    = INTERP_LINEAR( cell->weight_linear[iface], cell_c,    adjc_c    );

			// Cell and neighbor (surface / vol );			
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =   cell->S[iface][idim] * vars ->vol_inv;
				adjc_surf_over_vol[idim] = - cell->S[iface][idim] * neigh->vol_inv;
			}

			// Calculate cell_p and neigh_p
			PRECISION r = vars->q_old[0];
			PRECISION E = vars->q_old[DIM_CNT+1] / r;
			PRECISION vmag = 0.0;
			for(unsigned idim = 0; idim < DIM_CNT; idim++){
				vmag += (vars->q_old[idim+1]/r) * (vars->q_old[idim+1]/r);
			}
			PRECISION cell_p = r * m_dGammaMinusOne * (E-HALF*vmag);
			r = neigh->q_old[0];
			E = neigh->q_old[DIM_CNT+1] / r;
			vmag = 0.0;
			for(unsigned idim = 0; idim < DIM_CNT; idim++){
				vmag += (neigh->q_old[idim+1]/r) * (neigh->q_old[idim+1]/r);
			}
			PRECISION neigh_p = r * m_dGammaMinusOne * (E-HALF*vmag);

			// Average Pressure
			PRECISION p = INTERP_LINEAR( cell->weight_linear[iface], cell_p, neigh_p);

			// Cell & Neighbor Gradient
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				//p_grad
				vars->p_grad [idim] += p  * cell_surf_over_vol[idim];
				neigh->p_grad [idim] += p  * adjc_surf_over_vol[idim];

				// U_grad
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars->U_grad[jdim][idim] += UU[jdim] * cell_surf_over_vol[idim];
					neigh->U_grad[jdim][idim] += UU[jdim] * adjc_surf_over_vol[idim];
				}

				//rho_grad
				vars->rho_grad [idim] += rho  * cell_surf_over_vol[idim];
				neigh->rho_grad [idim] += rho  * adjc_surf_over_vol[idim];

				//rhoU_grad
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars->rhoU_grad[jdim][idim] += rhoU[jdim] * cell_surf_over_vol[idim];
					neigh->rhoU_grad[jdim][idim] += rhoU[jdim] * adjc_surf_over_vol[idim];
				}

				//rhoE_grad
				vars->rhoE_grad[idim] += rhoE * cell_surf_over_vol[idim];
				neigh->rhoE_grad[idim] += rhoE * adjc_surf_over_vol[idim];

				//Rpsi_grad
				vars->Rpsi_grad[idim] += Rpsi * cell_surf_over_vol[idim];
				neigh->Rpsi_grad[idim] += Rpsi * adjc_surf_over_vol[idim];

				//c_grad
				vars->c_grad   [idim] += c    * cell_surf_over_vol[idim];
				neigh->c_grad   [idim] += c    * adjc_surf_over_vol[idim];				
			}
		}
	}
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate needed gradients for cell to face interpolation schemes
// 	- only linear interpolation
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_gradients( MPI_env &mpi_env)
{
	PRECISION cell_surf_over_vol[DIM_CNT];
	PRECISION adjc_surf_over_vol[DIM_CNT];
	PRECISION rhoU[DIM_CNT];

	for( size_t icell=cells_cfd.start(); icell < cells_cfd.end(); icell++ ){
		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		PRECISION cell_Rpsi = compute_Rpsi( vars->q_old );
		PRECISION cell_c    = sqrt( m_dGamma * cell_Rpsi );

		// Check Each Face
		for( unsigned iface=0; iface < cell->m_nFaceCount; iface++ )
		{
			// Get Neighbour
			auto *neigh = cells_cfd[icell].neighs[iface];

			// Neighbor primitive vars
			PRECISION adjc_Rpsi = compute_Rpsi( neigh->q_old );
			PRECISION adjc_c    = sqrt( m_dGamma * adjc_Rpsi );

			// Fluxes
			PRECISION rho = INTERP_LINEAR( cell->weight_linear[iface],	vars ->q_old[0],
																		neigh->q_old[0] );
			

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoU[idim] = INTERP_LINEAR( cell->weight_linear[iface], vars ->q_old[idim+1],
																		neigh->q_old[idim+1] );
			}

			PRECISION rhoE = INTERP_LINEAR( cell->weight_linear[iface], vars ->q_old[DIM_CNT+1],
																		neigh->q_old[DIM_CNT+1] );

			PRECISION Rpsi = INTERP_LINEAR( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			PRECISION c    = INTERP_LINEAR( cell->weight_linear[iface], cell_c,    adjc_c    );

			// Cell and neighbor (surface / vol );			
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =   cell->S[iface][idim] * vars ->vol_inv;
				adjc_surf_over_vol[idim] = - cell->S[iface][idim] * neigh->vol_inv;
			}

			// Cell and neighbors' gradients
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				vars->rho_grad [idim] += rho  * cell_surf_over_vol[idim];
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars->rhoU_grad[jdim][idim] += rhoU[jdim] * cell_surf_over_vol[idim];
				}
				vars->rhoE_grad[idim] += rhoE * cell_surf_over_vol[idim];
				vars->Rpsi_grad[idim] += Rpsi * cell_surf_over_vol[idim];
				vars->c_grad   [idim] += c    * cell_surf_over_vol[idim];

				neigh->rho_grad [idim] += rho  * adjc_surf_over_vol[idim];
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					neigh->rhoU_grad[jdim][idim] += rhoU[jdim] * adjc_surf_over_vol[idim];
				}
				neigh->rhoE_grad[idim] += rhoE * adjc_surf_over_vol[idim];
				neigh->Rpsi_grad[idim] += Rpsi * adjc_surf_over_vol[idim];
				neigh->c_grad   [idim] += c    * adjc_surf_over_vol[idim];
			}
		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate Derivatives for Viscosity when using Smagorinsky
// 		- interpolation: linear or midPoint
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_VIS_Smagorinsky( MPI_env &mpi_env ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Viscosity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~			
	PRECISION cell_Rpsi, adjc_Rpsi, cell_rho_inv, adjc_rho_inv;
	PRECISION face_Rpsi, face_T, mu, diagSum;
	PRECISION face_U[DIM_CNT], cell_surf_over_vol[DIM_CNT], adjc_surf_over_vol[DIM_CNT];
	PRECISION cell_U[DIM_CNT], adjc_U[DIM_CNT];
	PRECISION tau[DIM_CNT][DIM_CNT];

	// Needed for Smagorinsky model
	PRECISION smag_constant;
	PRECISION Strain_Mag,S_ij, divu;
	

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		cell_Rpsi = compute_Rpsi( vars->q_old );

		// For Smagorinsky, non dynamic, constant -2*(Cs*Delta)^2
		smag_constant = -(PRECISION)2.0 * pow((PRECISION)0.16*pow(1.0 / vars->vol_inv,1./3.),2);
		
		for( unsigned iface=0; iface < cell->m_nFaceCount; iface++ )
		{
			auto *neigh = cells_cfd[icell].neighs[iface];

			// Neighbor primitive vars
			adjc_Rpsi = compute_Rpsi( neigh->q_old );

			// Fluxes
			cell_rho_inv = ONE / vars ->q_old[0];
			adjc_rho_inv = ONE / neigh->q_old[0];
			
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				face_U[idim] = INTERP_LINEAR( cell->weight_linear[iface], vars->q_old[idim+1] * cell_rho_inv, neigh->q_old[idim+1] * adjc_rho_inv );
			}

			face_Rpsi = INTERP_LINEAR( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			face_T    = face_Rpsi * m_dRgas_inv;

			// Cell and neighbor (surface vector / vol )
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =  cell->S[iface][idim] * vars ->vol_inv ;
				adjc_surf_over_vol[idim] = -cell->S[iface][idim] * neigh->vol_inv ;
			}

			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->dudx[idim][jdim] += face_U[idim] * cell_surf_over_vol[jdim];
					neigh->dudx[idim][jdim] += face_U[idim] * adjc_surf_over_vol[jdim];
				}
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				vars ->dTdx[idim] += face_T * cell_surf_over_vol[idim] ;
				neigh->dTdx[idim] += face_T * adjc_surf_over_vol[idim] ;
			}

			mu = m_dMu0;
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				cell_U[idim] = vars->q_old[idim+1] / vars->q_old[0];
				adjc_U[idim] = neigh->q_old[idim+1] / neigh->q_old[0];
			}

			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * vars->dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (vars->dudx[nD][nD2] + vars->dudx[nD2][nD]);
					tau[nD][nD] -= vars->dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}
			
			diagSum = 0;
			for (int nD = 0; nD < DIM_CNT; nD++)
				diagSum -= vars->dudx[nD][nD];
			diagSum *= mu * 2.0 / 3.0;

			Strain_Mag = 0.0;
			S_ij = 0.0;
			divu = 0.0;

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
				vars->sigmaU[idim] = dot_product( cell_U, tau[idim] );				
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->tauMC[idim][jdim] = mu * vars ->dudx[jdim][idim];
					// Computation of strain rate magnitude
					S_ij        = HALF*(vars->dudx[idim][jdim] + vars->dudx[jdim][idim]);
					Strain_Mag += pow(S_ij,2);
				}
				vars->tauMC[idim][idim] += diagSum;
				divu += vars->dudx[idim][idim];
			}
			Strain_Mag = sqrt(2.0*Strain_Mag);

			// Contribution of the Smagorinsky model added to tauMC
			for (unsigned idim = 0; idim < DIM_CNT; idim++){
				for(unsigned jdim = 0; jdim < DIM_CNT; jdim++){
					S_ij        = HALF*(vars->dudx[idim][jdim] + vars->dudx[jdim][idim]);
					vars->tauMC[idim][jdim] -= smag_constant * Strain_Mag * S_ij;
				}
				vars->tauMC[idim][idim] += smag_constant * Strain_Mag * divu / (PRECISION)3.0; // Compressible correction
			}



			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * neigh->dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (neigh->dudx[nD][nD2] + neigh->dudx[nD2][nD]);
					tau[nD][nD] -= neigh->dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}
			
			diagSum = 0;
			for (int nD = 0; nD < DIM_CNT; nD++)
				diagSum -= neigh->dudx[nD][nD];
			diagSum *= mu * 2.0 / 3.0;

			Strain_Mag = 0.0;
			S_ij = 0.0;
			divu = 0.0;

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
				neigh->sigmaU[idim] = dot_product( adjc_U, tau[idim] );				
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					neigh ->tauMC[idim][jdim] = mu * neigh ->dudx[jdim][idim];
					// Computation of strain rate magnitude
					S_ij        = HALF*(neigh->dudx[idim][jdim] + neigh->dudx[jdim][idim]);
					Strain_Mag += pow(S_ij,2);
				}
				neigh->tauMC[idim][idim] += diagSum;
				divu += neigh->dudx[idim][idim];
			}
			Strain_Mag = sqrt(2.0*Strain_Mag);

			// Contribution of the Smagorinsky model added to tauMC
			for (unsigned idim = 0; idim < DIM_CNT; idim++){
				for(unsigned jdim = 0; jdim < DIM_CNT; jdim++){
					S_ij        = HALF*(neigh->dudx[idim][jdim] + neigh->dudx[jdim][idim]);
					neigh->tauMC[idim][jdim] -= smag_constant * Strain_Mag * S_ij;
				}
				neigh->tauMC[idim][idim] += smag_constant * Strain_Mag * divu / (PRECISION)3.0; // Compressible correction
			}

		}
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Calculate Derivatives for Viscosity
// 		- interpolation: linear or midPoint
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: calc_VIS( MPI_env &mpi_env ){

	int myrank;
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Viscosity
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
	PRECISION cell_Rpsi, adjc_Rpsi, cell_rho_inv, adjc_rho_inv;
	PRECISION face_Rpsi, face_T, mu, diagSum;
	PRECISION face_U[DIM_CNT], cell_surf_over_vol[DIM_CNT], adjc_surf_over_vol[DIM_CNT];
	PRECISION cell_U[DIM_CNT], adjc_U[DIM_CNT];
	PRECISION tau[DIM_CNT][DIM_CNT];

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		auto *cell = &cells_cfd[icell];
		auto *vars = &cells_cfd[icell].vars;

		// Cell primitive vars
		cell_Rpsi = compute_Rpsi( vars->q_old );

		for( unsigned iface=0; iface < cell->m_nFaceCount; iface++ )
		{
			auto *neigh = cells_cfd[icell].neighs[iface];

			// Neighbor primitive vars
			adjc_Rpsi = compute_Rpsi( neigh->q_old );

			// Fluxes
			cell_rho_inv = ONE / vars ->q_old[0];
			adjc_rho_inv = ONE / neigh->q_old[0];
			
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				face_U[idim] = INTERP_LINEAR( cell->weight_linear[iface], vars->q_old[idim+1] * cell_rho_inv, neigh->q_old[idim+1] * adjc_rho_inv );
			}

			face_Rpsi = INTERP_LINEAR( cell->weight_linear[iface], cell_Rpsi, adjc_Rpsi );
			face_T    = face_Rpsi * m_dRgas_inv;

			// Cell and neighbor (surface vector / vol )
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				cell_surf_over_vol[idim] =  cell->S[iface][idim] * vars ->vol_inv ;
				adjc_surf_over_vol[idim] = -cell->S[iface][idim] * neigh->vol_inv ;
			}

			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->dudx[idim][jdim] += face_U[idim] * cell_surf_over_vol[jdim];
					neigh->dudx[idim][jdim] += face_U[idim] * adjc_surf_over_vol[jdim];
				}
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				vars ->dTdx[idim] += face_T * cell_surf_over_vol[idim] ;
				neigh->dTdx[idim] += face_T * adjc_surf_over_vol[idim] ;
			}

			mu = m_dMu0;
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				cell_U[idim] = vars->q_old[idim+1] / vars->q_old[0];
				adjc_U[idim] = neigh->q_old[idim+1] / neigh->q_old[0];
			}

			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * vars->dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (vars->dudx[nD][nD2] + vars->dudx[nD2][nD]);
					tau[nD][nD] -= vars->dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}
			
			diagSum = 0;
			for (int nD = 0; nD < DIM_CNT; nD++)
				diagSum -= vars->dudx[nD][nD];
			diagSum *= mu * 2.0 / 3.0;

			for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
				vars->sigmaU[idim] = dot_product( cell_U, tau[idim] );				
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					vars ->tauMC[idim][jdim] = mu * vars ->dudx[jdim][idim];
				}
				vars->tauMC[idim][idim] += diagSum;
			}

			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * neigh->dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (neigh->dudx[nD][nD2] + neigh->dudx[nD2][nD]);
					tau[nD][nD] -= neigh->dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}

			diagSum = 0;
			for (int nD = 0; nD < DIM_CNT; nD++)
				diagSum -= neigh->dudx[nD][nD];
			diagSum *= mu * 2.0 / 3.0;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				neigh->sigmaU[idim] = dot_product( adjc_U, tau[idim] );

				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					neigh->tauMC[idim][jdim] = mu * neigh->dudx[jdim][idim];
				}
				neigh->tauMC[idim][idim] += diagSum;
			}
		}
	}
}

template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: p5Pos(PRECISION M, PRECISION alpha){
	// Computation of M2
	PRECISION M2Pos = 0.25 * (M + 1.0) * (M + 1.0);
	PRECISION M2Neg = -0.25 * (M - 1.0) * (M - 1.0);
	// Computation of M1
	PRECISION M1Pos = 0.5 * (M + abs(M));
	PRECISION M1Neg = 0.5 * (M - abs(M));
	PRECISION p5 = 0.0;
	if(abs(M) < 1){
		p5 = M2Pos * (( 2.0-M) - 16.0*alpha*M*M2Neg);
	} else {
		p5 = M1Pos / M;
	}
	return p5;
}

template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: p5Neg(PRECISION M, PRECISION alpha){
	// Computation of M2
	PRECISION M2Pos = 0.25 * (M + 1.0) * (M + 1.0);
	PRECISION M2Neg = -0.25 * (M - 1.0) * (M - 1.0);
	// Computation of M1
	PRECISION M1Pos = 0.5 * (M + abs(M));
	PRECISION M1Neg = 0.5 * (M - abs(M));
	PRECISION p5 = 0.0;
	if(abs(M) < 1){
		p5 = M2Neg * ((-2.0-M) + 16.0*alpha*M*M2Pos);
	} else {
		p5 = M1Neg / M;
	}
	return p5;
}

template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: one_rk_step_M2( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES)
{
	int myrank = mpi_env.rank();

	// Allocate data on stack
	PRECISION rhoavg, rhoavg_inv;
	PRECISION rhoUavg[DIM_CNT], uavg[DIM_CNT];
	PRECISION cell_e,adjc_e, eavg,Rpsiavg, pavg, cavg, phiavg, rhoEavg;


	PRECISION cell_rho_inv, cell_Rpsi, cell_T, adjc_rho_inv;
	PRECISION adjc_Rpsi, adjc_T, diagSum;
	PRECISION cP, cN, S_mag;
	PRECISION mu, d_mag, dmag_inv, Cp, k, delta_mag, laplacianT, divSigmaU, weight, oneMinusWeight;
	PRECISION cell_Uvec[DIM_CNT], adjc_Uvec[DIM_CNT], delta[DIM_CNT], K[DIM_CNT];
	PRECISION phiUp[DIM_CNT], d_norm[DIM_CNT], dTdx[DIM_CNT], laplacianU[DIM_CNT];
	PRECISION sigmaU[DIM_CNT], U_f[DIM_CNT];
	PRECISION rhs[DIM_CNT + 2];
	PRECISION dudx[DIM_CNT][DIM_CNT], tau[DIM_CNT][DIM_CNT];
	PRECISION divTauMC[DIM_CNT], tauMC[DIM_CNT][DIM_CNT];

	PRECISION cell_H, adjc_H, Havg;

	// M2 AUSM variables
	PRECISION cell_divu, neigh_divu, cell_curlu, neigh_curlu, cell_theta, neigh_theta,theta_avg, r,E,cell_u[DIM_CNT];
	PRECISION cell_vmag, neigh_vmag,cell_p,neigh_u[DIM_CNT],neigh_p,rhoPos,rhoNeg,pPos,pNeg;
	PRECISION uPos,uNeg,unPos,unNeg,norm,Msq,MPos,MNeg,Minf,M0,fa,alpha;
	PRECISION phalf,pu;
	
	const bool bRkCheck = (RES && rk_step == 0) ? true : false;
	
	for (size_t icell=0; icell < cells_cfd.size(); icell++)
	{
		auto *cell  = &cells_cfd[icell].vars;
		auto *param = &cells_cfd[icell];

		// Cell primitive variables
		cell_rho_inv = ONE / cell->q_old[0];
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			cell_Uvec[idim] = cell->q_old[idim+1] * cell_rho_inv ;
		}

		cell_Rpsi = compute_Rpsi( cell->q_old );
		cell_T = cell_Rpsi * m_dRgas_inv;

		// Neighbors contribution
		for(unsigned iface=0; iface < param->m_nFaceCount; iface++ ){
		
			// Fetch neighbor solution vars
			auto *neigh = param->neighs[iface];

			weight = param->weight_linear[iface];
			oneMinusWeight = ONE - weight;

			// Neighbor primitive variables
			adjc_rho_inv = ONE / neigh->q_old[0];
			
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				adjc_Uvec[idim] = neigh->q_old[idim+1] * adjc_rho_inv ;
			}

			adjc_Rpsi = compute_Rpsi( neigh->q_old );
			adjc_T = adjc_Rpsi * m_dRgas_inv;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Convection
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Following Pirozzoli's approach: use mid-point averaging for cell values
			rhoavg = HALF * (cell->q_old[0] + neigh->q_old[0]);

			rhoavg_inv = ONE / rhoavg;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoUavg[idim] = HALF * (cell->q_old[idim+1] + neigh->q_old[idim+1]);
			}

			cell_e = cell->q_old[DIM_CNT+1] * cell_rho_inv;
			adjc_e = neigh->q_old[DIM_CNT+1] * adjc_rho_inv;
			for(int nD = 0; nD < DIM_CNT; nD++){
				cell_e -= 0.5*cell_Uvec[nD]*cell_Uvec[nD];
				adjc_e -= 0.5*adjc_Uvec[nD]*adjc_Uvec[nD];
			}

			eavg = HALF * (cell_e + adjc_e);

			Rpsiavg = HALF * (cell_Rpsi + adjc_Rpsi);

			pavg = rhoavg * Rpsiavg;

			cell_H = cell->q_old[DIM_CNT+1]/cell->q_old[0] + cell_Rpsi;
            adjc_H = neigh->q_old[DIM_CNT+1]/neigh->q_old[0] + adjc_Rpsi;
            Havg = HALF * (cell_H + adjc_H);

			cP = sqrt(m_dGamma * cell_Rpsi );
			cN = sqrt(m_dGamma * adjc_Rpsi );

			cavg = HALF * (cP + cN);

			phiavg = 0.0;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				uavg[idim] = rhoUavg[idim] * rhoavg_inv ;
				phiavg +=   uavg[idim] * param->S[iface][idim];
			}

			S_mag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				S_mag += pow( cells_cfd[icell].S[iface][idim], 2 );
			}
			S_mag = sqrt( S_mag );

			// Contribution of the Eulerian fluxes to RHS

			rhoEavg = eavg;
			for(int nD = 0; nD < DIM_CNT; nD++){
				rhoEavg += 0.5 * uavg[nD]*uavg[nD];
			}
			rhoEavg *= rhoavg;

			rhs[0] = - rhoavg * phiavg;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] = - (rhoUavg[idim] * phiavg + pavg * param->S[iface][idim]);
			}
			rhs[DIM_CNT+1] = -(rhoavg*Havg*phiavg);//- (rhoEavg * phiavg);
			
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Viscosity
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Derivatives
			mu = m_dMu0;

			d_mag = vector_mag( param->d[iface] );
			dmag_inv = ONE / d_mag ;
			
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				d_norm[idim] = param->d[iface][idim] * dmag_inv ;
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					dudx[idim][jdim] = INTERP_LINEAR( weight, cell->dudx[idim][jdim], neigh->dudx[idim][jdim] );
				}
				dTdx[idim] = INTERP_LINEAR( weight, cell->dTdx[idim], neigh->dTdx[idim] );
			}

			// explicit correction for physical boundaries (wall/inlet/outlet)		
			if (neigh->is_ghost)
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
						dudx[idim][jdim] = (adjc_Uvec[idim] - cell_Uvec[idim]) * d_norm[jdim] * dmag_inv;
					}
					dTdx[idim] = (adjc_T - cell_T) * d_norm[idim] * dmag_inv;
				}
			
			// Stress tensor
			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (dudx[nD][nD2] + dudx[nD2][nD]);
					tau[nD][nD] -= dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}

			// Head conductivity			
			Cp = m_dCp;
			k  = Cp * mu * m_dPr_inv;

			// div(tauMC)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					tauMC[idim][jdim] = INTERP_LINEAR( weight, cell->tauMC[idim][jdim], neigh->tauMC[idim][jdim] );
				}
			}

			if (neigh->is_ghost)
			{
				diagSum = 0;
				for (int nD = 0; nD < DIM_CNT; nD++)
					diagSum -= dudx[nD][nD];
				diagSum *= mu * 2.0 / 3.0;

				for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
					for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
						tauMC[idim][jdim] = mu * dudx[jdim][idim];
					}
					tauMC[idim][idim] += diagSum;
				}
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				divTauMC[idim] = dot_product( tauMC[idim], param->S[iface] );
			}

			// Laplacians
			delta_mag = 0.0;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				delta[idim] = param->d[iface][idim] * S_mag * S_mag / dot_product( param->S[iface], param->d[iface] );
				delta_mag  += delta[idim] * delta[idim];
				K[idim]     = param->S[iface][idim] - delta[idim];
			}
			delta_mag = sqrt(delta_mag);			

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				laplacianU[idim] = mu * ( delta_mag * (adjc_Uvec[idim] - cell_Uvec[idim]) * dmag_inv + dot_product( K, dudx[idim] ) );
			}

			laplacianT = k * ( delta_mag * (adjc_T - cell_T) * dmag_inv + dot_product( K, dTdx ) );

			// div(sigmaU)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sigmaU[idim] = INTERP_LINEAR( weight, cell->sigmaU[idim], neigh->sigmaU[idim] );
			}

			if (neigh->is_ghost)
			{
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					U_f[idim] = INTERP_LINEAR( weight, cell_Uvec[idim], adjc_Uvec[idim] );
				}
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					sigmaU[idim] = dot_product( U_f, tau[idim] );
				}
			}

			divSigmaU = dot_product( sigmaU, param->S[iface] );

			// RHS contribution
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] += divTauMC[idim] + laplacianU[idim];
			}
			rhs[DIM_CNT+1] += divSigmaU + laplacianT;

			if (bRkCheck)
			{
			 	for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){			 		
			 		cell ->RES[idim] += rhs[idim];
			 		neigh->RES[idim] -= rhs[idim];
			 	}
			}

			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
			{
				cell ->delta_q[idim] += dt * rhs[idim] * cell ->vol_inv;
				neigh->delta_q[idim] -= dt * rhs[idim] * neigh->vol_inv;
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Sponge Layer contribution
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cell->delta_q[0] += dt * param->sponge_sigma * ( m_drhoInf - cell->q_old[0] );
		for( unsigned idim=0; idim < DIM_CNT; idim++ ){
			cell->delta_q[idim+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dUInf[idim] - cell->q_old[idim+1] );
		}
		cell->delta_q[DIM_CNT+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dEInf - cell->q_old[DIM_CNT+1] );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// 6 - Compute next RK step
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
			cell->q_old[idim] += m_dBk[rk_step] * cell->delta_q[idim];
		}
	}

	if (bRkCheck)
		for (size_t icell=0; icell < cells_cfd.size(); icell++)
		{
			auto *cell  = &cells_cfd[icell].vars;
			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
				RES[idim] += cell->RES[idim] * cell->RES[idim];
		}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Perform numerical scheme for the Hpath and regular cells
// 	- interpolation: linear (convection) or midPoint (viscous)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: one_rk_step_M2AUSM( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES)
{
	int myrank = mpi_env.rank();

	// Allocate data on stack
	PRECISION rhoavg, rhoavg_inv;
	PRECISION rhoUavg[DIM_CNT], uavg[DIM_CNT];
	PRECISION cell_e,adjc_e, eavg,Rpsiavg, pavg, cavg, phiavg, rhoEavg;


	PRECISION cell_rho_inv, cell_Rpsi, cell_T, adjc_rho_inv;
	PRECISION adjc_Rpsi, adjc_T, diagSum;
	PRECISION cP, cN, S_mag;
	PRECISION mu, d_mag, dmag_inv, Cp, k, delta_mag, laplacianT, divSigmaU, weight, oneMinusWeight;
	PRECISION cell_Uvec[DIM_CNT], adjc_Uvec[DIM_CNT], delta[DIM_CNT], K[DIM_CNT];
	PRECISION phiUp[DIM_CNT], d_norm[DIM_CNT], dTdx[DIM_CNT], laplacianU[DIM_CNT];
	PRECISION sigmaU[DIM_CNT], U_f[DIM_CNT];
	PRECISION rhs[DIM_CNT + 2];
	PRECISION dudx[DIM_CNT][DIM_CNT], tau[DIM_CNT][DIM_CNT];
	PRECISION divTauMC[DIM_CNT], tauMC[DIM_CNT][DIM_CNT];
	PRECISION cell_H, adjc_H, Havg;

	// M2 AUSM variables
	PRECISION cell_divu, neigh_divu, cell_curlu, neigh_curlu, cell_theta, neigh_theta,theta_avg, r,E,cell_u[DIM_CNT];
	PRECISION cell_vmag, neigh_vmag,cell_p,neigh_u[DIM_CNT],neigh_p,rhoPos,rhoNeg,pPos,pNeg;
	PRECISION uPos,uNeg,unPos,unNeg,norm,Msq,MPos,MNeg,Minf,M0,fa,alpha;
	PRECISION phalf,pu;
	
	const bool bRkCheck = (RES && rk_step == 0) ? true : false;
	
	for (size_t icell=0; icell < cells_cfd.size(); icell++)
	{
		auto *cell  = &cells_cfd[icell].vars;
		auto *param = &cells_cfd[icell];

		// Cell primitive variables
		cell_rho_inv = ONE / cell->q_old[0];
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			cell_Uvec[idim] = cell->q_old[idim+1] * cell_rho_inv ;
		}

		cell_Rpsi = compute_Rpsi( cell->q_old );
		cell_T = cell_Rpsi * m_dRgas_inv;

		// Neighbors contribution
		for(unsigned iface=0; iface < param->m_nFaceCount; iface++ ){
		
			// Fetch neighbor solution vars
			auto *neigh = param->neighs[iface];

			weight = param->weight_linear[iface];
			oneMinusWeight = ONE - weight;

			// Neighbor primitive variables
			adjc_rho_inv = ONE / neigh->q_old[0];
			
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				adjc_Uvec[idim] = neigh->q_old[idim+1] * adjc_rho_inv ;
			}

			adjc_Rpsi = compute_Rpsi( neigh->q_old );
			adjc_T = adjc_Rpsi * m_dRgas_inv;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Convection
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

			// Following Pirozzoli's approach: use mid-point averaging for cell values
			rhoavg = HALF * (cell->q_old[0] + neigh->q_old[0]);

			rhoavg_inv = ONE / rhoavg;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoUavg[idim] = HALF * (cell->q_old[idim+1] + neigh->q_old[idim+1]);
			}

			cell_e = cell->q_old[DIM_CNT+1] * cell_rho_inv;
			adjc_e = neigh->q_old[DIM_CNT+1] * adjc_rho_inv;
			for(int nD = 0; nD < DIM_CNT; nD++){
				cell_e -= 0.5*cell_Uvec[nD]*cell_Uvec[nD];
				adjc_e -= 0.5*adjc_Uvec[nD]*adjc_Uvec[nD];
			}

			eavg = HALF * (cell_e + adjc_e);

			Rpsiavg = HALF * (cell_Rpsi + adjc_Rpsi);

			pavg = rhoavg * Rpsiavg;
			
			cell_H = cell->q_old[DIM_CNT+1]/cell->q_old[0] + cell_Rpsi;
            adjc_H = neigh->q_old[DIM_CNT+1]/neigh->q_old[0] + adjc_Rpsi;
            Havg = HALF * (cell_H + adjc_H);

			cP = sqrt(m_dGamma * cell_Rpsi );
			cN = sqrt(m_dGamma * adjc_Rpsi );

			cavg = HALF * (cP + cN);

			cell_divu = 0.; neigh_divu = 0.; cell_curlu = 0.; neigh_curlu = 0.; cell_vmag = 0.; neigh_vmag = 0.;
			unPos = 0.; unNeg = 0.; norm = 0.;

			for(int nD = 0; nD < DIM_CNT; nD++){
				cell_divu  +=  cell->dudx[nD][nD];
				neigh_divu += neigh->dudx[nD][nD];
				for(int nD2 = nD+1; nD2 < DIM_CNT; nD2++){
					cell_curlu  += (cell->dudx[nD][nD2] - cell->dudx[nD2][nD]) * (cell->dudx[nD][nD2] - cell->dudx[nD2][nD]);
					neigh_curlu += (neigh->dudx[nD][nD2] - neigh->dudx[nD2][nD]) * (neigh->dudx[nD][nD2] - neigh->dudx[nD2][nD]);
				}
				cell_u[nD] = cell->q_old[nD+1] / cell->q_old[0];
				cell_vmag += cell_u[nD] * cell_u[nD];
				neigh_u[nD] = neigh->q_old[nD+1] / neigh->q_old[0];
				neigh_vmag += neigh_u[nD] * neigh_u[nD];
				// MinMod of normal velocity
				uPos   = interp_minmod( cell_u[nD], neigh_u[nD], cell ->U_grad[nD], param->d[iface], param->weight_linear[iface],   ONE );
				uNeg   = interp_minmod( cell_u[nD], neigh_u[nD], neigh->U_grad[nD], param->d[iface], param->weight_linear[iface], M_ONE );
				unPos += uPos*param->S[iface][nD];
				unNeg += uNeg*param->S[iface][nD];
				// norm
				norm += param->S[iface][nD] * param->S[iface][nD];
			}
			// Sensor computation
			cell_theta = 	max(-(cell_divu/sqrt(cell_divu*cell_divu+cell_curlu+4e-2)),0.);
			neigh_theta = 	max(-(neigh_divu/sqrt(neigh_divu*neigh_divu+neigh_curlu+4e-2)),0.);
			theta_avg = HALF * (cell_theta + neigh_theta);
			// Local pressure
			r = cell->q_old[0];
			E = cell->q_old[DIM_CNT+1] / r;
			cell_p = r * m_dGammaMinusOne * (E-HALF*cell_vmag); // Hard coded for now for quick debug and testing
			r = neigh->q_old[0];
			E = neigh->q_old[DIM_CNT+1] / r;
			neigh_p = r * m_dGammaMinusOne * (E-HALF*neigh_vmag); // Hard coded for now for quick debug and testing
			// Minmod interpolation of density
			rhoPos = interp_minmod( cell->q_old[0], neigh->q_old[0], cell ->rho_grad, param->d[iface], param->weight_linear[iface],   ONE );
			rhoNeg = interp_minmod( cell->q_old[0], neigh->q_old[0], neigh->rho_grad, param->d[iface], param->weight_linear[iface], M_ONE );
			// Minmod interpolation of pressure
			pPos   = interp_minmod( cell_p, neigh_p, cell ->p_grad, param->d[iface], param->weight_linear[iface],   ONE );
			pNeg   = interp_minmod( cell_p, neigh_p, neigh->p_grad, param->d[iface], param->weight_linear[iface], M_ONE );
			// Computation of local Mach number
			Msq  = (unPos*unPos + unNeg*unNeg) / (2.*cavg*cavg*norm);
			MPos = unPos/(sqrt(norm)*cavg);
			MNeg = unNeg/(sqrt(norm)*cavg);
			Minf = 0.2; // Hard coded for now for quick debug and testing
			M0   = sqrt(min((PRECISION)1.0,max(Msq,Minf*Minf)));
			// Computation of fa
			fa = M0 * (2.0 - M0);
			alpha = 3.0 * (5.0*fa*fa - 4.0) / 16.0;
			// Computation of phalf
			phalf = pNeg*(p5Pos(MNeg,alpha)-p5Neg(MNeg,alpha)) - pPos*(p5Pos(MPos,alpha)-p5Neg(MPos,alpha));
			// Computation of pu
			pu = -0.75*p5Pos(MPos,alpha)*p5Neg(MNeg,alpha)*(rhoPos+rhoNeg)*cavg*fa*(unNeg-unPos);
			// Computation of pAUSM
			pavg += theta_avg*(pu - HALF*phalf);			

			phiavg = 0.0;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				uavg[idim] = rhoUavg[idim] * rhoavg_inv ;
				phiavg +=   uavg[idim] * param->S[iface][idim];
			}

			S_mag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				S_mag += pow( cells_cfd[icell].S[iface][idim], 2 );
			}
			S_mag = sqrt( S_mag );

			// Contribution of the Eulerian fluxes to RHS

			rhoEavg = eavg;
			for(int nD = 0; nD < DIM_CNT; nD++){
				rhoEavg += 0.5 * uavg[nD]*uavg[nD];
			}
			rhoEavg *= rhoavg;

			rhs[0] = - rhoavg * phiavg;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] = - (rhoUavg[idim] * phiavg + pavg * param->S[iface][idim]);
			}
			rhs[DIM_CNT+1] = -(rhoavg*Havg*phiavg);//- (rhoEavg * phiavg);
			
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Viscosity
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Derivatives			
			mu = m_dMu0;

			d_mag = vector_mag( param->d[iface] );
			dmag_inv = ONE / d_mag ;
			
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				d_norm[idim] = param->d[iface][idim] * dmag_inv ;
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					dudx[idim][jdim] = INTERP_LINEAR( weight, cell->dudx[idim][jdim], neigh->dudx[idim][jdim] );
				}
				dTdx[idim] = INTERP_LINEAR( weight, cell->dTdx[idim], neigh->dTdx[idim] );
			}

			// explicit correction for physical boundaries (wall/inlet/outlet)		
			if (neigh->is_ghost)
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
						dudx[idim][jdim] = (adjc_Uvec[idim] - cell_Uvec[idim]) * d_norm[jdim] * dmag_inv;
					}
					dTdx[idim] = (adjc_T - cell_T) * d_norm[idim] * dmag_inv;
				}
			
			// Stress tensor
			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (dudx[nD][nD2] + dudx[nD2][nD]);
					tau[nD][nD] -= dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}

			// Head conductivity			
			Cp = m_dCp;
			k  = Cp * mu * m_dPr_inv;

			// div(tauMC)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					tauMC[idim][jdim] = INTERP_LINEAR( weight, cell->tauMC[idim][jdim], neigh->tauMC[idim][jdim] );
				}
			}

			if (neigh->is_ghost)
			{
				diagSum = 0;
				for (int nD = 0; nD < DIM_CNT; nD++)
					diagSum -= dudx[nD][nD];
				diagSum *= mu * 2.0 / 3.0;

				for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
					for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
						tauMC[idim][jdim] = mu * dudx[jdim][idim];
					}
					tauMC[idim][idim] += diagSum;
				}
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				divTauMC[idim] = dot_product( tauMC[idim], param->S[iface] );
			}

			// Laplacians
			delta_mag = 0.0;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				delta[idim] = param->d[iface][idim] * S_mag * S_mag / dot_product( param->S[iface], param->d[iface] );
				delta_mag  += delta[idim] * delta[idim];
				K[idim]     = param->S[iface][idim] - delta[idim];
			}
			delta_mag = sqrt(delta_mag);			

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				laplacianU[idim] = mu * ( delta_mag * (adjc_Uvec[idim] - cell_Uvec[idim]) * dmag_inv + dot_product( K, dudx[idim] ) );
			}

			laplacianT = k * ( delta_mag * (adjc_T - cell_T) * dmag_inv + dot_product( K, dTdx ) );

			// div(sigmaU)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sigmaU[idim] = INTERP_LINEAR( weight, cell->sigmaU[idim], neigh->sigmaU[idim] );
			}

			if (neigh->is_ghost)
			{
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					U_f[idim] = INTERP_LINEAR( weight, cell_Uvec[idim], adjc_Uvec[idim] );
				}
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					sigmaU[idim] = dot_product( U_f, tau[idim] );
				}
			}

			divSigmaU = dot_product( sigmaU, param->S[iface] );

			// RHS contribution
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] += divTauMC[idim] + laplacianU[idim];
			}
			rhs[DIM_CNT+1] += divSigmaU + laplacianT;

			if (bRkCheck)
			{
			 	for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){			 		
			 		cell ->RES[idim] += rhs[idim];
			 		neigh->RES[idim] -= rhs[idim];
			 	}
			}

			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
			{
				cell ->delta_q[idim] += dt * rhs[idim] * cell ->vol_inv;
				neigh->delta_q[idim] -= dt * rhs[idim] * neigh->vol_inv;
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Sponge Layer contribution
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cell->delta_q[0] += dt * param->sponge_sigma * ( m_drhoInf - cell->q_old[0] );
		for( unsigned idim=0; idim < DIM_CNT; idim++ ){
			cell->delta_q[idim+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dUInf[idim] - cell->q_old[idim+1] );
		}
		cell->delta_q[DIM_CNT+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dEInf - cell->q_old[DIM_CNT+1] );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// 6 - Compute next RK step
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
			cell->q_old[idim] += m_dBk[rk_step] * cell->delta_q[idim];
		}

	}

	if (bRkCheck)
		for (size_t icell=0; icell < cells_cfd.size(); icell++)
		{
			auto *cell  = &cells_cfd[icell].vars;
			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
				RES[idim] += cell->RES[idim] * cell->RES[idim];
		}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Perform numerical scheme for the Hpath and regular cells
// 	- interpolation: linear (convection) or midPoint (viscous)
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: one_rk_step_M1( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES )
{
	int myrank = mpi_env.rank();

	// Allocate data on stack
	PRECISION cell_rho_inv, cell_Rpsi, cell_T, adjc_rho_inv;
	PRECISION adjc_Rpsi, adjc_T, rhoPos, rhoNeg, rhoPos_inv, rhoNeg_inv;
	PRECISION cell_e, adjc_e, ePos, eNeg, RpsiPos, RpsiNeg;
	PRECISION pPos, pNeg, cP, cN, cPos, cNeg, phiPos = 0.0, phiNeg = 0.0;
	PRECISION S_mag, psiPos, psiNeg, a0, a1, aPos, aNeg, phi, rhoEPos, rhoENeg, diagSum;
	PRECISION mu, d_mag, dmag_inv, Cp, k, delta_mag, laplacianT, divSigmaU, weight, oneMinusWeight;
	PRECISION cell_Uvec[DIM_CNT], adjc_Uvec[DIM_CNT], delta[DIM_CNT], K[DIM_CNT];
	PRECISION rhoUPos[DIM_CNT], rhoUNeg[DIM_CNT], uPos[DIM_CNT], uNeg[DIM_CNT];
	PRECISION phiUp[DIM_CNT], d_norm[DIM_CNT], dTdx[DIM_CNT], laplacianU[DIM_CNT];
	PRECISION sigmaU[DIM_CNT], U_f[DIM_CNT];
	PRECISION rhs[DIM_CNT + 2];
	PRECISION dudx[DIM_CNT][DIM_CNT], tau[DIM_CNT][DIM_CNT];
	PRECISION divTauMC[DIM_CNT], tauMC[DIM_CNT][DIM_CNT];

	const bool bRkCheck = (RES && rk_step == 0) ? true : false;
	
	for (size_t icell=0; icell < cells_cfd.size(); icell++)
	{
		auto *cell  = &cells_cfd[icell].vars;
		auto *param = &cells_cfd[icell];

		// Cell primitive variables
		cell_rho_inv = ONE / cell->q_old[0];
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
			cell_Uvec[idim] = cell->q_old[idim+1] * cell_rho_inv ;
		}

		cell_Rpsi = compute_Rpsi( cell->q_old );
		cell_T = cell_Rpsi * m_dRgas_inv;

		// Neighbors contribution
		for(unsigned iface=0; iface < param->m_nFaceCount; iface++ ){
		
			// Fetch neighbor solution vars
			auto *neigh = param->neighs[iface];

			weight = param->weight_linear[iface];
			oneMinusWeight = ONE - weight;

			// Neighbor primitive variables
			adjc_rho_inv = ONE / neigh->q_old[0];
			
			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				adjc_Uvec[idim] = neigh->q_old[idim+1] * adjc_rho_inv ;
			}

			adjc_Rpsi = compute_Rpsi( neigh->q_old );
			adjc_T = adjc_Rpsi * m_dRgas_inv;

			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Convection
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Fluxes
			rhoPos = INTERP_LINEAR(       weight, cell ->q_old[0], neigh->q_old[0] );
			rhoNeg = INTERP_LINEAR( oneMinusWeight, neigh->q_old[0], cell ->q_old[0] );

			rhoPos_inv = ONE / rhoPos;
			rhoNeg_inv = ONE / rhoNeg;


			#pragma omp simd
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoUPos[idim] = INTERP_LINEAR(       weight, cell ->q_old[idim+1], neigh->q_old[idim+1] );
				rhoUNeg[idim] = INTERP_LINEAR( oneMinusWeight, neigh->q_old[idim+1], cell ->q_old[idim+1] );
			}

			cell_e = 2 * cell ->q_old[DIM_CNT + 1] * cell_rho_inv;
			adjc_e = 2 * neigh->q_old[DIM_CNT + 1] * adjc_rho_inv;
			for (int nD = 0; nD < DIM_CNT; nD++)
			{
				cell_e -= cell_Uvec[nD] * cell_Uvec[nD];
				adjc_e -= adjc_Uvec[nD] * adjc_Uvec[nD];
			}
			cell_e *= 0.5;
			adjc_e *= 0.5;
			
			ePos = INTERP_LINEAR(       weight, cell_e, adjc_e );
			eNeg = INTERP_LINEAR( oneMinusWeight, adjc_e, cell_e );

			RpsiPos = INTERP_LINEAR(       weight, cell_Rpsi, adjc_Rpsi );
			RpsiNeg = INTERP_LINEAR( oneMinusWeight, adjc_Rpsi, cell_Rpsi );

			pPos = rhoPos * RpsiPos;
			pNeg = rhoNeg * RpsiNeg;

			cP = sqrt(m_dGamma * cell_Rpsi );
			cN = sqrt(m_dGamma * adjc_Rpsi );

			cPos = INTERP_LINEAR(       weight, cP, cN );
			cNeg = INTERP_LINEAR( oneMinusWeight, cN, cP );

			phiPos = 0.0;
			phiNeg = 0.0;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				uPos[idim] = rhoUPos[idim] * rhoPos_inv ;
				uNeg[idim] = rhoUNeg[idim] * rhoNeg_inv ;

				phiPos +=   uPos[idim] * param->S[iface][idim];
				phiNeg += - uNeg[idim] * param->S[iface][idim];
			}

			S_mag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				S_mag += pow( cells_cfd[icell].S[iface][idim], 2 );
			}
			S_mag = sqrt( S_mag );

			psiPos = max( { phiPos + cPos * S_mag, - phiNeg + cNeg * S_mag, ZERO } );
			psiNeg = min( { phiPos - cPos * S_mag, - phiNeg - cNeg * S_mag, ZERO } );

			// RHS contribution
			a0 = ONE / (psiPos - psiNeg );
			a1 = psiPos * psiNeg ;
			aPos = psiPos * phiPos ;
			aNeg = psiNeg * phiNeg ;
			phi = (aPos*rhoPos + aNeg*rhoNeg + (rhoNeg - rhoPos) * a1 ) * a0;

			rhoEPos = ePos;
			rhoENeg = eNeg;
			for (int nD = 0; nD < DIM_CNT; nD++)
			{
				rhoEPos += 0.5 * uPos[nD] * uPos[nD];
				rhoENeg += 0.5 * uNeg[nD] * uNeg[nD];
			}
			rhoEPos *= rhoPos;
			rhoENeg *= rhoNeg;

			
			rhs[0] = - (aPos*rhoPos + aNeg*rhoNeg + (rhoNeg - rhoPos) * a1 ) * a0;

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ )
			{
				phiUp[idim] = (aPos*rhoUPos[idim] + aNeg*rhoUNeg[idim] + (rhoUNeg[idim] - rhoUPos[idim]) * a1 ) * a0 + (pPos * psiPos - pNeg * psiNeg) * a0 * param->S[iface][idim];
				rhs[idim+1] = - phiUp[idim];
			}

			rhs[DIM_CNT+1] = - ( aPos*rhoEPos + aNeg*rhoENeg + (rhoENeg - rhoEPos) * a1 + (aPos*pPos + aNeg*pNeg) ) * a0;


			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Viscosity
			//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
			// Derivatives
			mu = m_dMu0;

			d_mag = vector_mag( param->d[iface] );
			dmag_inv = ONE / d_mag ;
			
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				d_norm[idim] = param->d[iface][idim] * dmag_inv ;
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					dudx[idim][jdim] = INTERP_LINEAR( weight, cell->dudx[idim][jdim], neigh->dudx[idim][jdim] );
				}
				dTdx[idim] = INTERP_LINEAR( weight, cell->dTdx[idim], neigh->dTdx[idim] );
			}

			// explicit correction for physical boundaries (wall/inlet/outlet)		
			if (neigh->is_ghost)
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
						dudx[idim][jdim] = (adjc_Uvec[idim] - cell_Uvec[idim]) * d_norm[jdim] * dmag_inv;
					}
					dTdx[idim] = (adjc_T - cell_T) * d_norm[idim] * dmag_inv;
				}
			
			// Stress tensor
			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (dudx[nD][nD2] + dudx[nD2][nD]);
					tau[nD][nD] -= dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}

			// Head conductivity			
			Cp = m_dCp;
			k  = Cp * mu * m_dPr_inv;

			// div(tauMC)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					tauMC[idim][jdim] = INTERP_LINEAR( weight, cell->tauMC[idim][jdim], neigh->tauMC[idim][jdim] );
				}
			}

			if (neigh->is_ghost)
			{
				diagSum = 0;
				for (int nD = 0; nD < DIM_CNT; nD++)
					diagSum -= dudx[nD][nD];
				diagSum *= mu * 2.0 / 3.0;

				for( unsigned idim=0; idim < DIM_CNT; idim++ ){				
					for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
						tauMC[idim][jdim] = mu * dudx[jdim][idim];
					}
					tauMC[idim][idim] += diagSum;
				}
			}

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				divTauMC[idim] = dot_product( tauMC[idim], param->S[iface] );
			}

			// Laplacians
			delta_mag = 0.0;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				delta[idim] = param->d[iface][idim] * S_mag * S_mag / dot_product( param->S[iface], param->d[iface] );
				delta_mag  += delta[idim] * delta[idim];
				K[idim]     = param->S[iface][idim] - delta[idim];
			}
			delta_mag = sqrt(delta_mag);			

			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				laplacianU[idim] = mu * ( delta_mag * (adjc_Uvec[idim] - cell_Uvec[idim]) * dmag_inv + dot_product( K, dudx[idim] ) );
			}

			laplacianT = k * ( delta_mag * (adjc_T - cell_T) * dmag_inv + dot_product( K, dTdx ) );

			// div(sigmaU)
			#pragma omp simd
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sigmaU[idim] = INTERP_LINEAR( weight, cell->sigmaU[idim], neigh->sigmaU[idim] );
			}

			if (neigh->is_ghost)
			{
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					U_f[idim] = INTERP_LINEAR( weight, cell_Uvec[idim], adjc_Uvec[idim] );
				}
				for( unsigned idim=0; idim < DIM_CNT; idim++ ){
					sigmaU[idim] = dot_product( U_f, tau[idim] );
				}
			}

			divSigmaU = dot_product( sigmaU, param->S[iface] );

			// RHS contribution
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				rhs[idim+1] += divTauMC[idim] + laplacianU[idim];
			}
			rhs[DIM_CNT+1] += divSigmaU + laplacianT;

			if (bRkCheck)
			{
			 	for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){			 		
			 		cell ->RES[idim] += rhs[idim];
			 		neigh->RES[idim] -= rhs[idim];
			 	}
			}

			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
			{
				cell ->delta_q[idim] += dt * rhs[idim] * cell ->vol_inv;
				neigh->delta_q[idim] -= dt * rhs[idim] * neigh->vol_inv;
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Sponge Layer contribution
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		cell->delta_q[0] += dt * param->sponge_sigma * ( m_drhoInf - cell->q_old[0] );
		for( unsigned idim=0; idim < DIM_CNT; idim++ ){
			cell->delta_q[idim+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dUInf[idim] - cell->q_old[idim+1] );
		}
		cell->delta_q[DIM_CNT+1] += dt * param->sponge_sigma * ( m_drhoInf * m_dEInf - cell->q_old[DIM_CNT+1] );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// 6 - Compute next RK step
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		#pragma omp simd
		for( unsigned idim = 0; idim < DIM_CNT+2; idim++ ){
			cell->q_old[idim] += m_dBk[rk_step] * cell->delta_q[idim];
		}
	}

	if (bRkCheck)
		for (size_t icell=0; icell < cells_cfd.size(); icell++)
		{
			auto *cell  = &cells_cfd[icell].vars;
			for( unsigned idim = 0; idim < DIM_CNT+2; idim++ )
				RES[idim] += cell->RES[idim] * cell->RES[idim];
		}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Exchange ghosts cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: exchange_ghost_cells( MPI_env &mpi_env ){

	// Serial case
	if( mpi_env.is_serial() ){		
		if ( !ghost_mpi.empty() && !ghost_mpi[0].empty() ) {
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				ghost_mpi[0][icell].vars = cells_cfd[icell].vars;
			}
		}
		return;
	}

	// Depending on framework used...
	switch( mpi_env.get_halo_comm_type() ){
		case MPI_ONESIDED_NONB:
			// mpi_env.barrier();
			mpi_env.irget( ghost_mpi );
			break;
		case MPI_TWOSIDED_NONB:
			mpi_env.isendrecv( cells_cfd, ghost_mpi );
			break;
		case MPI_TWOSIDED_BLCK:
			mpi_env.sendrecv( cells_cfd, ghost_mpi );
			break;
		case MPI_ONESIDED_BLCK:
			// mpi_env.barrier();
			mpi_env.rget( ghost_mpi );
			break;
		case MPI_TWOSIDED_PERS:
			mpi_env.startAllPersistantSendRecv(0);
			break;
		case MPI_ONESIDED_PERS:
			break;
		case MPI_NEIGHCOLL_BLCK:			
			break;
		case MPI_NEIGHCOLL_NONB:			
			break;
		case MPI_NEIGHCOLL_PERS:			
			break;
		case MPI_NONE:
			cout << "MPI halo communication type has not been set. Terminate simulation." << endl;
			MPI_Abort( MPI_COMM_WORLD, 910 );
		default:
			break;
	}
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: compute_cfl( PRECISION dt ){

	PRECISION cfl, cfl_max = 0.;
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		cfl = 0;
		PRECISION rho_inv = ONE / cells_cfd[icell].vars.q_old[0];

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			PRECISION flux = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				flux += cells_cfd[icell].vars.q_old[idim+1] * cells_cfd[icell].S[iface][idim];
			}
			cfl += abs(flux);
		}
		cfl *= rho_inv * cells_cfd[icell].vars.vol_inv * dt;

		if( cfl_max < cfl )
			cfl_max = cfl;
	}

	return cfl_max;
}

//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
PRECISION CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: compute_dt( PRECISION cflMax ){

	PRECISION sumFlux, dt, dt_min=100.0;
	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		sumFlux = 0;
		PRECISION rho = cells_cfd[icell].vars.q_old[0];

		for( unsigned iface=0; iface < FACE_CNT; iface++ ){
			PRECISION sMag = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				sMag += cells_cfd[icell].S[iface][idim] * cells_cfd[icell].S[iface][idim];
			}
			sMag = sqrt( sMag );

			PRECISION flux = 0.;
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				flux += cells_cfd[icell].vars.q_old[idim+1] * cells_cfd[icell].S[iface][idim] / sMag;
			}
			sumFlux += abs(flux);
		}
		PRECISION cells_cfd_vol = 1.0 / cells_cfd[icell].vars.vol_inv;
		dt = cflMax * rho * 2.0 * cells_cfd_vol / sumFlux;

		if( dt_min > dt )
			dt_min = dt;
	}

	return dt_min;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: gather_info(
									MPI_env									 &mpi_env,
									struct t_domain_info<PRECISION, DIM_CNT> &domain_info ){

	domain_info.rho_min = cells_cfd[0].vars.q_old[0];
	domain_info.rho_max = cells_cfd[0].vars.q_old[0];

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){

		if( domain_info.rho_min > cells_cfd[icell].vars.q_old[0] )
			domain_info.rho_min = cells_cfd[icell].vars.q_old[0];

		if( domain_info.rho_max < cells_cfd[icell].vars.q_old[0] )
			domain_info.rho_max = cells_cfd[icell].vars.q_old[0];
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::updateAverageField(const int nTimeStep)
{
	m_pMeshReader->updateScalarField(m_pMeshReader->getScalarFieldIndex("pAvg"), m_dpAVG, m_nSubmeshIndex, 1.0 / nTimeStep);
	m_pMeshReader->updateScalarFieldSqrt(m_pMeshReader->getScalarFieldIndex("pRMS"), m_dpRMS, m_nSubmeshIndex, 1.0 / nTimeStep);
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::updateSolutionResidual()
{
	// Field Index
	const int nRhoIndex = m_pMeshReader->getScalarFieldIndex("rhoRES");
	const int nUVWIndex[3] = {m_pMeshReader->getScalarFieldIndex("uRES"),
							  m_pMeshReader->getScalarFieldIndex("vRES"),
							  m_pMeshReader->getScalarFieldIndex("wRES")};
	const int nEIndex = m_pMeshReader->getScalarFieldIndex("eRES");
	
	// Update Fields
	for (int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
	{
		// Cell Index
		const int nCellIndex = cells_cfd[nCIndex].m_nCellIndex;
		
		// Update Fields
		m_pMeshReader->updateScalarField(nRhoIndex, cells_cfd[nCIndex].vars.RES[0], nCellIndex);
		for (int nD = 0; nD < DIM_CNT; nD++)
			m_pMeshReader->updateScalarField(nUVWIndex[nD], cells_cfd[nCIndex].vars.RES[nD + 1], nCellIndex);
		m_pMeshReader->updateScalarField(nEIndex, cells_cfd[nCIndex].vars.RES[DIM_CNT + 1], nCellIndex);
	}
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::updateSolutionPrimitives()
{
	// Field Index
	const int nRhoIndex = m_pMeshReader->getScalarFieldIndex("rho");
	const int nEIndex = m_pMeshReader->getScalarFieldIndex("E");
	const int nPIndex = m_pMeshReader->getScalarFieldIndex("p");
	const int nUIndex = m_pMeshReader->getVectorFieldIndex("U");
	
	// Update Fields
	PRECISION r, Umag2, E, p;
	PRECISION uvw[3] = {0};
	for (int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
	{
		// Cell Index
		const int nCellIndex = cells_cfd[nCIndex].m_nCellIndex;
		
		// Calculate primitives
		r = cells_cfd[nCIndex].vars.q_old[0];
		Umag2 = 0;
		for (int nD = 0; nD < DIM_CNT; nD++)
		{
			uvw[nD] = cells_cfd[nCIndex].vars.q_old[nD + 1] / r;
			Umag2 += uvw[nD] * uvw[nD];
		}
		E = cells_cfd[nCIndex].vars.q_old[DIM_CNT + 1] / r;		
		p = r * m_dGammaMinusOne * (E - 0.5 * Umag2);

		// Update Fields
		m_pMeshReader->updateScalarField(nRhoIndex, r, nCellIndex);
		m_pMeshReader->updateScalarField(nEIndex, E, nCellIndex);
		m_pMeshReader->updateScalarField(nPIndex, p, nCellIndex);
		m_pMeshReader->updateVectorField(nUIndex, uvw, nCellIndex);
	}
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> ::updateSolutionBlendFactor()
{
	// Field Index
	const int nSigmaIndex = m_pMeshReader->getScalarFieldIndex("sigma0");
	
	// Update Fields
	for (int nCIndex = 0; nCIndex < cells_cfd.size(); nCIndex++)
	{
		// Cell Index
		const int nCellIndex = cells_cfd[nCIndex].m_nCellIndex;
		m_pMeshReader->updateScalarField(nSigmaIndex, cells_cfd[nCIndex].sponge_sigma, nCellIndex);
	}
}


//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Post-Process functions
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: prePostProc( MPI_env &mpi_env,
																bool haveProbes,
																bool haveSampling,
																bool haveAverage,
																bool haveForces,
																vector<probe> &x_probes,
																vector<probe> &x_samples,
																int sm ){
	
	if( haveAverage ){
		m_dpAVG.resize( cells_cfd.size() );
		m_dpRMS.resize( cells_cfd.size() );

		for( size_t icell=0; icell < cells_cfd.size(); icell++ ) {
			m_dpAVG[icell] = ZERO;
			m_dpRMS[icell] = ZERO;
		}
	}

	if( haveForces ){
		m_dFpre.resize( m_nWallBCList.size() );
		m_dFvis.resize( m_nWallBCList.size() );
	}
}

// Average fields
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: postProcAverage( int time_step )
{

	if(time_step > 0){
		for( size_t icell=0; icell < cells_cfd.size(); icell++ ) 
		{
			PRECISION r = cells_cfd[icell].vars.q_old[0];
			PRECISION velMag = 0;
			for (int nD = 1; nD <= DIM_CNT; nD++)
			{
				const PRECISION vel = cells_cfd[icell].vars.q_old[nD] / r;
				velMag += (vel * vel);
			}
			PRECISION E = cells_cfd[icell].vars.q_old[DIM_CNT + 1] / r;
			PRECISION p = ( r * m_dGammaMinusOne ) * (E - 0.5 * velMag);

			m_dpAVG[icell] += p;
			m_dpRMS[icell] += pow(p - m_dpAVG[icell] / time_step, 2);
		}
	}
}

// Average fields OUTPUT
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: postProcAverageOutput( gmsh_mesh *mesh, string filename, PRECISION time, int time_step )
{
	ofstream fout;
	fout.open( filename );

	fout << "TITLE = \"Example: 2D Finite Element Data\"" << endl;
	fout << "VARIABLES = \"X\", \"Y\", \"m_dpAVG\", \"m_dpRMS\"" << endl;

	// Assume triangular elements for now
	if( FACE_CNT == 3 ){
		fout << "ZONE NODES = " << setw(5) << mesh->nodes.size()
		 	 << ", ELEMENTS = " << setw(5) << cells_cfd.size()
		 	 << ", DATAPACKING=BLOCK"
			 << ", SOLUTIONTIME = " << setw(5) << time
		 	 << ", ZONETYPE=FETRIANGLE"
			 << ", VarLocation=([3,4]=CellCentered)"
		 	 << endl;
	}else{
		fout << "ZONE NODES = " << setw(5) << mesh->nodes.size()
		 	 << ", ELEMENTS = " << setw(5) << cells_cfd.size()
		 	 << ", DATAPACKING=BLOCK"
			 << ", SOLUTIONTIME = " << setw(5) << time
		 	 << ", ZONETYPE=FEQUADRILATERAL"
			 << ", VarLocation=([3,4]=CellCentered)"
		 	 << endl;
	}

	// Node coordinates
	for( auto inode=mesh->nodes.begin(); inode < mesh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[0];
	}
	fout << endl;

	for( auto inode=mesh->nodes.begin(); inode < mesh->nodes.end(); inode++ ){
		fout << setw(15) << inode->xn[1];
	}
	fout << endl;

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << m_dpAVG[icell] / time_step;
	}
	fout << endl;

	for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
		fout << setw(15) << std::fixed << std::setprecision(10) << sqrt( m_dpRMS[icell] / time_step );
	}
	fout << endl;

	// Mesh connectivity	
	for( size_t icell = 0; icell < cells_cfd.size(); icell++ )
	{		
		for( unsigned inode = 0; inode < mesh->cells[icell].nodeCount; inode++ )
		{
			fout << setw(15) << mesh->cell2node[icell][inode] + 1;
		}
		fout << endl;
	}
	fout.close();
}

// Calculate forces
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: postProcForces( MPI_env &mpi_env, PRECISION time, MPI_Datatype mpi_precision ){

	m_dtimeVec.push_back( time );

	for ( size_t itype = 0; itype < m_nWallBCList.size(); itype++) {
		array<PRECISION, DIM_CNT> Fpre_local = { 0.0 };
		array<PRECISION, DIM_CNT> Fvis_local = { 0.0 };
		array<PRECISION, DIM_CNT> Fpre_global, Fvis_global;
		for ( size_t icell = 0; icell < boundaries[m_nWallBCList[itype]].size(); icell++ ) {
			int cell = boundaries[m_nWallBCList[itype]][icell][0];
			int face = boundaries[m_nWallBCList[itype]][icell][1];

			auto *param = &cells_cfd[cell];
			auto *vars  = &cells_cfd[cell].vars;
			auto *neigh = cells_cfd[cell].neighs[face];

			PRECISION r = vars->q_old[0];
			PRECISION rU[DIM_CNT] = {0.0};
			for (unsigned idim = 0; idim < DIM_CNT; ++idim){
				rU[idim] = vars->q_old[idim+1];
			}
			PRECISION rE = vars->q_old[DIM_CNT+1];

			PRECISION Uvec[DIM_CNT];
			PRECISION Umag = 0.0;
			PRECISION Smag = 0.0;
			for (unsigned idim = 0; idim < DIM_CNT; ++idim){
				Uvec[idim] = rU[idim] / r;
				Umag += Uvec[idim] * Uvec[idim];
				Smag += cells_cfd[cell].S[face][idim] * cells_cfd[cell].S[face][idim];
			}
			Umag = sqrt( Umag );
			Smag = sqrt( Smag );
			PRECISION E = rE / r;
			PRECISION p = ( r * m_dGammaMinusOne ) * (E - 0.5 * Umag * Umag );

			for (int nD = 0; nD < DIM_CNT; nD++)
				Fpre_local[nD] += cells_cfd[cell].S[face][nD] * (p - m_dpInf);
			
			PRECISION dudx[DIM_CNT][DIM_CNT];

			// explicit correction for physical boundaries (wall/inlet/outlet)
			PRECISION d_mag = vector_mag( param->d[face] );
			PRECISION dmag_inv = ONE / d_mag ;
			PRECISION d_norm[DIM_CNT];
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				d_norm[idim] = param->d[face][idim] * dmag_inv ;
			}
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				PRECISION cell_Uvec = vars ->q_old[idim+1] / vars ->q_old[0];
				PRECISION adjc_Uvec = neigh->q_old[idim+1] / neigh->q_old[0];
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					dudx[idim][jdim] = (adjc_Uvec - cell_Uvec) * d_norm[jdim] * dmag_inv;
				}
			}
			
			PRECISION mu = m_dMu0;
			PRECISION tau[DIM_CNT][DIM_CNT];
			for (int nD = 0; nD < DIM_CNT; nD++){
				tau[nD][nD] = 2.0 * dudx[nD][nD];
				for (int nD1 = nD + 1; nD1 < DIM_CNT + nD; nD1++){
					const int nD2 = nD1 % DIM_CNT;
					tau[nD][nD2] = mu * (dudx[nD][nD2] + dudx[nD2][nD]);
					tau[nD][nD] -= dudx[nD2][nD2];
				}
				tau[nD][nD] *= 2.0 / 3.0 * mu;
			}

			for (int nD1 = 0; nD1 < DIM_CNT; nD1++)
				for (int nD2 = 0; nD2 < DIM_CNT; nD2++)					
					Fvis_local[nD1] -= tau[nD1][nD2] * cells_cfd[cell].S[face][nD2];				
		}

		MPI_Allreduce( &Fpre_local, &Fpre_global, DIM_CNT, mpi_precision, MPI_SUM, MPI_COMM_WORLD );
		MPI_Allreduce( &Fvis_local, &Fvis_global, DIM_CNT, mpi_precision, MPI_SUM, MPI_COMM_WORLD );

		if( mpi_env.is_master() ){
			m_dFpre[itype].push_back( Fpre_global );
			m_dFvis[itype].push_back( Fvis_global );
		}
	}
}

// Forces OUTPUT
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: postProcForcesOutput( MPI_env &mpi_env, ofstream &fout, size_t itype )
{

	for (size_t i = 0; i < m_dtimeVec.size(); i++ ) 
	{
		fout << std::scientific << std::setprecision(16);
		fout << setw(25) << m_dtimeVec[i];
		for (int nD = 0; nD < DIM_CNT; nD++)
			 fout << setw(25) << m_dFpre[itype][i][nD];
		for (int nD = 0; nD < DIM_CNT; nD++)
			 fout << setw(25) << m_dFvis[itype][i][nD];
		for (int nD = 0; nD < DIM_CNT; nD++)
			 fout << setw(25) << m_dFpre[itype][i][nD] + m_dFvis[itype][i][nD];
		fout << endl;
	}

	if ( itype == m_nWallBCList.size() - 1 ) {
		m_dtimeVec.clear();
	}
	m_dFpre[itype].clear();
	m_dFvis[itype].clear();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// New exchange ghost cells:
// 	- pack data into contiguous array for each neighbor
// 	- send/recv contiguous array
// 	- unpack data into ghost cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: new_exchange_ghost( MPI_env &mpi_env ){

	// Serial case
	if( mpi_env.is_serial() ){		
		if ( !ghost_mpi.empty() && !ghost_mpi[0].empty() ) {
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				ghost_mpi[0][icell].vars = cells_cfd[icell].vars;
			}
		}
		return;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// non-MPI packing of variables (as opposed to using MPI_Pack/Unpack functions)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Wait for any open epoch to complete (does nothing is the window was not posted)
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB)
		mpi_env.waitWindow(0);

	int offset = 0;	
	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {
		for (size_t icell=0; icell < local_cells_to_send[ineigh].size(); icell++) {

			int gcell = local_cells_to_send[ineigh][icell];

			// q_old
			for (unsigned idim=0; idim < DIM_CNT+2; idim++)
				m_SendBufList[offset++] = cells_cfd[gcell].vars.q_old[idim];

			// dudx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					m_SendBufList[offset++] = cells_cfd[gcell].vars.dudx[idim][jdim];

			// dTdx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				m_SendBufList[offset++] = cells_cfd[gcell].vars.dTdx[idim];

			// tauMC
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					m_SendBufList[offset++] = cells_cfd[gcell].vars.tauMC[idim][jdim];

			// sigmaU
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				m_SendBufList[offset++] = cells_cfd[gcell].vars.sigmaU[idim];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI communication
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB || mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.postWindow(0);
		mpi_env.startWindow(0);
	}

	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_PERS)
		mpi_env.startAllPersistantSendRecv(0);
	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) 
	{
		switch (mpi_env.get_halo_comm_type()) {
			case MPI_TWOSIDED_NONB:
			case MPI_TWOSIDED_BLCK:
				mpi_env.new_isendrecv<PRECISION>(ineigh, m_nSendBufCountList[ineigh], &m_SendBufList[m_nSendBufOffsetList[ineigh]],
														 m_nRecvBufCountList[ineigh], &m_RecvBufList[m_nRecvBufOffsetList[ineigh]], 0);
				break;				
			case MPI_ONESIDED_NONB:
			case MPI_ONESIDED_BLCK:
				mpi_env.new_irget(ineigh, m_nRecvBufCountList[ineigh], &m_RecvBufList[m_nRecvBufOffsetList[ineigh]], m_nRecvWindowOffsetList[ineigh], 0);
				break;
			case MPI_TWOSIDED_PERS:
				//mpi_env.startPersistantSendRecv(0, ineigh);
				break;
			default:
				cout << "Packed MPI communication currently only supports two-sided and one-sided non-blocking halo comm modes." << endl;
				MPI_Abort( MPI_COMM_WORLD, 909 );
		}
	}

	// Blocking: Wait all
	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_BLCK)
		mpi_env.waitAll(0);		

	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.completeWindow(0);
		mpi_env.waitWindow(0);
	}	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: exchange_ghost_solvars( MPI_env &mpi_env ){

	// Serial case
	if( mpi_env.is_serial() ){		
		if ( !ghost_mpi.empty() && !ghost_mpi[0].empty() ) {
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				ghost_mpi[0][icell].vars = cells_cfd[icell].vars;
			}
		}
		return;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// non-MPI packing of variables (as opposed to using MPI_Pack/Unpack functions)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	// Wait for any open epoch to complete
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB)
		mpi_env.waitWindow(0);

	int offset = 0;
	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {		
		for (size_t icell=0; icell < local_cells_to_send[ineigh].size(); icell++) {

			int gcell = local_cells_to_send[ineigh][icell];

			// q_old
			for (unsigned idim=0; idim < DIM_CNT+2; idim++)
				m_SendBufSolvarsList[offset++] = cells_cfd[gcell].vars.q_old[idim];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI communication
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB || mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.postWindow(0);
		mpi_env.startWindow(0);
	}
	
	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_PERS)
		mpi_env.startAllPersistantSendRecv(0);

	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {
		switch (mpi_env.get_halo_comm_type()) {
			case MPI_TWOSIDED_NONB:
			case MPI_TWOSIDED_BLCK:
				mpi_env.new_isendrecv<PRECISION>(ineigh, m_nSendBufSolvarCountList[ineigh], &m_SendBufSolvarsList[m_nSendBufSolvarOffsetList[ineigh]],
														 m_nRecvBufSolvarCountList[ineigh], &m_RecvBufSolvarsList[m_nRecvBufSolvarOffsetList[ineigh]], 0);
				break;
			case MPI_ONESIDED_NONB:
			case MPI_ONESIDED_BLCK:
				mpi_env.new_irget(ineigh, m_nRecvBufSolvarCountList[ineigh], &m_RecvBufSolvarsList[m_nRecvBufSolvarOffsetList[ineigh]], m_nRecvSolvarsWindowOffsetList[ineigh], 0);
				break;
			case MPI_TWOSIDED_PERS:
				//mpi_env.startPersistantSendRecv(0, ineigh);
				break;
			default:
				cout << "Packed MPI communication currently only supports two-sided and one-sided non-blocking halo comm modes." << endl;
				MPI_Abort( MPI_COMM_WORLD, 909 );
		}
	}

	// Blocking: Wait all
	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_BLCK)
		mpi_env.waitAll(0);

	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.completeWindow(0);
		mpi_env.waitWindow(0);
	}	
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: exchange_ghost_viscous( MPI_env &mpi_env ){

	// Serial case
	if( mpi_env.is_serial() ){		
		if ( !ghost_mpi.empty() && !ghost_mpi[0].empty() ) {
			for( size_t icell=0; icell < cells_cfd.size(); icell++ ){
				ghost_mpi[0][icell].vars = cells_cfd[icell].vars;
			}
		}
		return;
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// non-MPI packing of variables (as opposed to using MPI_Pack/Unpack functions)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Wait window
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB)
		mpi_env.waitWindow(1);

	int offset = 0;
	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {
		for (size_t icell=0; icell < local_cells_to_send[ineigh].size(); icell++) {

			int gcell = local_cells_to_send[ineigh][icell];

			// dudx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					m_SendBufViscousList[offset++] = cells_cfd[gcell].vars.dudx[idim][jdim];

			// dTdx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				m_SendBufViscousList[offset++] = cells_cfd[gcell].vars.dTdx[idim];

			// tauMC
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					m_SendBufViscousList[offset++] = cells_cfd[gcell].vars.tauMC[idim][jdim];

			// sigmaU
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				m_SendBufViscousList[offset++] = cells_cfd[gcell].vars.sigmaU[idim];
		}
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// MPI communication
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_NONB || mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.postWindow(1);
		mpi_env.startWindow(1);
	}
	
	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_PERS)
		mpi_env.startAllPersistantSendRecv(1);

	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {
		switch (mpi_env.get_halo_comm_type()) {
			case MPI_TWOSIDED_NONB:
			case MPI_TWOSIDED_BLCK:
				mpi_env.new_isendrecv<PRECISION>(ineigh, m_nSendBufViscousCountList[ineigh], &m_SendBufViscousList[m_nSendBufViscousOffsetList[ineigh]],
														 m_nRecvBufViscousCountList[ineigh], &m_RecvBufViscousList[m_nRecvBufViscousOffsetList[ineigh]], 1);
				break;
			case MPI_ONESIDED_NONB:
			case MPI_ONESIDED_BLCK:
				mpi_env.new_irget(ineigh, m_nRecvBufViscousCountList[ineigh], &m_RecvBufViscousList[m_nRecvBufViscousOffsetList[ineigh]], m_nRecvViscousWindowOffsetList[ineigh], 1);
				break;
			case MPI_TWOSIDED_PERS:
				//mpi_env.startPersistantSendRecv(1, ineigh);
				break;
			default:
				cout << "Packed MPI communication currently only supports two-sided and one-sided non-blocking halo comm modes." << endl;
				MPI_Abort( MPI_COMM_WORLD, 909 );
		}
	}

	// Blocking: Wait all
	if (mpi_env.get_halo_comm_type() == MPI_TWOSIDED_BLCK)
		mpi_env.waitAll(1);

	if (mpi_env.get_halo_comm_type() == MPI_ONESIDED_BLCK)
	{
		mpi_env.completeWindow(1);
		mpi_env.waitWindow(1);
	}	
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: mpi_communication(MPI_env &mpi_env, int comm_step){

	// Original MPI comm: send full bnd Hpath to all neighbors
	if (mpi_env.m_nCommType == MPI_EXCHANGE_FULL_BND) {
		exchange_ghost_cells(mpi_env);
		return;
	}

	// Packed MPI comm: send only necessary cells to each neighbor
	if (mpi_env.m_nCommType == MPI_EXCHANGE_PACKED) {
		new_exchange_ghost(mpi_env);
		return;
	}

	// Split MPI comm into two steps:
	// 	- exchange only solvars at end of time-step
	// 	- exchange only viscous vars after viscous computation
	if (mpi_env.m_nCommType == MPI_EXCHANGE_SPLIT) {
		if (comm_step == 0 )
			exchange_ghost_solvars(mpi_env);
		else
			exchange_ghost_viscous(mpi_env);
		return;
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: mpi_wait(MPI_env &mpi_env, int comm_step)
{
	// WinIndex
	const int nWinIndex = (mpi_env.m_nCommType == MPI_EXCHANGE_SPLIT) ? comm_step : 0;

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Wait for ghost cells
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	mpi_env.halo_comm_wait(nWinIndex);

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Unpack Data
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Packed MPI exchange
	if (mpi_env.m_nCommType == MPI_EXCHANGE_PACKED)
	{		
		unpack_mpi_data(mpi_env);
	}

	// Split MPI exchange
	if (mpi_env.m_nCommType == MPI_EXCHANGE_SPLIT) {
		if (comm_step == 0)
			unpack_mpi_data_solvars(mpi_env);
		else
			unpack_mpi_data_viscous(mpi_env);
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: unpack_mpi_data(MPI_env &mpi_env){

	int offset = 0;

	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {		
		for (size_t icell=0; icell < neigh_cells_to_recv[ineigh].size(); icell++) {

			int gcell = neigh_cells_to_recv[ineigh][icell];


			// q_old
			for (unsigned idim=0; idim < DIM_CNT+2; idim++)
				ghost_mpi[ineigh][gcell].vars.q_old[idim] = m_RecvBufList[offset++];

			// dudx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
				 	ghost_mpi[ineigh][gcell].vars.dudx[idim][jdim] = m_RecvBufList[offset++];

			// dTdx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				ghost_mpi[ineigh][gcell].vars.dTdx[idim] = m_RecvBufList[offset++];

			// tauMC
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					ghost_mpi[ineigh][gcell].vars.tauMC[idim][jdim] = m_RecvBufList[offset++];

			// sigmaU
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				ghost_mpi[ineigh][gcell].vars.sigmaU[idim] = m_RecvBufList[offset++];
		}
	}
}

//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: unpack_mpi_data_solvars(MPI_env &mpi_env){

	int offset = 0;

	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {
		for (size_t icell=0; icell < neigh_cells_to_recv[ineigh].size(); icell++) 
		{
			int gcell = neigh_cells_to_recv[ineigh][icell];

			// q_old
			for (unsigned idim=0; idim < DIM_CNT+2; idim++)
				ghost_mpi[ineigh][gcell].vars.q_old[idim] = m_RecvBufSolvarsList[offset++];
		}
	}
}
//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
void CFDv0_solver<PRECISION, DIM_CNT, FACE_CNT> :: unpack_mpi_data_viscous(MPI_env &mpi_env){

	int offset = 0;

	for (size_t ineigh=0; ineigh < mpi_env.mpi_neighbors.size(); ineigh++) {		
		for (size_t icell=0; icell < neigh_cells_to_recv[ineigh].size(); icell++) {

			int gcell = neigh_cells_to_recv[ineigh][icell];

			// dudx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
				 	ghost_mpi[ineigh][gcell].vars.dudx[idim][jdim] = m_RecvBufViscousList[offset++];

			// dTdx
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				ghost_mpi[ineigh][gcell].vars.dTdx[idim] = m_RecvBufViscousList[offset++];

			// tauMC
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				for (unsigned jdim=0; jdim < DIM_CNT; jdim++)
					ghost_mpi[ineigh][gcell].vars.tauMC[idim][jdim] = m_RecvBufViscousList[offset++];

			// sigmaU
			for (unsigned idim=0; idim < DIM_CNT; idim++)
				ghost_mpi[ineigh][gcell].vars.sigmaU[idim] = m_RecvBufViscousList[offset++];
		}
	}
}
