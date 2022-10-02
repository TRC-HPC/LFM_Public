#include <iomanip>

#include "api/mpi_env.h"
#include <cstddef>
#include <cmath>

#include "api/templates_mpienv.h"

#ifdef OPEN_MPI			// check for OMPI
#include "mpi-ext.h"
#if OMPI_MAJOR_VERSION < 5	// redirect to experimental neighborhood collectives on OMPI < 5
#endif
#else				// neigborhood collectives only in MPI
	#pragma message("Warning: Neighborhood collectives are only available in OpenMPI. Compiling without.")
#endif

//***************************************************************************************************
void MPI_env::initialize( bool is_solver_sim )
{
	int init;
	MPI_Initialized(&init);
	if(!init)
		MPI_Init(NULL, NULL);

	// Initializ global rank and number of processors
	MPI_Comm_rank( MPI_COMM_WORLD, &myrank );
	MPI_Comm_size( MPI_COMM_WORLD, &nProcs );

	// Check if we are the master (in global)
	master = (myrank == 0);

	if( is_solver_sim ){
		return;
	}

	// Get node master
	char local_proc_name[MPI_MAX_PROCESSOR_NAME];
	int   name_length;
	MPI_Get_processor_name( local_proc_name, &name_length );

	// Get processor name of each node
	std::string name( local_proc_name );
	char *global_proc_name = new char[ this->size() * MPI_MAX_PROCESSOR_NAME ];
	MPI_Allgather( &local_proc_name, name_length, MPI_CHAR, &global_proc_name[0], name_length, MPI_CHAR, MPI_COMM_WORLD );

	// Convert processor name to std::string
	std::vector<std::string> all_proc_names( this->size(), "" );
	for( int irank=0; irank < this->size(); irank++ ){
		for( int ichar=0; ichar < name_length; ichar++ ){
			all_proc_names[irank] += global_proc_name[ichar + name_length*irank];
		}
	}

	// Initialize invalid master node (among nodes with the same processor)
	int node_master_rank = this->size();
	std::vector<bool> is_node_local( this->size(), false );

	for( int irank=0; irank < this->size(); irank++ ){
		if( all_proc_names[this->rank()].compare( all_proc_names[irank] ) == 0 ){

			// Add node-local MPI ranks
			node_local_ranks.push_back( irank );
			is_node_local[irank] = true;
			
			// Identify node master
			if( node_master_rank > irank ){
				node_master_rank = irank;
			}
			node_master_rank = node_master_rank > irank ? irank : node_master_rank;
		}
	}

	// Check if we are the master node within the same-processor nodes
	if( myrank == node_master_rank ){
		node_master = true;
	}else{
		node_master = false;
	}

	// Prepare node-local communicator
	// Split according to the node master rank (i.e. all nodes with same local master will share this communicator)
	MPI_Comm node_comm;
	MPI_Comm_split( MPI_COMM_WORLD, node_master_rank, this->rank(), &node_comm );

	// Get number of nodes and local rank within the local communicator
	int node_size, node_rank;
	MPI_Comm_size( node_comm, &node_size );
	MPI_Comm_rank( node_comm, &node_rank );

	this->node_comm.set_comm( node_comm );
	this->node_comm.set_size( node_size );
	this->node_comm.set_rank( node_rank );
}

//***************************************************************************************************
void MPI_env::finalize(){
	
	switch( halo_comm_type ){

		case MPI_ONESIDED_NONB:
		case MPI_ONESIDED_BLCK:
			for (int nI = 0; nI < m_WindowList.size(); nI++)
		    	MPI_Win_free(&m_WindowList[nI]);
			break;

		default:
			break;
	}

	MPI_Finalize();
}

//***************************************************************************************************
void MPI_env :: set_halo_comm_type_onesided(enum t_mpi_halo_comm type_mpi_comm,
										std::vector<void*>& pBufferList, std::vector<size_t>& nBufferElementCountList, 
										std::vector<size_t>& nBufferElementSizeList)

{
	halo_comm_type = type_mpi_comm;	
	m_WindowList.resize(pBufferList.size());
	for (int nIndex = 0; nIndex < pBufferList.size(); nIndex++)
	{
		// Create Window
		m_nCellSize = nBufferElementSizeList[nIndex];
		MPI_Aint dataSize = nBufferElementCountList[nIndex] * nBufferElementSizeList[nIndex];
		MPI_Win_create(pBufferList[nIndex], dataSize, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &m_WindowList[nIndex]);
		m_bWindowLockedList.push_back(false);
		m_bWindowIsPostList.push_back(false);
		m_bWindowIsStartList.push_back(false);
	}

	// Create group with all neighbours
    MPI_Group global_group;
	MPI_Comm_group(MPI_COMM_WORLD,&global_group);	
    MPI_Group_incl(global_group, mpi_neighbors.size(), mpi_neighbors.data(), &m_NeighborGroup);	
	MPI_Group_free(&global_group);
}

//***************************************************************************************************
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env :: set_halo_comm_type( enum t_mpi_halo_comm type_mpi_comm,
									fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>& cells_cfd,
									fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>>& ghost_cells ){

	halo_comm_type = type_mpi_comm;

	m_nCellSize = sizeof(cells_cfd[0]);

	// Neighborhood collectives
	switch( halo_comm_type )
	{
		case MPI_ONESIDED_NONB:
		case MPI_ONESIDED_BLCK:
		{			
			m_WindowList.resize(1);

			MPI_Aint dataSize = cells_cfd.size() * m_nCellSize;
		    MPI_Win_create(&cells_cfd[0], dataSize, 1, MPI_INFO_NULL, MPI_COMM_WORLD, &m_WindowList[0]);
			m_bWindowLockedList.push_back(false);
		}
			break;

		default:
			break;
	}

}

//***************************************************************************************************
// MPI_Wait/Waitall
// 		- can be used with any non-blocking and persistent calls
// 		- when blocking calls are used, simply returns without doing anything
//***************************************************************************************************
void MPI_env:: halo_comm_wait(const int nWinIndex)
{
	// Serial simulations
	if( this->is_serial() ) return;

	switch( halo_comm_type ){

		case MPI_TWOSIDED_BLCK:
		case MPI_ONESIDED_BLCK:
		case MPI_NEIGHCOLL_BLCK:
			return;

		case MPI_TWOSIDED_NONB:
		case MPI_TWOSIDED_PERS:
		case MPI_ONESIDED_PERS:
			waitAll(nWinIndex);
			break;
		case MPI_ONESIDED_NONB:
			completeWindow(nWinIndex);			
			break;

		case MPI_NEIGHCOLL_NONB:
			MPI_Waitall( 1, &m_CommRequestList[nWinIndex][0], MPI_STATUS_IGNORE );
			break;
		case MPI_NEIGHCOLL_PERS:
			MPI_Waitall( 1, &neighcoll_req, MPI_STATUS_IGNORE );
			break;
		case MPI_NONE:
			cout << "Error with halo communication type" << endl;
			MPI_Abort( MPI_COMM_WORLD, 928 );
	}
}

//***************************************************************************************************
// MPI_Wait/Waitall
// 		- can be used with any non-blocking and persistent calls
// 		- when blocking calls are used, simply returns without doing anything
//***************************************************************************************************
void MPI_env:: halo_comm_test( int *flag, const int nWinIndex){

	// Serial simulations
	if( this->is_serial() ) return;

	// Depending on communication type...
	switch( halo_comm_type ){

		case MPI_TWOSIDED_BLCK:
		case MPI_ONESIDED_BLCK:
		case MPI_NEIGHCOLL_BLCK:
			return;

		case MPI_ONESIDED_PERS:
		case MPI_ONESIDED_NONB:
		case MPI_TWOSIDED_PERS:
		case MPI_TWOSIDED_NONB:

			MPI_Testall( m_CommRequestList[nWinIndex].size(), &m_CommRequestList[nWinIndex][0], flag, MPI_STATUS_IGNORE );
			break;

		case MPI_NEIGHCOLL_NONB:
			MPI_Waitall( 1, &m_CommRequestList[nWinIndex][0], MPI_STATUS_IGNORE );
			break;
		case MPI_NEIGHCOLL_PERS:
			MPI_Waitall( 1, &neighcoll_req, MPI_STATUS_IGNORE );
			break;
		case MPI_NONE:
			cout << "Error with halo communication type" << endl;
			MPI_Abort( MPI_COMM_WORLD, 928 );
	}
}

//***************************************************************************************************
void MPI_env:: init_graph( const std::vector<int> mpi_neighbors ){

	// Initialize graph topology
	int  src_size  = (int) mpi_neighbors.size();			// Number of neighbors
	int* src_ids   = new int[ mpi_neighbors.size() ];
	int* src_wghts = new int[ mpi_neighbors.size() ];

	int  des_size  = (int) mpi_neighbors.size();			// Number of neighbors
	int* des_ids   = new int[ mpi_neighbors.size() ];
	int* des_wghts = new int[ mpi_neighbors.size() ];

	bool allow_reorder = true;

	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){
		src_ids[ineigh] = mpi_neighbors[ineigh];
		des_ids[ineigh] = mpi_neighbors[ineigh];

		src_wghts[ineigh] = 1;
		des_wghts[ineigh] = 1;
	}

	MPI_Dist_graph_create_adjacent( MPI_COMM_WORLD,
									src_size, src_ids, src_wghts,
									des_size, des_ids, des_wghts,
									MPI_INFO_NULL, allow_reorder, &comm_graph );

	int* test_src_ids   = new int[ mpi_neighbors.size() ];
	int* test_src_wghts = new int[ mpi_neighbors.size() ];
	int* test_des_ids   = new int[ mpi_neighbors.size() ];
	int* test_des_wghts = new int[ mpi_neighbors.size() ];

	MPI_Dist_graph_neighbors( comm_graph, src_size, test_src_ids, test_src_wghts,
										  des_size, test_des_ids, test_des_wghts );

	delete[] src_ids;
	delete[] src_wghts;
	delete[] des_ids;
	delete[] des_wghts;

	delete[] test_src_ids;
	delete[] test_src_wghts;
	delete[] test_des_ids;
	delete[] test_des_wghts;
}

//***************************************************************************************************
//***************************************************************************************************
// 									CFD cells
//***************************************************************************************************
//***************************************************************************************************

// irget
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: irget(fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells )
{
	// Start Window
	startWindow(0);
	
	// MPI_Win_lock_all(MPI_MODE_NOCHECK, window);
	m_bWindowLockedList[0] = true;
	MPI_Request* request = &m_CommRequestList[0][0];
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ )
	{
		int src  = mpi_neighbors[ineigh];
		int neigh_size = ghost_cells[ineigh].size() * m_nCellSize;
		MPI_Rget( &ghost_cells[ineigh][0], neigh_size, MPI_CHAR, src, 0, neigh_size, MPI_CHAR, m_WindowList[0], request++);
	}
}

// rget
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: rget(fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells )
{
	// Start Windows
	startWindow(0);
	
	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int src  = mpi_neighbors[ineigh];
		int neigh_size = ghost_cells[ineigh].size() * m_nCellSize;
		MPI_Get( &ghost_cells[ineigh][0], neigh_size, MPI_CHAR, src, 0, neigh_size, MPI_CHAR, m_WindowList[0]);
	}

	// Complete Window
	completeWindow(0);
}

// sendrecv
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: sendrecv(
						  fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_cfd,
				fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells ){

	const int million = 1000;

	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int src  = mpi_neighbors[ineigh];
		int dest = mpi_neighbors[ineigh];

		int send_tag = 0;
		int recv_tag = 0;

		int local_size = (int) cells_cfd.size();
		int neigh_size = (int) ghost_cells[ineigh].size();

		MPI_Isend( &cells_cfd[0]          , local_size, type_cfdcell, dest, send_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh] );
		MPI_Irecv( &ghost_cells[ineigh][0], neigh_size, type_cfdcell, src , recv_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh+mpi_neighbors.size()] );
	}
	
	MPI_Waitall( m_CommRequestList[0].size(), &m_CommRequestList[0][0], MPI_STATUS_IGNORE );
	for( size_t ineigh=0; ineigh < m_CommRequestList[0].size(); ineigh++ )
		m_CommRequestList[0][ineigh] = MPI_REQUEST_NULL; 

}

// Isendrecv
template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: isendrecv(
						  fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_cfd,
					      fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells ){



	for( size_t ineigh=0; ineigh < mpi_neighbors.size(); ineigh++ ){

		int src  = mpi_neighbors[ineigh];
		int dest = mpi_neighbors[ineigh];

		int send_tag = 0;
		int recv_tag = 0;

		int local_size = (int) cells_cfd.size() * m_nCellSize;
		int neigh_size = (int) ghost_cells[ineigh].size() * m_nCellSize;
		
		MPI_Isend( &cells_cfd[0]          , local_size, MPI_CHAR, dest, send_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh] );
		MPI_Irecv( &ghost_cells[ineigh][0], neigh_size, MPI_CHAR, src , recv_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh+mpi_neighbors.size()] );
	}
}

template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
void MPI_env:: initializeCompletePersistantSendRecv(int ineigh,
						  fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_cfd,
					      fm_vector<fm_vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_cells ){

	int src = mpi_neighbors[ineigh];
	int dest = mpi_neighbors[ineigh];

	int send_tag = 0;
	int recv_tag = 0;

	int local_size = (int)cells_cfd.size() * m_nCellSize;
	int neigh_size = (int)ghost_cells[ineigh].size() * m_nCellSize;

	MPI_Send_init(&cells_cfd[0], local_size, MPI_CHAR, dest, send_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh]);
	MPI_Recv_init(&ghost_cells[ineigh][0], neigh_size, MPI_CHAR, src, recv_tag, MPI_COMM_WORLD, &m_CommRequestList[0][ineigh + mpi_neighbors.size()]);
}


// packed Isendrecv
template <typename PRECISION>
void MPI_env:: new_isendrecv(int ineigh, int send_size, PRECISION *send_buf, int recv_size, PRECISION *recv_buf, const int nWinIndex){

	MPI_Datatype mpi_precision;
	if (sizeof(PRECISION) == sizeof(float))
		mpi_precision = MPI_FLOAT;
	else if (sizeof(PRECISION) == sizeof(double))
		mpi_precision = MPI_DOUBLE;

	int src  = mpi_neighbors[ineigh];
	int dest = mpi_neighbors[ineigh];

	int send_tag = 0;
	int recv_tag = 0;

	MPI_Isend(send_buf, send_size, mpi_precision, dest, send_tag, MPI_COMM_WORLD, &m_CommRequestList[nWinIndex][ineigh]);
	MPI_Irecv(recv_buf, recv_size, mpi_precision, src , recv_tag, MPI_COMM_WORLD, &m_CommRequestList[nWinIndex][ineigh+mpi_neighbors.size()]);
}

template <typename PRECISION>
void MPI_env:: initializePersistantSendRecv(int ineigh, int send_size, PRECISION *send_buf, int recv_size, PRECISION *recv_buf, const int nWinIndex){

	MPI_Datatype mpi_precision;
	if (sizeof(PRECISION) == sizeof(float))
		mpi_precision = MPI_FLOAT;
	else if (sizeof(PRECISION) == sizeof(double))
		mpi_precision = MPI_DOUBLE;

	int src  = mpi_neighbors[ineigh];
	int dest = mpi_neighbors[ineigh];

	int send_tag = nWinIndex;
	int recv_tag = nWinIndex;

	MPI_Send_init(send_buf, send_size, mpi_precision, dest, send_tag, MPI_COMM_WORLD, &m_CommRequestList[nWinIndex][ineigh]);
	MPI_Recv_init(recv_buf, recv_size, mpi_precision, src , recv_tag, MPI_COMM_WORLD, &m_CommRequestList[nWinIndex][ineigh + mpi_neighbors.size()]);
}

void MPI_env::startAllPersistantSendRecv(const int nWinIndex)
{
	MPI_Startall(m_CommRequestList[nWinIndex].size(), &m_CommRequestList[nWinIndex][0]);
}

void MPI_env::startPersistantSendRecv(const int nWinIndex, const int ineigh)
{
	MPI_Start(&m_CommRequestList[nWinIndex][ineigh]);
	MPI_Start(&m_CommRequestList[nWinIndex][ineigh + mpi_neighbors.size()]);
}

void MPI_env::waitAll(const int nWinIndex)
{
	const int nCount = m_CommRequestList[nWinIndex].size();
	MPI_Waitall(nCount, &m_CommRequestList[nWinIndex][0], MPI_STATUS_IGNORE );

	if (halo_comm_type != MPI_TWOSIDED_PERS)
	for( size_t ineigh = 0; ineigh < nCount; ineigh++ )
		m_CommRequestList[nWinIndex][ineigh] = MPI_REQUEST_NULL; 
}

// packed Rget
template <typename PRECISION>
void MPI_env:: new_irget(int ineigh, int recv_size, PRECISION *recv_buf, int nWindowOffset, int nWinIndex)
{	
	int src = mpi_neighbors[ineigh];
	const int nRecvByte = recv_size * sizeof(PRECISION);
	MPI_Get(recv_buf, nRecvByte, MPI_CHAR, src, nWindowOffset * sizeof(PRECISION), nRecvByte, MPI_CHAR, m_WindowList[nWinIndex]);
}


void MPI_env::postWindow(const int nWinIndex)
{
	if (!m_bWindowIsPostList[nWinIndex])
	{
		MPI_Win_post(m_NeighborGroup, MPI_MODE_NOPUT, m_WindowList[nWinIndex]);
		m_bWindowIsPostList[nWinIndex] = true;
	}	
}

void MPI_env::startWindow(const int nWinIndex)
{
	if (!m_bWindowIsStartList[nWinIndex])
	{
		MPI_Win_start(m_NeighborGroup, 0, m_WindowList[nWinIndex]);
		m_bWindowIsStartList[nWinIndex] = true;
	}
	
}

void MPI_env::waitWindow(const int nWinIndex)
{	
	if (m_bWindowIsPostList[nWinIndex])
	{
		MPI_Win_wait(m_WindowList[nWinIndex]);
		m_bWindowIsPostList[nWinIndex] = false;
	}    
}

void MPI_env::completeWindow(const int nWinIndex)
{
	if (m_bWindowIsStartList[nWinIndex])
	{
		MPI_Win_complete(m_WindowList[nWinIndex]);
		m_bWindowIsStartList[nWinIndex] = false;
	}
	
}

// One-sided communication: lock windows
void MPI_env:: lock_windows(const int nWinIndex)
{
	for (size_t ineigh = 0; ineigh < mpi_neighbors.size(); ineigh++)
	{
		int src = mpi_neighbors[ineigh];
		MPI_Win_lock(MPI_LOCK_SHARED, src, MPI_MODE_NOCHECK, m_WindowList[nWinIndex]);
	}

	m_bWindowLockedList[nWinIndex] = true;
}
