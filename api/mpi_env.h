#ifndef MPI_ENV_H
#define MPI_ENV_H

#define MAX_MPI_RANKS 		10000
#define MAX_MPI_NEIGHBORS 	1000
#define MAX_CELL_FACE_COUNT 10

#define GETTIME(T) clock_gettime(CLOCK_MONOTONIC, &T)
#define DIFFTIME(Te,Ts) ((Te.tv_sec)-(Ts.tv_sec)+1e-9*((Te.tv_nsec)-(Ts.tv_nsec)))

#ifdef EXTRAE_PROFILING
	#include "extrae.h"
	#define EXTRAE_RESTART Extrae_restart();
	#define EXTRAE_EVENT(eventId,eventValue) Extrae_event((eventId),(eventValue));
	#define EXTRAE_FINI Extrae_fini();
	#define EXTRAE_SHUTDOWN Extrae_shutdown();
#else
	#define EXTRAE_RESTART
	#define EXTRAE_EVENT(eventId,eventValue)
	#define EXTRAE_FINI
	#define EXTRAE_SHUTDOWN
#endif

#include "mpi.h"
#include "elements.h"
#include "mpi_env.h"
#include "cfdv0_elements.h"

struct t_mpi_neighbors{
	int id_local;
	int id_global;
};

enum t_mpi_comm_type {
	MPI_EXCHANGE_FULL_BND = 0,	// Exchange full bnd Hpath will all neighbors
	MPI_EXCHANGE_PACKED   = 1,	// Exchange only necessary cells with each neighbor
	MPI_EXCHANGE_SPLIT    = 2,	// Exchange only necessary vars with each neighbor: two sorts of MPI comm type
	MPI_EXCHANGE_COUNT 	  = 3,
};

enum t_mpi_halo_comm {
    MPI_TWOSIDED_BLCK  = 0, // Plain cell
    MPI_TWOSIDED_NONB  = 1, // Neighbor is in a neighboring submesh
    MPI_TWOSIDED_PERS  = 2, // Persistante two-sided 
    MPI_ONESIDED_BLCK  = 3, // Plain cell
    MPI_ONESIDED_NONB  = 4, // Neighbor is in a neighboring submesh
    MPI_ONESIDED_PERS  = 5, // Neighbor is in a neighboring submesh
    MPI_NEIGHCOLL_BLCK = 6,
    MPI_NEIGHCOLL_NONB = 7,
    MPI_NEIGHCOLL_PERS = 8,
    MPI_NONE           = 9, // Error
	MPI_TYPE_COUNT     = 10,
};


class t_mpi_comm{

	public:
		// Constructors
		t_mpi_comm(){}
		t_mpi_comm( MPI_Comm tmp_comm, int tmp_size, int tmp_rank ){
			m_comm = tmp_comm;
			m_size = tmp_size;
			m_rank = tmp_rank;
		}

		// Communicator
		void set_comm( MPI_Comm tmp_comm ){
			m_comm = tmp_comm;
		}
		MPI_Comm self(){
			return m_comm;
		}

		// Size
		void set_size( int tmp_size ){
			m_size = tmp_size;
		}
		int size(){
			return m_size;
		}

		// Rank
		void set_rank( int tmp_rank ){
			m_rank = tmp_rank;
		}
		int rank(){
			return m_rank;
		}

		// Master
		int master(){
			return 0;
		}

	private:
		MPI_Comm m_comm;
		int m_size;
		int m_rank;
};

//***************************************************************************************************
// MPI environment class
// 		- handles all calls to MPI functions/routines
// 		- setup and wrapping of MPI environment
// 		- MPI datatypes, buffer communication, collectives...
//***************************************************************************************************
class MPI_env {

	public:

		// Initialize MPI environment
		void initialize( bool is_solver_sim );
		void finalize();

		//MPI_Datatype get_regular_mesh_MPI_datatype( t_mesh *mesh );
		//void test_mpi_exchange( t_mesh *mesh, MPI_Datatype type_node );

		int size(){
			return nProcs;
		}

		int rank(){
			return myrank;
		}

		int get_master(){
			return 0;
		}

		bool is_master(){
			return master;
		}

		int get_node_master(){
			return node_master_rank;
		}

		bool is_node_master(){
			return node_master;
		}

		bool is_serial(){
			return (nProcs == 1);
		}

		void set_cell_type( MPI_Datatype tmp_cell_type ){ type_cell = tmp_cell_type; };
		MPI_Datatype cell_type()      { return type_cell;       }
		MPI_Comm     graph_comm()     { return comm_graph;      }

		void barrier(void)          { MPI_Barrier( MPI_COMM_WORLD ); }
		void abort( int error_code ){ MPI_Abort  ( MPI_COMM_WORLD, error_code ); }
		void broadcast_partitions( std::vector<int> partitions );

		// Initialize MPI graph topology
		void init_graph( std::vector<int> mpi_neighbors );

		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void set_halo_comm_type( enum t_mpi_halo_comm type_mpi_comm,
									std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>& cells_cfd,
									std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>>& ghost_cells );

		void set_halo_comm_type_onesided(enum t_mpi_halo_comm type_mpi_comm,
										std::vector<void*>& pBufferList, std::vector<size_t>& nBufferElementCountList, 
										std::vector<size_t>& nBufferElementSizeList);

		enum t_mpi_halo_comm get_halo_comm_type(){ return halo_comm_type; }

		void set_mpi_neighbors( std::vector<int> m_mpi_neighbors ){
			mpi_neighbors = m_mpi_neighbors;
		}
		void set_neigh_cells(std::vector<std::vector<int>> m_mpi_neigh_cells_send, std::vector<std::vector<int>> m_mpi_neigh_cells_recv) {
			neigh_cells_send = m_mpi_neigh_cells_send;
			neigh_cells_recv = m_mpi_neigh_cells_recv;
		}

		void print_from_master( std::string print_statement ){
			barrier();
			if( is_master() ) cout << print_statement << endl;
		}

		void set_mesh_size_type( MPI_Datatype type_mesh_size ){ type_meshsize = type_mesh_size; }
		MPI_Datatype get_type_meshsize(){ return type_meshsize; }

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Wait functions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		void halo_comm_wait(const int nWinIndex);
		void halo_comm_test( int *flag, const int nWinIndex);

		// Neighborhood collective
		void save_neighborhood_request( MPI_Request request_hpath, MPI_Request request_deadend ){
			neighcoll_req_hpath   = request_hpath   ;
			neighcoll_req_deadend = request_deadend ;
		}


		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// CFD Solver
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// One-sided, blocking
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void rget(std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath );

		// One-sided, non-blocking
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void irget(std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath );

		// Two-sided, blocking
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void sendrecv(          std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_hpath,
                       std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath );

		// Two-sided, non-blocking
		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void isendrecv(          std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_hpath,
                       std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath );

		template <typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT>
		void initializeCompletePersistantSendRecv(const int ineigh, std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>  &cells_hpath,
                       											    std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> &ghost_hpath);
		
		// New functions
		inline void set_type_cfdcell( MPI_Datatype t_cfdcell ){ type_cfdcell = t_cfdcell; }
		inline void set_halo_comm_reqs( int n )
		{
			m_CommRequestList.resize(2);
			for (int nIndex = 0; nIndex < 2; nIndex++)
			{
				m_CommRequestList[nIndex].resize(n, MPI_REQUEST_NULL);
			}
		};

		std::vector<int> mpi_neighbors;
		MPI_Group m_NeighborGroup;
		std::vector<std::vector<int>> neigh_cells_send;
		std::vector<std::vector<int>> neigh_cells_recv;
		// std::vector<MPI_Request> comm_reqs;		
		// std::vector<MPI_Status> comm_status;
		std::vector<std::vector<MPI_Request> > m_CommRequestList;
		std::vector<std::vector<MPI_Status> > m_CommStatusList;
		std::vector<int> node_local_ranks;

		// Node communicator
		class t_mpi_comm node_comm;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// New MPI Comm type:
		// 	- true: exchange packed data + only required cells to reach neighbor
		// 	- false: exchange whole cell + all boundary cells to all neighbors
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		enum t_mpi_comm_type m_nCommType = MPI_EXCHANGE_SPLIT;

		// Two-sided, non-blocking
		template <typename PRECISION>
		void new_isendrecv(int ineigh, int send_size, PRECISION *send_buf, int recv_size, PRECISION *recv_buf, const int nWinIndex);

		// Initialize persistante two-sided communication
		template <typename PRECISION>
		void initializePersistantSendRecv(int ineigh, int send_size, PRECISION *send_buf, int recv_size, PRECISION *recv_buf, const int nWinIndex);
		
		void startAllPersistantSendRecv(const int nWinIndex);
		void startPersistantSendRecv(const int nWinIndex, const int ineigh);

		// One-sided
		void postWindow(const int nWinIndex);
		void startWindow(const int nWinIndex);
		void waitWindow(const int nWinIndex);
		void completeWindow(const int nWinIndex);

		void lock_windows(const int nWinIndex);

		template <typename PRECISION>
		void new_irget(int ineigh, int recv_size, PRECISION *recv_buf, int nWindowOffset, int nWinIndex);

		void waitAll(const int nWinIndex);
	private:
		int  nProcs;
		int  myrank;
		bool master;
		int  node_master_rank;
		bool node_master;

		MPI_Datatype type_cell;
		MPI_Comm comm_graph;
		MPI_Request neighcoll_req_hpath;
		MPI_Request neighcoll_req_deadend;

		MPI_Datatype type_cfdcell;

		std::vector<MPI_Request> twosided_send_req_hpath, twosided_send_req_reg;
		std::vector<MPI_Request> twosided_recv_req_hpath, twosided_recv_req_reg;

		int          *sendcounts, *recvcounts;
		MPI_Aint     *senddisps , *recvdisps ;
		MPI_Datatype *sendtypes , *recvtypes ;

		int          *sendcounts_reg, *recvcounts_reg;
		MPI_Aint     *senddisps_reg , *recvdisps_reg ;
		MPI_Datatype *sendtypes_reg , *recvtypes_reg ;

		MPI_Request neighcoll_req;

		std::vector<MPI_Win> m_WindowList;
		std::vector<bool> m_bWindowLockedList;
		std::vector<bool> m_bWindowIsPostList;
		std::vector<bool> m_bWindowIsStartList;

		enum t_mpi_halo_comm halo_comm_type = MPI_NONE;

		MPI_Aint m_nCellSize;
		MPI_Datatype type_meshsize;
};

#endif







































