#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <vector>
#include <iostream>
#include "mpi.h"
#define FASHMESH_MAX_REQUESTS (6)
#define PACKED        __attribute__((__packed__))
#define SIMD_ALIGNED  __attribute__((__aligned__(8)))
#define CACHE_ALIGNED __attribute__((__aligned__(64)))
#define CACHE_ALIGNED_LONG __attribute__((__aligned__(128)))
#define pow2(x)       ((x) * (x))

using namespace std;

const int INDEX_BND_SUBMESH = 0;

// Preprocess type (cells renumbering function)
enum reorder_cells_type{
	REORDER_NONE	= 0,
	REORDER_HPATH	= 1,
};

// Periodic faces
enum periodic_neigh_type{
	PERIODIC_NEIGH_NONE    = 0,
	PERIODIC_NEIGH_REGULAR = 1,
	PERIODIC_NEIGH_MPI     = 2,
};

struct fastmesh_periodic_face{
	enum periodic_neigh_type type;			// Possible options:
	int orig_id;
	int index;
	int rank ;
	int cell_id ;
};

enum neigh_type {
	CELL_REGULAR  		= 0,			// Neighbor is a regular cell begonging to the same submesh
	CELL_SUBMESH  		= 1,			// Neighbor is a cell belonging to a neighboring interior submesh
	CELL_MPI     		= 2,			// Neighbor cell belongs to another MPI rank
	CELL_PERIODIC		= 3,			// Periodic local neighbor
	CELL_PERIODIC_MPI	= 4,			// Periodic MPI   neighbor
	CELL_NONE     		= 5,			// Neighbor cell does not exist (out of domain)
};

enum bound_type {
	BOUND_NONE		= 0,
	BOUND_WALL		= 1,			// Wall boundary conditions
	BOUND_INLET		= 2,			// Inlet boundary conditions
	BOUND_OUTLET	= 3,			// Outlet boundary conditions
	BOUND_FARFIELD	= 4,			// Far-field boundary conditions
};

enum interpolation_scheme {
	LINEAR_INTERPOLATION		= 0,			// Linear cell to face interpolation
	UPWIND_INTERPOLATION		= 1,			// Upwind cell to face interpolation
	MINMOD_INTERPOLATION		= 2,			// Minmod cell to face interpolation
};

enum case_type {
	TUNNEL						= 0,			// Invented case to do tests for NSE solver
	ISENTROPIC_VORTEX			= 1,			// Isentropic Vortex (Euler)
	DOUBLE_SHEAR_LAYER_1		= 2,			// Double Shear Layer (Euler and NS)
	DOUBLE_SHEAR_LAYER_2		= 3,			// Double Shear Layer (Euler and NS)
	TAYLOR_GREEN_VORTEX			= 4,			// Taylor-Green Vortex (NS)
	CYLINDER_TUNNEL				= 5,			// Cylinder in tunnel (NS)
	CYLINDER					= 6,			// Cylinder case (NS)
	CYLINDER_with_VORTEX		= 7,			// Cylinder with vortex to invoke the shading vortex (NS)
	TANDEM						= 8,			// Two Squares in tandem configuration
	TANDEM_with_VORTEX			= 9,			// Two Squares in tandem configuration with vortex (NS)
	SBS_with_VORTEX_INPHASE		= 10,			// Two Squares in side-by-side configuration with vortex and in-phase results (NS)
	SBS_with_VORTEX_ANTIPHASE	= 11,			// Two Squares in side-by-side configuration with vortex and anti-phase results (NS)
	UNKNOWN_CASE				= 100,			// Unknown case
};

//
class gmsh_neighbor {
	public:
		// types of neighbors will be determined in GMSH physical groups, including periodic boundaries
		int id;
		int type;

		// Default constructor
		gmsh_neighbor() {
			id = -1;
			type = 0;
		}

		gmsh_neighbor & operator = ( const gmsh_neighbor rhs_neigh ){

			id   = rhs_neigh.id;
			type = rhs_neigh.type;

			return *this;
		}
};

// Contains information about neighboring cell
struct t_gmsh_neighbor {
	enum neigh_type type;		// Type of cell
	int  id;					// neighbor id
	int  sm;					// Submesh  id
	int  proc;					// MPI rank
	int  bound;

	t_gmsh_neighbor & operator = ( const t_gmsh_neighbor rhs_neigh ){

		id    = rhs_neigh.id;
		type  = rhs_neigh.type;
		sm    = rhs_neigh.sm  ;
		proc  = rhs_neigh.proc;
		bound = rhs_neigh.bound;

		return *this;
	}
};

// Node
struct gmsh_node{
	int id;
	int id_global;
	double xn[3];
	bool is_bnode;

	gmsh_node & operator = ( const gmsh_node rhs_node ){

		id        = rhs_node.id;
		id_global = rhs_node.id_global;
		xn[0]     = rhs_node.xn[0];
		xn[1]     = rhs_node.xn[1];
		xn[2]     = rhs_node.xn[2];
		is_bnode  = rhs_node.is_bnode;

		return *this;
	}
};

// Face
struct gmsh_face{
	int id;
	int id_global;
	double xf[3];

	int bound;
	bool is_bface;
	bool is_periodic;
	struct fastmesh_periodic_face twin_face;

	gmsh_face & operator = ( const gmsh_face rhs_face ){
		id          = rhs_face.id;
		id_global   = rhs_face.id_global;
		xf[0]       = rhs_face.xf[0];
		xf[1]       = rhs_face.xf[1];
		xf[2]       = rhs_face.xf[2];
		bound       = rhs_face.bound;
		is_bface    = rhs_face.is_bface;
		is_periodic = rhs_face.is_periodic;
		twin_face   = rhs_face.twin_face;

		return *this;
	}
};

// Cell
struct gmsh_cell{

	int    id, id_global;
	int    faceCount;
	int	   nodeCount;
	double xc[3];
	bool   is_bcell;

	gmsh_cell & operator = ( const gmsh_cell rhs_cell ){
		id        = rhs_cell.id;
		id_global = rhs_cell.id_global;
		faceCount = rhs_cell.faceCount;
		nodeCount = rhs_cell.nodeCount;
		xc[0]     = rhs_cell.xc[0];
		xc[1]     = rhs_cell.xc[1];
		xc[2]     = rhs_cell.xc[2];
		is_bcell  = rhs_cell.is_bcell;

		return *this;
	}
};

struct gmsh_mesh{

	int id;

	// Mesh center
	double x[3];

	// Elements
	unsigned short elmnt_type;
	int nodes_per_cell;
	int faces_per_cell;
	int nodes_per_face;
	int cells_per_face;

	std::vector<gmsh_node> nodes;
	std::vector<gmsh_face> faces;
	std::vector<gmsh_cell> cells;

	// Hpath vars
	std::vector<int> hpath;
	std::vector<int> deadend_cells;

	// Connectivity
	std::vector<std::vector<int>> face2node;
	std::vector<std::vector<int>> face2cell;
	std::vector<std::vector<int>> cell2face;
	std::vector<std::vector<int>> cell2node;
	std::vector<std::vector<t_gmsh_neighbor>> cell2neigh;

	// MPI connectivity
	std::vector<int> mpi_neighbors;

	// Geometrical parameters
	std::vector<double> cells_vol;
	std::vector<double> faces_area;

	std::vector<double> delta_c;
	std::vector<double> delta_f;
	std::vector<double> dot_prod;
	std::vector<std::vector<double>> v_tang;
	std::vector<std::vector<double>> v_norm;
	std::vector<std::vector<double>> v_neigh;
	std::vector<std::vector<double>> v_sign;
};

struct t_mesh_size{

	int n_nodes, n_faces, n_cells;

	int nodes_per_face;
	int cells_per_face;
	int nodes_per_cell;
	int faces_per_cell;
};

struct probe {
	double x[2];
	vector<int> sm;
	vector<int> nodes;
};

//***************************************************************************************************
// Fastmesh vector
//***************************************************************************************************
template < typename T >
class fm_vector{

	public:
		// Default constructor
		fm_vector(){
			array_is_allocated = false;
			m_array = NULL;
		}

		// Return the array size
		fm_vector( int n_elements ){

			m_array = new T[n_elements];
			array_size = (size_t) n_elements;
			array_is_allocated = true ;
		}

		~fm_vector(){
		}

		// Return the array size
		size_t size(){
			return array_size;
		}

		// Return TRUE if array empty and FALSE if not
		bool empty(){
			return array_size == 0 ? true : false;
		}

		// Resize function: allocate array
		void resize( int n_elements ){

			if( array_is_allocated ){
				std::cout << "fm_vector is already allocated, exiting simulation..." << std::endl;
				MPI_Abort( MPI_COMM_WORLD, 202 );
			}

			m_array = new T[n_elements];

			array_size = (size_t) n_elements;
			array_is_allocated = true;
		}

		T operator [] ( int64_t index ) const {
			return m_array[index];
		}
		T & operator [] ( int64_t index ) {
			return m_array[index];
		}

		void set_array_to_unallocated(){
			array_is_allocated = false;
		}

		void   clear(){
			if(m_array)
			 	delete[] m_array;
			m_array = NULL;
			array_is_allocated = false;
		}

		size_t start(){
 			return m_start; }
		size_t end  (){ 
			return m_end; }

		void set_start( int start ){ m_start = start; }
		void set_end  ( int end   ){ m_end   = end  ; }

	private:
		bool array_is_allocated = false;
		T* m_array;
		size_t array_size;

		size_t m_start;
		size_t m_end;
};

#endif








































