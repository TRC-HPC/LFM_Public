#ifndef CFDV0_H
#define CFDV0_H

#include <memory>
#include <functional>
#include <array>
#include <sstream>
#include "elements.h"
#include "cfdv0_elements.h"
#include "mpi_env.h"
#include "api/meshReader.h"
#include "api/inputReader.h"

#define INTERP_LINEAR(weight, a, b) ((weight) * (a) + (1.0 - (weight)) * (b))

using namespace std;

class ISolver
{
	// Public Variables
	public:
		// Mesh Reader
		IMeshReader* m_pMeshReader;

		// Submesh Id
		int m_nSubmeshIndex;

		// Case Type
		case_type case_t;
		std::string case_name;

		// Boundary Conditions
		vector<int> m_nWallBCList, m_nInletBCList, m_nOutletBCList, m_nFarfieldBCList;

		// Ghost boundary cells
		std::vector<std::vector<array<int, 2>>> boundaries;

		// Mapping from rank ID to local MPI neighbor ID
		std::vector<int> rank2local;

		// Input Reader
		CInputReader* m_pInput;
		
	// Public Virtual Methods
	public:
		// Virtual Methods
		virtual void prepare_for_timestep() = 0;		
		virtual void prepare_for_RKstep( int rk_step, int nSolver ) = 0;
		virtual void calc_gradients      ( MPI_env &mpi_envi ) = 0;
		virtual void calc_gradients_M2AUSM      ( MPI_env &mpi_envi) = 0;
		virtual void calc_VIS            ( MPI_env &mpi_env ) = 0;
		virtual void calc_VIS_Smagorinsky( MPI_env &mpi_env ) = 0;
		virtual void one_rk_step_M1      ( int , double , MPI_env &, double * ) {}
		virtual void one_rk_step_M1      ( int , float , MPI_env &, float * ) {}
		virtual void one_rk_step_M2      ( int , double , MPI_env &, double * ) {}
		virtual void one_rk_step_M2      ( int , float , MPI_env &, float *) {}
		virtual void one_rk_step_M2AUSM  ( int , double , MPI_env &, double * ) {}
		virtual void one_rk_step_M2AUSM  ( int , float , MPI_env &, float *) {}
		virtual void exchange_ghost_cells( MPI_env &mpi_env ) = 0;
		virtual void gather_info         ( MPI_env &, struct t_domain_info<double, 2> & ) {}
		virtual void gather_info         ( MPI_env &, struct t_domain_info<double, 3> & ) {}
		virtual void gather_info         ( MPI_env &, struct t_domain_info<float, 2> & ) {}
		virtual void gather_info         ( MPI_env &, struct t_domain_info<float, 3> & ) {}

		// New Output Methods
		virtual void updateAverageField(const int nTimeStep) = 0;
		virtual void updateSolutionResidual() = 0;
		virtual void updateSolutionPrimitives() = 0;
		virtual void updateSolutionBlendFactor() = 0;

		virtual void prePostProc( MPI_env &mpi_env, bool haveProbes, bool haveSampling, bool haveAverage, bool haveForces, vector<probe> &x_probes, vector<probe> &x_samples, int) = 0;
		virtual void postProcAverage( int time_step ) = 0;
		virtual void postProcAverageOutput( gmsh_mesh *, string , float , int  ) {}
		virtual void postProcAverageOutput( gmsh_mesh *, string , double , int  ) {}
		virtual void postProcForces( MPI_env &, float , MPI_Datatype  ) {}
		virtual void postProcForces( MPI_env &, double , MPI_Datatype  ) {}
		virtual void postProcForcesOutput( MPI_env &mpi_env, ofstream &fout, size_t itype ) = 0;
		
		virtual void set_boundary_conditions( ) = 0;

		// gmsh
		virtual void allocate       ( gmsh_mesh *pSubMesh, std::vector<void*> &cells_ptr ) = 0;
		virtual void deallocate     () = 0;
		virtual void allocate_ghost_cells( MPI_env &mpi_env,
								   gmsh_mesh *pSubMesh,
								   std::vector<std::vector<void*>> &ghost_ptr,
								   std::vector<std::vector<void*>> &ghost_ptr_bnd ) = 0;
		virtual void assign_pointers( gmsh_mesh *pSubMesh,
							  MPI_env &mpi_env,
							  std::vector<std::vector<void*>> cells_ptr,
							  std::vector<std::vector<void*>> ghost_ptr,
							  std::vector<std::vector<void*>> ghost_bnd ) = 0;
		virtual void init_params    ( gmsh_mesh *pSubMesh ) = 0;
		virtual void initialize     ( MPI_env &mpi_env, gmsh_mesh *pSubMesh, int nSolver) = 0;
		virtual void initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *pSubMesh) = 0;
		virtual void reorder_faces() = 0;

		virtual double compute_cfl(double  ) {return 0;}
		virtual float compute_cfl(float  ) {return 0;}
		virtual double compute_dt(double  ) {return 0;}
		virtual float compute_dt(float  ) {return 0;}
		
		virtual void mpi_communication(MPI_env &mpi_env, int comm_step) = 0;
		virtual void mpi_wait(MPI_env &mpi_env, int comm_step) = 0;
};

//***************************************************************************************************
//***************************************************************************************************
//
//***************************************************************************************************
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
class CFDv0_solver : public ISolver
{
	public:
		// Runge-Kutta variables
		PRECISION ti, tf, dt;

		vector<PRECISION> m_dAk;
		vector<PRECISION> m_dBk;

		// Cells
		std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> cells_cfd;

		// Dummy cell: used for NEIGHBOR_NONE access
		std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>> dummy_cell;
		t_solution_vars<PRECISION, DIM_CNT> dummy_var;

		std::vector<std::vector<PRECISION>> xc;

		// Ghost cells
		std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_mpi;
		std::vector<std::vector<CFDv0_cell<PRECISION, DIM_CNT, FACE_CNT>>> ghost_bnd;

		// Constants
		const PRECISION ZERO  =  0.0;
		const PRECISION HALF  =  0.5;
		const PRECISION ONE   =  1.0;
		const PRECISION M_ONE = -1.0;

		PRECISION m_dPr, m_dRe, m_dGamma, m_dL, m_dM, m_dGammaMinusOne, m_dPr_inv, m_dGamma_m_one_inv, m_dAoA;
		PRECISION m_dMu0, m_dT0, m_dS;

		PRECISION k, m_dCp, m_dCv;
		const PRECISION m_dRuniversal = 8.31447 ;		// Universal gas constant. TODO: one header for all constants
		PRECISION m_dMolWeight ;						// Molar Weight of the air in g/mol.
		PRECISION m_dRgas, m_dRgas_inv ;					// Specific gas constant.

		// Free stream coefficients
		PRECISION m_dpInf, m_dTInf, m_drhoInf, m_dUInf[DIM_CNT], m_dEInf;

		// For Post-Process Averaging
		vector<PRECISION> m_dpAVG, m_dpRMS;

		// For Post-Process Forces
		vector<PRECISION> m_dtimeVec;
		vector<vector<array<PRECISION, DIM_CNT>>> m_dFpre;
		vector<vector<array<PRECISION, DIM_CNT>>> m_dFvis;

		// new MPI protocol
		std::vector<std::vector<int>> local_cells_to_send;
		std::vector<std::vector<int>> neigh_cells_to_recv;

		//std::vector<PRECISION> send_buf;
		std::vector<int> m_nSendBufOffsetList, m_nRecvBufOffsetList;
		std::vector<int> m_nSendBufCountList, m_nRecvBufCountList;
		std::vector<int> m_nRecvWindowOffsetList;
		std::vector<PRECISION> m_SendBufList;
		std::vector<PRECISION> m_RecvBufList;

		std::vector<int> m_nSendBufSolvarOffsetList, m_nRecvBufSolvarOffsetList;
		std::vector<int> m_nSendBufSolvarCountList, m_nRecvBufSolvarCountList;
		std::vector<int> m_nSendBufViscousOffsetList, m_nRecvBufViscousOffsetList;
		std::vector<int> m_nSendBufViscousCountList, m_nRecvBufViscousCountList;
		std::vector<int> m_nRecvSolvarsWindowOffsetList, m_nRecvViscousWindowOffsetList;
		std::vector<PRECISION> m_SendBufSolvarsList, m_SendBufViscousList;
		std::vector<PRECISION> m_RecvBufSolvarsList, m_RecvBufViscousList;

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Functions
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		virtual void prepare_for_timestep() override;
		virtual void prepare_for_RKstep( int rk_step, int nSolver ) override;
		virtual void calc_gradients      ( MPI_env &mpi_env) override;
		virtual void calc_gradients_M2AUSM( MPI_env &mpi_env) override;
		virtual void calc_VIS            ( MPI_env &mpi_env ) override;
		virtual void calc_VIS_Smagorinsky( MPI_env &mpi_env ) override;
		virtual void one_rk_step_M1      ( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES ) override;
		virtual void one_rk_step_M2      ( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES ) override;
		virtual void one_rk_step_M2AUSM  ( int rk_step, PRECISION dt, MPI_env &mpi_env, PRECISION *RES ) override;
		// M2 functions for AUSM dissipation
		PRECISION    p5Pos(PRECISION M,PRECISION alpha);
		PRECISION    p5Neg(PRECISION M,PRECISION alpha);
		virtual void exchange_ghost_cells( MPI_env &mpi_env ) override;
		virtual void gather_info         ( MPI_env &mpi_env, struct t_domain_info<PRECISION, DIM_CNT> &domain_info );

		// New Output Methods
		virtual void updateAverageField(const int nTimeStep);
		virtual void updateSolutionResidual();
		virtual void updateSolutionPrimitives();
		virtual void updateSolutionBlendFactor();

		virtual void prePostProc( MPI_env &mpi_env, bool haveProbes, bool haveSampling, bool haveAverage, bool haveForces, vector<probe> &x_probes, vector<probe> &x_samples, int ) override;
		virtual void postProcAverage( int time_step ) override;
		virtual void postProcAverageOutput( gmsh_mesh *mesh, string filename, PRECISION time, int time_step ) override;
		virtual void postProcForces( MPI_env &mpi_env, PRECISION time, MPI_Datatype mpi_precision ) override;
		virtual void postProcForcesOutput( MPI_env &mpi_env, ofstream &fout, size_t itype ) override;
		
		virtual void set_boundary_conditions( ) override;

		// gmsh
		virtual void allocate       ( gmsh_mesh *pSubMesh, std::vector<void*> &cells_ptr ) override;
		virtual void deallocate     () override;
		virtual void allocate_ghost_cells( MPI_env &mpi_env,
								   gmsh_mesh *pSubMesh,
								   std::vector<std::vector<void*>> &ghost_ptr,
								   std::vector<std::vector<void*>> &ghost_ptr_bnd ) override;
		virtual void assign_pointers( gmsh_mesh *pSubMesh,
							  MPI_env &mpi_env,
							  std::vector<std::vector<void*>> cells_ptr,
							  std::vector<std::vector<void*>> ghost_ptr,
							  std::vector<std::vector<void*>> ghost_bnd ) override;
		virtual void init_params    ( gmsh_mesh *pSubMesh ) override;
		virtual void initialize     ( MPI_env &mpi_env, gmsh_mesh *pSubMesh, int nSolver ) override;
		virtual void initialize_mpi_env( MPI_env &mpi_env, gmsh_mesh *pSubMesh) override;
		virtual void reorder_faces() override;

		virtual PRECISION compute_cfl( PRECISION dt ) override;
		virtual PRECISION compute_dt( PRECISION cflMax ) override;
		
		virtual void mpi_communication(MPI_env &mpi_env, int comm_step) override;
		virtual void mpi_wait(MPI_env &mpi_env, int comm_step) override;

		// New packed MPI protocol
		void new_exchange_ghost    (MPI_env &mpi_env);
		void exchange_ghost_solvars(MPI_env &mpi_env);
		void exchange_ghost_viscous(MPI_env &mpi_env);
		void unpack_mpi_data(MPI_env &mpi_env);
		void unpack_mpi_data_solvars(MPI_env &mpi_env);
		void unpack_mpi_data_viscous(MPI_env &mpi_env);
	private:

		MPI_Datatype type_cfdcell;

		void init_mpi_types( MPI_env &mpi_env );
		void set_init_conditions( );

		void init_boundary_conditions( );

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline PRECISION compute_Rpsi( PRECISION q_sol[DIM_CNT+2] ){

			PRECISION rho_inv  = ONE / q_sol[0];
			PRECISION rhoU_sqr = q_sol[1] * q_sol[1];
			for( unsigned idim = 1; idim < DIM_CNT; idim++ ){
				rhoU_sqr += q_sol[idim+1] * q_sol[idim+1];
			}

			return m_dGammaMinusOne * (q_sol[DIM_CNT+1] - 0.5 * rhoU_sqr * rho_inv) * rho_inv;
		}

		inline PRECISION dot_product( PRECISION lvec[DIM_CNT], PRECISION rvec[DIM_CNT] ){

			PRECISION dot_prod = lvec[0] * rvec[0];
			for( unsigned idim=1; idim < DIM_CNT; idim++ ){
				dot_prod += lvec[idim] * rvec[idim];
			}
			return dot_prod;
		}

		inline PRECISION vector_mag( PRECISION vec[DIM_CNT] ){

			PRECISION vec_mag = vec[0] * vec[0];
			for( unsigned idim=1; idim < DIM_CNT; idim++ ){
				vec_mag += vec[idim] * vec[idim];
			}
			return sqrt( vec_mag );
		}

		inline void init_array( PRECISION array[DIM_CNT] ){
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				array[idim] = 0.;
			}
		}
		inline void init_array( PRECISION array[DIM_CNT][DIM_CNT] ){
			for( unsigned idim=0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim=0; jdim < DIM_CNT; jdim++ ){
					array[idim][jdim] = 0.;
				}
			}
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Linear interpolation scheme
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		inline PRECISION interp_linear( PRECISION weight, PRECISION cell_phi, PRECISION adjc_phi ){
			return weight * cell_phi + (ONE - weight) * adjc_phi;
		}

		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// Minmod interpolation scheme
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

		// Get sign of variable a
		inline PRECISION sign(PRECISION a) {
			return a < 0.0 ? -1.0 : 1.0 ;
		}

		// Compute "r"
		inline PRECISION calc_r( PRECISION phiP, PRECISION phiN, PRECISION phiGrad[DIM_CNT], PRECISION d[DIM_CNT] ){

			PRECISION gradf  = phiN - phiP + 1.0e-30;
			PRECISION gradcf = ZERO;

			for (unsigned idim=0; idim<DIM_CNT; idim++) {
				gradcf += d[idim] * phiGrad[idim];
			}

			if (abs(gradcf) >= 1000.0 * abs(gradf)) {
				return 2.0 * 1000.0 * sign(gradcf) * sign(gradf) - ONE;
			} else {
				return 2.0 * (gradcf / gradf) - ONE;
			}
		}

		// Compute "r" for Vector
		inline PRECISION calc_rV( PRECISION phiP[DIM_CNT], PRECISION phiN[DIM_CNT], PRECISION phiGrad[DIM_CNT][DIM_CNT], PRECISION d[DIM_CNT] ){

			PRECISION gradfV[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				gradfV[idim] = phiN[idim] - phiP[idim];
			}
			PRECISION gradf = 1.0e-30;
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				gradf += gradfV[idim] * gradfV[idim];
			}
			
			PRECISION gradcfV[DIM_CNT];
			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				for( unsigned jdim = 0; jdim < DIM_CNT; jdim++ ){
					gradcfV[idim] = d[idim] * phiGrad[jdim][idim];
				}
			}

			PRECISION gradcf = ZERO;
			for (unsigned idim=0; idim<DIM_CNT; idim++) {
				gradcf += gradfV[idim] * gradcfV[idim];
			}

			if (abs(gradcf) >= 1000.0 * abs(gradf)) {
				return 2.0 * 1000.0 * sign(gradcf) * sign(gradf) - ONE;
			} else {
				return 2.0 * (gradcf / gradf) - ONE;
			}
		}

		// Minmod interpolation
		inline PRECISION interp_minmod( PRECISION cell_phi, PRECISION adjc_phi,
										PRECISION grad_phi[DIM_CNT], PRECISION d[DIM_CNT],
										PRECISION weight_linear, PRECISION flux ){

			PRECISION r       = calc_r(cell_phi, adjc_phi, grad_phi, d);
			PRECISION limiter = max(min(r, ONE), ZERO);
			PRECISION weight  = limiter * weight_linear + ( ONE - limiter ) * (ONE + flux) * HALF;
			//PRECISION weight = limiter * weight_linear + ( ONE - limiter );
			//weight = HALF * (ONE - flux) + flux * weight;

			return weight * cell_phi + (ONE - weight) * adjc_phi;
		}

		// Minmod Vector interpolation
		inline void interp_minmodV( PRECISION cell_phi[DIM_CNT], PRECISION adjc_phi[DIM_CNT],
										 PRECISION grad_phi[DIM_CNT][DIM_CNT], PRECISION d[DIM_CNT],
										 PRECISION weight_linear, PRECISION flux,
										 PRECISION* &rhoU ){

			PRECISION r       = calc_rV(cell_phi, adjc_phi, grad_phi, d);
			PRECISION limiter = max(min(r, ONE), ZERO);
			PRECISION weight  = limiter * weight_linear + ( ONE - limiter ) * (ONE + flux) * HALF;
			//PRECISION weight = limiter * weight_linear + ( ONE - limiter );
			//weight = HALF * (ONE - flux) + flux * weight;

			for( unsigned idim = 0; idim < DIM_CNT; idim++ ){
				rhoU[idim] = weight * cell_phi[idim] + (ONE - weight) * adjc_phi[idim];
			}
			//return weight * cell_phi + (ONE - weight) * adjc_phi;
		}
};

#ifdef SINGLE_PRECISION
template class CFDv0_solver<float ,2, 3>;
template class CFDv0_solver<float ,2, 4>;
template class CFDv0_solver<float ,3, 3>;
template class CFDv0_solver<float ,3, 4>;
template class CFDv0_solver<float ,3, 5>;
template class CFDv0_solver<float ,3, 6>;
#endif

template class CFDv0_solver<double ,2, 3>;
template class CFDv0_solver<double ,2, 4>;
template class CFDv0_solver<double ,3, 3>;
template class CFDv0_solver<double ,3, 4>;
template class CFDv0_solver<double ,3, 5>;
template class CFDv0_solver<double ,3, 6>;

#endif
