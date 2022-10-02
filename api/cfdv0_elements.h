#ifndef CFDV0_ELEMENTS_H
#define CFDV0_ELEMENTS_H


//***************************************************************************************************
// Contains the solution variables stored at the cells center
// 		- separated from the rest for better neighbor access
// 		- only variables require by neighbors: fetch this array instead of the whole cell
// Struct size:
// 		- 2D = 96 bytes
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT >
struct t_solution_vars{

	// Flow variables

	//PRECISION rhs  [DIM_CNT+2];		// RHS      vector at t_{n}
	PRECISION q_old[DIM_CNT+2];			// Solution vector at t_{n}
	PRECISION delta_q[DIM_CNT+2];		// Temporary solution vector at t_{n}, remove

	// Gradients for r calculations interpolation
	PRECISION rho_grad [DIM_CNT];
	PRECISION rhoU_grad[DIM_CNT][DIM_CNT];
	PRECISION rhoE_grad[DIM_CNT];
	PRECISION Rpsi_grad[DIM_CNT];
	PRECISION c_grad   [DIM_CNT];

	// Adding gradients for M2 AUSM dissipation
	PRECISION p_grad[DIM_CNT];
	PRECISION U_grad[DIM_CNT][DIM_CNT];

	//// Artificial Diffusion
	//PRECISION D2U[DIM_CNT+2];

	// Viscous Part
	PRECISION dudx[DIM_CNT][DIM_CNT];
	PRECISION dTdx[DIM_CNT];

	PRECISION tauMC[DIM_CNT][DIM_CNT];
	PRECISION sigmaU[DIM_CNT];

	PRECISION vol_inv;								// Cell volume

	PRECISION RES[DIM_CNT+2];						// Cell's residual

	// Debugging purposes
	unsigned int	id;
	unsigned short	sm_id;
	bool			is_ghost;
};

//***************************************************************************************************
// Cell class for RK4-5 scheme in Caafoam
// 	- tri   = 4 + 96 + 3*8 + 3*2*8 + 3*8 + 3*8 + 3*8 + 4 ====> 246, rounded to 256
// 	- quads = 4 + 96 + 4*8 + 4*2*8 + 4*8 + 4*8 + 4*8 + 4 ====> 296, rounded to 320
//***************************************************************************************************
template < typename PRECISION, unsigned DIM_CNT, unsigned FACE_CNT >
struct CFDv0_cell{

	//unsigned int id;							// Local cell ID
	t_solution_vars<PRECISION, DIM_CNT> vars;
	t_solution_vars<PRECISION, DIM_CNT> *neighs[FACE_CNT];

	// Geomerical parameters

	// Required
	PRECISION S       [FACE_CNT][DIM_CNT];		// Face surface vector
	PRECISION d       [FACE_CNT][DIM_CNT];		// Vector btw cell center to neighbor's center
	PRECISION weight_linear[FACE_CNT];			// Coefficients for interpolation
	PRECISION sponge_sigma;						// Coefficients for sponge layer

	int m_nCellIndex;
	int m_nFaceCount;
};

template < typename PRECISION, unsigned DIM_CNT >
struct t_domain_info{
	PRECISION rho_min, rho_max;
	PRECISION u_min[DIM_CNT], u_max[DIM_CNT];
	PRECISION p_min, p_max;
};

#endif

