//***************************************************************************************************
// Function template declaration
//***************************************************************************************************

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// CFD cells
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//
template void MPI_env::set_halo_comm_type<float , 2, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,2,3>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 2, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,2,4>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 2, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,2,3>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 2, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,2,4>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>>& ghost_hpath );

template void MPI_env::set_halo_comm_type<float , 3, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,3>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,3>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,4>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,4>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 5>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,5>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,5>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 6>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,6>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,6>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 7>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,7>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,7>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 8>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,8>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,8>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<float , 3, 9>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<float ,3,9>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,9>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 3, 3>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,3>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,3>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 3, 4>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,4>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,4>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 3, 5>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,5>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,5>>>& ghost_hpath );								
template void MPI_env::set_halo_comm_type<double, 3, 6>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,6>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,6>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 3, 7>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,7>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,7>>>& ghost_hpath );								
template void MPI_env::set_halo_comm_type<double, 3, 8>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,8>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,8>>>& ghost_hpath );
template void MPI_env::set_halo_comm_type<double, 3, 9>(
										  enum t_mpi_halo_comm halo_comm_type,
										  fm_vector<CFDv0_cell<double,3,9>>&  cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,9>>>& ghost_hpath );

// Send/recv, one-sided, RMA, blocking
template void MPI_env::rget<float , 2, 3>(fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::rget<float , 2, 4>(fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::rget<double, 2, 3>(fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::rget<double, 2, 4>(fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 3>(fm_vector<fm_vector<CFDv0_cell<float ,3,3>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 4>(fm_vector<fm_vector<CFDv0_cell<float ,3,4>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 5>(fm_vector<fm_vector<CFDv0_cell<float ,3,5>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 6>(fm_vector<fm_vector<CFDv0_cell<float ,3,6>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 7>(fm_vector<fm_vector<CFDv0_cell<float ,3,7>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 8>(fm_vector<fm_vector<CFDv0_cell<float ,3,8>>> &ghost_hpath );
template void MPI_env::rget<float , 3, 9>(fm_vector<fm_vector<CFDv0_cell<float ,3,9>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 3>(fm_vector<fm_vector<CFDv0_cell<double,3,3>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 4>(fm_vector<fm_vector<CFDv0_cell<double,3,4>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 5>(fm_vector<fm_vector<CFDv0_cell<double,3,5>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 6>(fm_vector<fm_vector<CFDv0_cell<double,3,6>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 7>(fm_vector<fm_vector<CFDv0_cell<double,3,7>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 8>(fm_vector<fm_vector<CFDv0_cell<double,3,8>>> &ghost_hpath );
template void MPI_env::rget<double, 3, 9>(fm_vector<fm_vector<CFDv0_cell<double,3,9>>> &ghost_hpath );

// Send/recv, one-sided, RMA, non-blocking
template void MPI_env::irget<float , 2, 3>(fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::irget<float , 2, 4>(fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::irget<double, 2, 3>(fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::irget<double, 2, 4>(fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 3>(fm_vector<fm_vector<CFDv0_cell<float ,3,3>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 4>(fm_vector<fm_vector<CFDv0_cell<float ,3,4>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 5>(fm_vector<fm_vector<CFDv0_cell<float ,3,5>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 6>(fm_vector<fm_vector<CFDv0_cell<float ,3,6>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 7>(fm_vector<fm_vector<CFDv0_cell<float ,3,7>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 8>(fm_vector<fm_vector<CFDv0_cell<float ,3,8>>> &ghost_hpath );
template void MPI_env::irget<float , 3, 9>(fm_vector<fm_vector<CFDv0_cell<float ,3,9>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 3>(fm_vector<fm_vector<CFDv0_cell<double,3,3>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 4>(fm_vector<fm_vector<CFDv0_cell<double,3,4>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 5>(fm_vector<fm_vector<CFDv0_cell<double,3,5>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 6>(fm_vector<fm_vector<CFDv0_cell<double,3,6>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 7>(fm_vector<fm_vector<CFDv0_cell<double,3,7>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 8>(fm_vector<fm_vector<CFDv0_cell<double,3,8>>> &ghost_hpath );
template void MPI_env::irget<double, 3, 9>(fm_vector<fm_vector<CFDv0_cell<double,3,9>>> &ghost_hpath );

// Send/recv, two-sided, blocking
template void MPI_env::sendrecv<float , 2, 3>(
										  fm_vector<CFDv0_cell<float ,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 2, 4>(
										  fm_vector<CFDv0_cell<float ,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 2, 3>(
										  fm_vector<CFDv0_cell<double,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 2, 4>(
										  fm_vector<CFDv0_cell<double,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );

template void MPI_env::sendrecv<float , 3, 3>(
										  fm_vector<CFDv0_cell<float ,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,3>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 4>(
										  fm_vector<CFDv0_cell<float ,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,4>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 5>(
										  fm_vector<CFDv0_cell<float ,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,5>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 6>(
										  fm_vector<CFDv0_cell<float ,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,6>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 7>(
										  fm_vector<CFDv0_cell<float ,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,7>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 8>(
										  fm_vector<CFDv0_cell<float ,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,8>>> &ghost_hpath );
template void MPI_env::sendrecv<float , 3, 9>(
										  fm_vector<CFDv0_cell<float ,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,9>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 3>(
										  fm_vector<CFDv0_cell<double,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,3>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 4>(
										  fm_vector<CFDv0_cell<double,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,4>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 5>(
										  fm_vector<CFDv0_cell<double,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,5>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 6>(
										  fm_vector<CFDv0_cell<double,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,6>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 7>(
										  fm_vector<CFDv0_cell<double,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,7>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 8>(
										  fm_vector<CFDv0_cell<double,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,8>>> &ghost_hpath );
template void MPI_env::sendrecv<double, 3, 9>(
										  fm_vector<CFDv0_cell<double,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,9>>> &ghost_hpath );

// Send/recv, two-sided, non-blocking
template void MPI_env::isendrecv<float , 2, 3>(
										  fm_vector<CFDv0_cell<float ,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 2, 4>(
										  fm_vector<CFDv0_cell<float ,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 2, 3>(
										  fm_vector<CFDv0_cell<double,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 2, 4>(
										  fm_vector<CFDv0_cell<double,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );

template void MPI_env::isendrecv<float , 3, 3>(
										  fm_vector<CFDv0_cell<float ,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,3>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 4>(
										  fm_vector<CFDv0_cell<float ,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,4>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 5>(
										  fm_vector<CFDv0_cell<float ,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,5>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 6>(
										  fm_vector<CFDv0_cell<float ,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,6>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 7>(
										  fm_vector<CFDv0_cell<float ,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,7>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 8>(
										  fm_vector<CFDv0_cell<float ,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,8>>> &ghost_hpath );
template void MPI_env::isendrecv<float , 3, 9>(
										  fm_vector<CFDv0_cell<float ,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,9>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 3>(
										  fm_vector<CFDv0_cell<double,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,3>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 4>(
										  fm_vector<CFDv0_cell<double,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,4>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 5>(
										  fm_vector<CFDv0_cell<double,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,5>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 6>(
										  fm_vector<CFDv0_cell<double,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,6>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 7>(
										  fm_vector<CFDv0_cell<double,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,7>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 8>(
										  fm_vector<CFDv0_cell<double,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,8>>> &ghost_hpath );
template void MPI_env::isendrecv<double, 3, 9>(
										  fm_vector<CFDv0_cell<double,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,9>>> &ghost_hpath );



// initializeCompletePersistantSendRecv
template void MPI_env::initializeCompletePersistantSendRecv<float , 2, 3>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,3>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 2, 4>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,2,4>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 2, 3>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,2,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,3>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 2, 4>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,2,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,2,4>>> &ghost_hpath );

template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 3>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,3>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 4>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,4>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 5>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,5>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 6>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,6>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 7>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,7>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 8>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,8>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<float , 3, 9>(const int ineigh, 
										  fm_vector<CFDv0_cell<float ,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<float ,3,9>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 3>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,3>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,3>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 4>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,4>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,4>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 5>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,5>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,5>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 6>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,6>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,6>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 7>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,7>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,7>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 8>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,8>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,8>>> &ghost_hpath );
template void MPI_env::initializeCompletePersistantSendRecv<double, 3, 9>(const int ineigh, 
										  fm_vector<CFDv0_cell<double,3,9>>  &cells_hpath,
								fm_vector<fm_vector<CFDv0_cell<double,3,9>>> &ghost_hpath );

// Packing isend/recv
template void MPI_env::new_isendrecv<float >(int ineigh, int send_size, float  *send_buf, int recv_size, float  *recv_buf, const int nWinIndex);
template void MPI_env::new_isendrecv<double>(int ineigh, int send_size, double *send_buf, int recv_size, double *recv_buf, const int nWinIndex);

// Persistant Initialize
template void MPI_env::initializePersistantSendRecv(int ineigh, int send_size, float *send_buf, int recv_size, float *recv_buf, const int nWinIndex);
template void MPI_env::initializePersistantSendRecv(int ineigh, int send_size, double *send_buf, int recv_size, double *recv_buf, const int nWinIndex);


template void MPI_env::new_irget<float >(int ineigh, int recv_size, float  *recv_buf, int nWindowOffset, int nWinIndex);
template void MPI_env::new_irget<double>(int ineigh, int recv_size, double *recv_buf, int nWindowOffset, int nWinIndex);


