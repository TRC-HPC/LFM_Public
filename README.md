# libFastMesh - Finite Volume CFD Solver

## Directories & Files

libFastMesh files are divided into the following directories:
 - api: contains the header files (.h) for all the classes
 - src: contains the source files (.cpp)
 - info: contains the main function lfm_solve.cpp responsible to parse the input arguments and run libFastMesh
 - dictReaderOF, polyMeshReaderOF, runTimeManagerOF: contains the extension modules for OpenFOAM which are used within libFastMesh

The main classes within libFastMesh are:
- MPI_env (api/mpi_env.h, src/mpi_env.cpp): responsible for the MPI communication, for both one-sided and two-sided communication
- IMesh (api/fastmesh.h, src/mesh_reader.cpp, src/mesh_solver.cpp): Interface representing the complete computational mesh assigned to the processor. Creates instances of ISolver which are responsible to solve a submesh  (boundary / interior).
 Responsible to execture the iterations and utilize MPI_env and ISolver instances to communicate the halo cells between adjacent neighboring processors.
- ISolver (api/cfdv0_solver.h, src/cfd_v0.cpp): Interface representing a submesh of the complete mesh. Include data structures which contains the cells geometry and conservatives.
- polyMeshReaderOF: Responsible to read the polyMesh in OpenFOAM format and extract the geometry and boundary conditions. Enables creation of additional scalar and vector variables.
- dictReaderOF: Allows extraction of arguments from OpenFOAM dictionary files
- runTimeManagerOF: Responsible to read the controlDict parameters referring to I/O frequency and format, and write the scalar/vector fields within the current timestamp directory.
-
The classes are written in templates format to minimize the memory footprint based on executed case - allowing both single/double precision, 2D/3D problems and geometry aware data structures.

## Installation

### Prerequisities

#### Solver requirements
- GNU make
- GNU autotools
- MPI compiler (development)
- OpenFOAM v2112 or higher

#### OpenFOAM environment
The OpenFOAM environment must be initialized prior to LFM building, via:
```[bash]
source ${FOAM_INSTALLDIR}/etc/bashrc
```

#### OpenFOAM utilities
The following OpenFOAM utilities are required by libFastMesh:
 - libdictReaderOF.so: Reads OpenFOAM dictionaries
 - libpolyMeshReaderOF.so: Reads OpenFOAM mesh
 - librunTimeManagerOF.so: Support OpenFOAM runTime capabilities

These utilities are found under runTimeManagerOF, polyMeshReaderOF and dictReaderOF directories, respectively. They can be built using:
```[bash]
cd runTimeManagerOF && wmake && cd ..
cd polyMeshReaderOF && wmake && cd ..
cd dictReaderOF && wmake && cd ..
```

#### Hamiltonian Path Renumbering (Optional)
In order to use the hpath renumbering, go to `${FOAM_INSTALLDIR}/src/renumber/renumberMethods`. 
Then, place the hpathRenumber folder inside it, and modify the Make/files file to include `hpathRenumber/hpathRenumber.C`.
Finally, execute `wmake` inside renumberMethods folder. renumberMesh should now support "hpath" as a renumbering method.
Please note that this method works on 2D meshes only.

### Building LFM
To build LFM, the following commands need to be executed:
```[bash]
autoreconf -i
./configure --prefix=${INSTALL_DIR} CC=mpicc CXX=mpic++ CXXFLAGS=-O3
make -j 8
make install
```

for a debug build, replace the configure command with:
```[bash]
autoreconf -i
./configure --prefix=${INSTALL_DIR} CC=mpicc CXX=mpic++ CXXFLAGS="-g -O0 -Wall -DDEBUG"
make -j 8
make install
```
 
${INSTALL_DIR} is the location where the user wants to place the solver binaries.
For more details on the installation process, please refer to the **INSTALL** document.

After installation, it is required to copy the following library files (found inside lib directory) to ${INSTALL_DIR}/lib directory:


It is also recommended to add the ${INSTALL_DIR}/bin directory to PATH.


## Launching a Simluation

Several example test cases are stored under the 'examples' directory. Consult examples/README.md for a list of test cases.

### Input Files
LFM input files consists of a polyMesh, and dictionary files containing the solver parameters

 - `constant/polyMesh`: This directory must contain the polyMesh in the same format read by OpenFOAM. The boundary file must be with the following restrictions:
	 - Patch type must be one of the following: [wall, patch, empty, processor]
	 - Inlet / outlet patches are set using type 'patch', and a name 'inlet' / 'outlet'
 - `0/U, 0/p, 0/T, 0/alpha`: Scalar and vector fields dictionaries with initial values of the velocity, pressure, temperature and alpha
	 - alpha scalar field contains the cell distance from the outlet boundary, and is used only when sponge layer is applied (Ls > 0)
 - `constant/spongeDict`: dictionary file
	 - *Ls* sponge-layer length
	 - *M* freestream mach number
	 - *gamma* specific heat ratio
	 - *pinf* freestream pressure
	 - *Tinf* freestream temperature
	 - *Uinf_x,Uinf_y,Uinf_z* freestream velocity
	 - *k* sponge layer argument [default 4.34294481903252e-01]
 - `constant/thermophysicalProperties`: dictionary file
	 - *mixture/thermodynamics/Cp* gas Cp value
	 - *mixture/specie/molWeight* molecular weight
	 - *mixture/transport/mu* gas Mu0 value
	 - *mixture/transport/Pr* Prandtl value
 - `system/fvSchemes`: dictionary file
	 - *lfm/solver* solver method 
		 - *0* : Use M1 scheme
		 - *1*: Use M2 scheme
	 - *lfm/dimension* problem dimension (2/3)
	 - *lfm/rkOrder* number of runge-kutta steps [default 5]
	 - *lfm/minmodExists* perform minmod calculate [default false]
	 - *lfm/constantTimeStep* use constant time step [default true]
 - `system/controlDict`: dictionary file
	 - Read all the simulation execution prameters required by runTime class of OpenFOAM
	 - *lfm/commType* control the MPI exchanged data
		 - 0: exchange all boundary cells with each neighbor
		 - 1: exchange only necessary cells with each neighbor
		 - 2: exchange only necessary vars with each neighbor. split the communication to two steps, sending only required data for each step
	 - *lfm/haloCommType* control MPI communication pattern [default 1]
		 - 0: two-sided, blocking communication. Typically used only for overlap benchmarking
		 - 1: two-sided, non-blocking communication. Potentially enables communication-computation overlap
		 - 2: two-sided persistent, non-blocking communication. Used for repetitive communication pattern between neighboring ranks
		 - 3: one-sided (RMA), blocking communication. Typically used only for overlap benchmarking
		 - 4: one-sided (RMA), non-blocking communication. Potentially enables communication-computation overlap
		 - 5: one-sided (RMA) persistent, non-blocking communication. Used for repetitive communication pattern between neighboring ranks
		 - 6: neighborhood collectives, blocking
		 - 7: neighborhood collectives, non-blocking
		 - 8: neighborhood collectives, persistent
	 - *lfm/doublePrecision* use double precision flag [default true]
	 - *lfm/post/haveProbes* enable probes [TBD]
	 - *lfm/post/haveSampling* enable sampling [TBD]
	 - *lfm/post/haveAverage* write pressure average and RMS scalar fields [default false]
	 - *lfm/post/haveForces* write forces on wall boundaries [default false]
 	 - *lfm/post/saveForcesStep* write forces frequency [default 1]
	 - *lfm/post/haveResiduals* write residual for each iteration
	 - *lfm/post/printInfoFreq* write info to terminal frequencty [default 1]
	 - *lfm/post/saveResiduals* write residual scalar field [default false]
	 - *lfm/post/saveRank* write processor rank scalar field
	 - *lfm/post/saveCellOrder* write cell order scalar field
	 - *maxCO* maximum CFL value (variable local time-step)


### Mesh Partitioning
In order to partition the mesh and prepare it for simulation with LFM, you may use 'decomposePar' from OpenFOAM to partition the polyMesh into to the desired number of subdomains.
Information about decomposePar usage may be found in the following link: https://openfoamwiki.net/index.php/DecomposePar

### Mesh Renumber
After the mesh partitioning, you may change the cells order to improve performance using the 'renumberMesh' tool from OpenFOAM. 
Users may set the renumberMethod inside renumberMeshDict, as detailed in the following link: https://openfoamwiki.net/index.php/RenumberMesh.
Please note that renumberMesh must be invoked with '-parallel' and '-overwrite' flags to apply the change for each submesh generated during the mesh partitioning step.
 
### Execution
After the preprocessing stage completes, the simulation is ready for running. The command to begin the simulation in parallel:

```[bash]
mpirun -np NP lfm_solve -p
```
The *'-p'* argument tells LFM to run in parallel based on the partitioned sub-meshes found inside the directories '`processor#rank`' (created during the **decomposePar** step)

For serial execution of LFM, the *-p'* argument should be omitted, in which case LFM will read the mesh from `constant/polyMesh` directory.
  
### Post-processing
The results are saved in the same format as OpenFOAM solvers. Users may reconstruct the solution from all processors using *reconstructPar* tool of OpenFOAM, and examine it using VisIt or Paraview.

The integral force coefficients that are summed over each wall boundary are written to the 'Output' directory.
