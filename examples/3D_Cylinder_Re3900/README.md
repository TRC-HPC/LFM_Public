# Case

The case here is the 3D turbulent flow over a cylinder. Unsteady time-marching Large Eddy Simulation is performed.

The domain size is as followed:
  - 38D in radial direction (13D are employed for the non-reflective zone)
  - pi*D in span-wise direction
  
# Case setup

The case consists of a cylinder of unit diameter at the center of the computational domain. The inlet Mach number is 0.2, and the flow properties are such that the Reynolds number is 3900. Three different mesh sizes are provided:
  - Small (S): around 1 million total cells
  - Medium (M): around 12.5 million total cells
  - Large (L): around 80 million total cells  
  
Because of the large size of some of the mesh files, git-lfs was used to store some of the generated polyMesh files.
  
To select communication type (blocking / non-blocking), the 'haloCommType' should be modified (found inside system/controlDict) where:
	0 = two-sided blocking
	1 = two-sided non-blocking
  3 = one-sided (RMA) blocking
  4 = one-sided (RMA) non-blocking
The commType controls how the boundary data is sent and should remain 2

# Run the cases

The run.sh script untars the heavy mesh files, performs the domain decomposition, renumbers of the mesh, and runs the test case. Just run the script as:
```
./run.sh X
```
where X is the desired number of processes.

# Results and reference data

Data are compared with experimental data of [Parnaudeau et. al](https://aip.scitation.org/doi/abs/10.1063/1.2957018)
