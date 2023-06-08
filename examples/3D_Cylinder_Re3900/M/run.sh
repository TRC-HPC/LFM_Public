# Untar data
cd 0
tar -xvzf alpha.tgz
rm alpha.tgz
cd ..
cd constant/polyMesh
tar -xvzf cellZones.tgz.tgz
tar -xvzf cellZones.tgz
rm cellZones.tgz.tgz
rm cellZones.tgz
tar -xvzf faces.tgz.tgz
tar -xvzf faces.tgz
rm faces.tgz.tgz
rm faces.tgz
tar -xvzf neighbour.tgz.tgz
tar -xvzf neighbour.tgz
rm neighbour.tgz.tgz
rm neighbour.tgz
tar -xvzf owner.tgz.tgz
tar -xvzf owner.tgz
rm owner.tgz.tgz
rm owner.tgz
tar -xvzf points.tgz.tgz
tar -xvzf points.tgz
rm points.tgz.tgz
rm points.tgz
cd ../..
# Decompose Par
sed -i "s/numberOfSubdomains .*;/numberOfSubdomains $1;/g" system/decomposeParDict
decomposePar -force
# Renumber Mesh
mpirun -np $1 renumberMesh -dict constant/renumberMeshDict -overwrite -parallel
# Run case
mpirun -np $1 caafoam-m2-l -parallel
