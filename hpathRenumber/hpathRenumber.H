/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::hpathRenumber

Description
    Renumber the cells using a Hamiltonian path.


SourceFiles
    hpathRenumber.C

\*---------------------------------------------------------------------------*/

#ifndef hpathRenumber_H
#define hpathRenumber_H

#include "renumberMethod.H"




namespace Foam
{
class VisitMarker;
/*---------------------------------------------------------------------------*\
                      Class hpathRenumber Declaration
\*---------------------------------------------------------------------------*/

enum HpathAlgoType {
    ALGO_STANDARD = 0,
    ALGO_TRUE_SPIRAL = 1
};

class hpathRenumber
:
    public renumberMethod
{
    // Private data

    // Private Member Functions

        //- No copy construct
        hpathRenumber(const hpathRenumber&) = delete;

        //- No copy assignment
        void operator=(const hpathRenumber&) = delete;

public:


    //- Runtime type information
    TypeName("hpath");

    // Constructors

    //- Construct given the renumber dictionary
    hpathRenumber(const dictionary& renumberDict);


    //- Destructor
    virtual ~hpathRenumber() = default;


    // Member Functions

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  This is only defined for geometric renumberMethods.
        virtual labelList renumber(const pointField&) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  Use the mesh connectivity (if needed)
        virtual labelList renumber
        (
            const polyMesh& mesh,
            const pointField& cc
        ) const;

        //- Return the order in which cells need to be visited, i.e.
        //  from ordered back to original cell label.
        //  The connectivity is equal to mesh.cellCells() except
        //  - the connections are across coupled patches
        virtual labelList renumber
        (
            const labelListList& cellCells,
            const pointField& cc
        ) const;
    protected:

        //Given a cell and a face index, compute its neighbour.
        // - If facing the boundary, returns -1.
        // - Otherwise, returns FaceOwner/FaceNeighbour[nFaceIdx], the one that's different
        //   from nCell.
        int getNei(const polyMesh& mesh, int nCell, int nFaceIdx) const;

        //Check if neighbour is deadend.
        bool isNeighbourDeadend(const polyMesh& mesh, int nCurCell, int nNeighCell, std::vector<bool>& bWalkedCellsList, std::vector<int>& nCellSubmeshIndexList) const;


        //Divide the mesh to submeshes and return for each cell its index.
        //Currently assigns each cell whether it's a boundary cell(0) or an interior cell(1).
        void getCellSubmeshIndex(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList) const;

        //This function computes possible starting cells and checks them with
        //findBoundaryHpath. Returns whether successful or not.
        bool getBoundaryHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nBoundaryIndex, std::vector<int>& resultingHpath) const;
        
        //Try to find a boundary hpath starting from that cell, using a bounded-distance BFS.
        //and store it in the resulting hpath. returns whether successful or not.
        bool findBoundaryHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nBoundaryCells, int nBoundaryIndex, int nStartingCell, std::vector<int>& resultingHpath) const;

        //Find next cell on the boundary hpath.
        int getNextBoundaryHpathCell(const polyMesh& mesh, VisitMarker& visitMarker, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bWalkedCellsList, std::vector<bool>& bValidFaceList, int nStartingCell, int nBoundaryIndex, int nDeadendNeighbour) const;

        //Get a deadend neighbour on boundary hpath (if exists) (currently looking at one deadend neighbour).
        int getDeadendNeighbourBnd(const polyMesh& mesh, int nStartingCell, int nLastCell, std::vector<bool>& bWalkedCellsList, std::vector<int>& nCellSubmeshIndexList) const; 

	    //Compute the center of the mesh (using average of cell centers)
        point getCenter(const polyMesh& mesh) const;

	    //Compute squared euclidean distance
        double computeSquaredDistance(const point& P, const point& Q) const;

        //Get possible starting interior cells
        std::vector<int> getInteriorStartingCells(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex) const;        

        //In the interior algorithm, get (one) deadend cell neighbour (if it exists)
        int getInteriorDeadendNeighbour(const polyMesh& mesh, int nThisCell, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bWalkedCellsList, std::vector<bool>& bValidFaceList) const; 

        //In the interior standard algo, get the next hpath cell
        int getInteriorNextHpathCell(const polyMesh& mesh, int nThisCell, std::vector<int>& nCellSubmeshIndexList, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, point& pMeshCenter) const;

        //In the interior standard algo, get the next interior hpath cell
        int getInteriorNextHpathCellBnd(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, int nThisCell, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, bool& bUpdateHpathAlgo) const;

        //In the interior spiral algo, get the next hpath cell
        int getInteriorNextHpathCellSpiral(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, int nThisCell, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, std::vector<int>& nHpathIdxList, bool& bUpdateHpathAlgo) const;

	    //A wrapper function for GetInteriorNextHpathCell(./Bnd/Spiral) using the enum AlgoType and the boolean bUpdateHpathAlgo.
        bool getInteriorNextHpathCellSwitch(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, std::vector<bool>& bWalkedCellsList, point& pMeshCenter, int nThisCell, int nDeadendNeighbour, bool& bUpdateHpathAlgo, std::vector<int>& nHpathIdxList, HpathAlgoType eAlgoType, int& nNextCell) const;

	    //Tries to find hpath given a starting cell.
        void findInteriorHpathNoDeadend(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex, int nStartingCell, HpathAlgoType eAlgoType, int& nTotalHpathCells, std::vector<int>& nHpathList, std::vector<int>& nDeadendCellsList) const;

        //Compute the interior hpath of a specific submesh
        bool getInteriorHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex, std::vector<int>& nResultHpathList, std::vector<int>& nResultDeadendCellsList) const; 

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
