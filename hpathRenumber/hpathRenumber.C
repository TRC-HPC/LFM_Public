/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/
#include <queue>
#include <set>
#include "hpathRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Visit Marker  * * * * * * * * * * * * * * * //

/*
 * In order to implement BFS, we will need a visited array.
 * However, initializing such array every BFS costs O(N),
 * where most BFSs will terminate after a constant number of processed nodes.
 * This class supports O(1) marking nodes as visited and checking whether they're visited,
 * and clearing the class costs O(# of nodes processed).
 * It is implemented as a boolean vector and a stack of nodes to clear.
 */
class Foam::VisitMarker {
private:
    std::vector<bool> bVisArray;
    std::stack<int> nVisStack;
public:
    void init(int N) {
        bVisArray.resize(N);
        bVisArray.assign(N,false);
    }
    VisitMarker(int n) { init(n); }

    bool check(int i) {
        if(i<0||i>=int(bVisArray.size())) {
                Info << "out of bounds error at VisitMarker: i="<<i<<". Returning false." << nl << nl;
                return false;
            }
        return bVisArray[i];
    }

    void mark(int i) {
        if(bVisArray[i]) return;
        bVisArray[i] = true;
        nVisStack.push(i);
    }
    void clear() {
        while(!nVisStack.empty()) {
            bVisArray[nVisStack.top()] = false;
            nVisStack.pop();
        }
    }
};
// ************************************************************************* //

namespace Foam
{
    defineTypeNameAndDebug(hpathRenumber, 0);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        hpathRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hpathRenumber::hpathRenumber(const dictionary& renumberDict)
:
    renumberMethod(renumberDict) {}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::hpathRenumber::renumber
(
    const pointField& points
) const
{
    Info << "Note: Hpath renumber cannot be used with this signature. Only the version with the polymesh works." << nl;
    Info << "Returning identity renumbering." << nl << nl;
    labelList orderList(identity(points.size()));
    return orderList;
}

Foam::labelList Foam::hpathRenumber::renumber
(
    const labelListList& cellCells,
    const pointField& points
) const
{
    Info << "Note: Hpath renumber cannot be used with this signature. Only the version with the polymesh works." << nl;
    Info << "Returning identity renumbering." << nl << nl;
    labelList orderList(identity(points.size()));
    return orderList;
}



Foam::labelList Foam::hpathRenumber::renumber
(
    const polyMesh& mesh,
    const pointField& points
) const
{
    // Thats where the magic happens!!!!
    // Cell Count
    const int nCellCount = mesh.cells().size();
    Info << "******************************************************\n";
    Info << "Start Cell Renumbering: " << nCellCount << nl;
    std::vector<bool> bValidFaceList(mesh.nFaces(), true);
    //Mark empty boundary patches as invalid faces.
    const polyBoundaryMesh& bndMesh = mesh.boundaryMesh();
    for(int nBndType = 0; nBndType < bndMesh.size(); ++nBndType) {
        if(bndMesh[nBndType].type().compare("empty")==0) {
            labelRange range = bndMesh.patchRanges()[nBndType];
            for(int nFaceIdx = range.min(); nFaceIdx <= range.max(); ++nFaceIdx) bValidFaceList[nFaceIdx] = false;
        }
    }
    //Get cells' submesh indices, then find boundary hpath, then interior hpath.
    labelList cellOrder(identity(nCellCount));
    std::vector<int> nBoundaryHpathList, nCellSubmeshIndexList;
    getCellSubmeshIndex(mesh, nCellSubmeshIndexList, bValidFaceList);
    int nCurCellOrderIdx = 0;
    if(getBoundaryHpath(mesh, nCellSubmeshIndexList, bValidFaceList, 0, nBoundaryHpathList)) {
        for(int x : nBoundaryHpathList) cellOrder[nCurCellOrderIdx++] = x;
    } else {
        Info << "Get boundary hpath failed. Using identity renumbering instead." << nl;
        for(int nCellIdx = 0; nCellIdx < mesh.nCells(); ++nCellIdx) {
            if(nCellSubmeshIndexList[nCellIdx] == 0) cellOrder[nCurCellOrderIdx++] = nCellIdx;
        }
    }
    Info << "boundary hpath is OK!" << nl;
    std::vector<int> nInteriorHpathList, nInteriorDeadendList;
    if(getInteriorHpath(mesh,nCellSubmeshIndexList,bValidFaceList,1,nInteriorHpathList, nInteriorDeadendList)) {
        for(int nCellIdx : nInteriorHpathList) 
            cellOrder[nCurCellOrderIdx++] = nCellIdx;
        
        for(int nCellIdx : nInteriorDeadendList)
            cellOrder[nCurCellOrderIdx++] = nCellIdx;
        
    } else {
        Info << "Get interior hpath failed. Using identity renumbering instead." << nl << nl;
        for(int nCellIdx = 0; nCellIdx < mesh.nCells(); ++nCellIdx) {
            if(nCellSubmeshIndexList[nCellIdx] == 1) cellOrder[nCurCellOrderIdx++] = nCellIdx;
        }
    }
    Info << "******************************************************\n";

    return cellOrder;

}


int Foam::hpathRenumber::getNei(const polyMesh& mesh, int nCell, int nFaceIdx) const {
    if(nFaceIdx >= mesh.faceNeighbour().size()) return -1;
    int nOwner = mesh.faceOwner()[nFaceIdx], nNei = mesh.faceNeighbour()[nFaceIdx];
    if(nOwner == nCell) return nNei;
    return nOwner;
}

bool Foam::hpathRenumber::getBoundaryHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nBoundaryIndex, std::vector<int>& resultingHpath) const {
    //The idea for starting cells:
    // - The first cell with a boundary face.
    // - Include any quad cells with more than 1 non-boundary face.
    const int NUM_STARTING_CELLS = 10; 
    std::vector<int> nStartingCellsList(NUM_STARTING_CELLS, -1);
    int nCurStartingCells = 0;
    for(int nCellIdx = 0; nCellIdx < mesh.nCells(); ++nCellIdx) {
        if(nCellSubmeshIndexList[nCellIdx]==nBoundaryIndex) {
            nStartingCellsList[nCurStartingCells++] = nCellIdx;
            break;
        }
    }
    int nFaceNeighbourSize = mesh.faceNeighbour().size();
    //We are looking for quads on the boundary with more than two (boundary) neighbors
    for(int nCellIdx = 0; nCellIdx < mesh.nCells() && nCurStartingCells < NUM_STARTING_CELLS; ++nCellIdx) {
        if(nCellSubmeshIndexList[nCellIdx]!=nBoundaryIndex) continue;
        auto& cellData = mesh.cells()[nCellIdx];
        int validFaces = 0;
        for(int nFaceIdx : cellData) if(bValidFaceList[nFaceIdx]) ++validFaces;
        if(validFaces < 4) continue;
        int nGoodNeighbors = 0;
        for(int nFaceIdx : cellData) {
            if(!bValidFaceList[nFaceIdx]) continue;
            if(nFaceIdx < nFaceNeighbourSize && nCellSubmeshIndexList[getNei(mesh, nCellIdx, nFaceIdx)]==nBoundaryIndex) {
                if(++nGoodNeighbors > 1) break;
            }
        }
        if(nGoodNeighbors > 1) nStartingCellsList[nCurStartingCells++] = nCellIdx;
    }
 
    //Try building an hpath from each one of the starting cells.
    if(nCurStartingCells==0) {
        Info << "Could not find any good starting cells for the boundary hpath." << nl;
        return false;
    }

    //count total amount of boundary cells (required for findBoundaryHpath)
    int nBoundaryCells = 0;
    for(int nSubmeshIndex : nCellSubmeshIndexList) if(nSubmeshIndex == nBoundaryIndex) ++nBoundaryCells;

    //Call findBoundaryHpath
    for(int startCell : nStartingCellsList) {
        if(startCell==-1) break;
        if(findBoundaryHpath(mesh, nCellSubmeshIndexList, bValidFaceList, nBoundaryIndex, startCell, nBoundaryCells, resultingHpath))
            return true;
    }
    
    Info << "Could not find any suitable boundary hpath." << nl;
    return false;
}


int Foam::hpathRenumber::getNextBoundaryHpathCell(const polyMesh& mesh, VisitMarker& visitMarker, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bWalkedCellsList, std::vector<bool>& bValidFaceList, int nStartingCell, int nBoundaryIndex, int nDeadendNeighbour) const {
    //We want to find the neareast yet-to-be-walked boundary node.
    //We will use a BFS using a queue and the VisitMarker.
    //The visitmarker helps us keep track of which nodes where visited during our BFS.

    //Every set of nodes with the same distance will end with a -1.
    //Visiting a -1 means that we've finishes processing all nodes with current distance.
    std::queue<int> BFSQueue;
    visitMarker.clear();
    int nCurCell, nCurNei;
    BFSQueue.push(nStartingCell);
    BFSQueue.push(-1);
    int nCurDistance = 0;
    while(!BFSQueue.empty() && nCurDistance <= 10) {
        nCurCell = BFSQueue.front();
        BFSQueue.pop();
        if(nCurCell==-1) { //distance increases
            nCurDistance++;
            BFSQueue.push(-1);
            continue;
        }
        if(visitMarker.check(nCurCell)) continue;
        visitMarker.mark(nCurCell);
        if(nCurCell != nStartingCell && !bWalkedCellsList[nCurCell] && nCurCell != nDeadendNeighbour && nCellSubmeshIndexList[nCurCell] == nBoundaryIndex) {
            //This cell is a boundary cell
            //it is also not walked and not deadend
            //it is also one of the closest nodes that we can return, so let's return it
            return nCurCell;
        }
        for(int nFaceIdx : mesh.cells()[nCurCell]) {
            nCurNei = getNei(mesh,nCurCell,nFaceIdx);
            if(nCurNei > -1 && !visitMarker.check(nCurNei)) BFSQueue.push(nCurNei);
        }
    }
    //No boundary cell found at all.
    return -1;
}
bool Foam::hpathRenumber::isNeighbourDeadend(const polyMesh& mesh, int nCurCell, int nNeighCell, std::vector<bool>& bWalkedCellsList, std::vector<int>& nCellSubmeshIndexList) const {
    for(int nFaceIdx : mesh.cells()[nNeighCell]) {
        int nCurNei = getNei(mesh, nNeighCell, nFaceIdx);
        //skip no-nei, our cell, and walked cells.
        if(nCurNei < 0 || nCurNei == nCurCell || bWalkedCellsList[nCurNei] || nCellSubmeshIndexList[nCurNei] != nCellSubmeshIndexList[nCurCell]) continue;
        return false;
    }
    return true;
}

int Foam::hpathRenumber::getDeadendNeighbourBnd(const polyMesh& mesh, int nCurCell, int nLastCell, std::vector<bool>& bWalkedCellsList, std::vector<int>& nCellSubmeshIndexList) const {
    int nCellsLeft = mesh.cells()[nCurCell].nFaces();
    int nNeighCell;
    bool bNeighDeadend = false;
    for(int nFaceIdx: mesh.cells()[nCurCell]) {
        nNeighCell = getNei(mesh, nCurCell, nFaceIdx);
        if(nNeighCell < 0 || bWalkedCellsList[nNeighCell] || nCellSubmeshIndexList[nNeighCell] != nCellSubmeshIndexList[nCurCell]) --nCellsLeft;
    }
    //If only one path left, there is no choice anyway
    //If we are deadend, there is also no choice anyway
    if(nCellsLeft <= 1) return -1;
    //Otherwise, check if one of them is deadend
    for(int nFaceIdx : mesh.cells()[nCurCell]) {
        nNeighCell = getNei(mesh, nCurCell, nFaceIdx);
        if(nNeighCell < 0 || bWalkedCellsList[nNeighCell] || nNeighCell == nLastCell) continue;
        if(nCellSubmeshIndexList[nNeighCell] != nCellSubmeshIndexList[nCurCell]) continue;
        bNeighDeadend = isNeighbourDeadend(mesh, nCurCell, nNeighCell, bWalkedCellsList, nCellSubmeshIndexList);
        if(bNeighDeadend) break;
    }
    return bNeighDeadend ? nNeighCell : -1;
}

bool Foam::hpathRenumber::findBoundaryHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nBoundaryIndex, int nStartingCell, int nBoundaryCells, std::vector<int>& resultingHpath) const {
    int nTotalCells = mesh.nCells();
    std::vector<int> nTmpHpathList(nTotalCells,-1);
    std::vector<bool> isDeadendList(nTotalCells,false);
    VisitMarker visitMarker(nTotalCells);
    std::vector<bool> bWalkedCellsList(nTotalCells, false);
    int nLastCell = -1;
    int nFaceNeighbourSize = mesh.faceNeighbour().size();
    int nTmpCell = getNextBoundaryHpathCell(mesh, visitMarker, nCellSubmeshIndexList, bWalkedCellsList, bValidFaceList, nStartingCell, nBoundaryIndex, -1);
    int nValidFaces = 0;
    for(int nFaceIdx : mesh.cells()[nStartingCell]) if(bValidFaceList[nFaceIdx]) ++nValidFaces;
    for(int nFaceIdx : mesh.cells()[nStartingCell]) {
        if(nFaceIdx >= nFaceNeighbourSize) continue;
        int nCurNeighbour = getNei(mesh, nStartingCell, nFaceIdx);
        if(nCellSubmeshIndexList[nCurNeighbour]!=nBoundaryIndex) continue;
        if(nValidFaces==3) {
            if(nCurNeighbour!=nTmpCell) nLastCell = nCurNeighbour;
        } else {
            if(nCellSubmeshIndexList[nCurNeighbour]==nBoundaryIndex) nLastCell = nCurNeighbour;
        }
    }
    //Initialize counters and deadend list 
    int nThisCell;
    int nNextCell = nStartingCell;
    int walkedBoundaryCells = 0;
    int nTotalWalkedCells = 0;
    int nDeadendCount = 0;
    std::vector<int> nTmpDeadendList;
    nTmpDeadendList.reserve(nTotalCells);
    do {
        nThisCell = nNextCell;
        bWalkedCellsList[nThisCell] = true;
        if(nCellSubmeshIndexList[nThisCell] == nBoundaryIndex) walkedBoundaryCells++;
        nTmpHpathList[nTotalWalkedCells++] = nThisCell;

        int nDeadendNeighbour = getDeadendNeighbourBnd(mesh, nThisCell, nLastCell, bWalkedCellsList, nCellSubmeshIndexList);
        if(nDeadendNeighbour > -1) {
            nTmpDeadendList.push_back(nDeadendNeighbour);
            nDeadendCount++;
            isDeadendList[nDeadendNeighbour] = true;
        }
        nNextCell = getNextBoundaryHpathCell(mesh, visitMarker, nCellSubmeshIndexList, bWalkedCellsList, bValidFaceList, nThisCell, nBoundaryIndex, nDeadendNeighbour);
    } while(nNextCell != -1);
    nTmpHpathList.resize(nTotalWalkedCells);

    //Check if boundary cells are missing and if so, add to deadends
    for(int nCellIdx = 0; nCellIdx < nTotalCells; ++nCellIdx) {
        if(nCellSubmeshIndexList[nCellIdx] != nBoundaryIndex) continue;
        if(bWalkedCellsList[nCellIdx] || isDeadendList[nCellIdx]) continue;
        nTmpDeadendList.push_back(nCellIdx);
        nDeadendCount++;
    }
    resultingHpath.clear();
    resultingHpath.reserve(nTmpHpathList.size()+nTmpDeadendList.size());
    
    for(int nCell : nTmpHpathList) resultingHpath.push_back(nCell);
    
    for(int nCell : nTmpDeadendList)
        if(!bWalkedCellsList[nCell]) resultingHpath.push_back(nCell);
    
    //resultingHpath is now initialized correctly.
    int nActualBoundaryCells = std::count(nCellSubmeshIndexList.begin(),nCellSubmeshIndexList.end(),0);
    if(int(resultingHpath.size()) != nActualBoundaryCells) {
        Info << "Error: resulting hpath size is " << int(resultingHpath.size()) << ", actual number of cells " << nActualBoundaryCells << nl;
        return false;
    }
    return true;
}


void Foam::hpathRenumber::getCellSubmeshIndex(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList) const {

    //There are two definitions of boundary - 
    //1. Having a face on the boundary
    //2. Having a point of a face on the boundary.
    //The first part of the code marks cells with face on the boundary. The second part of the code marks the rest.
 
    //currently there are only two submeshes - interior/boundary.
    //boundary are assigned with index 0, interior with 1.
    //Finding boundary cells - the idea:
        // - Given boundary faces (faces with no neighbour), mark their points as boundary.
        // - For each cell, check every face if it has a boundary endpoint.
    

    std::vector<bool> isBoundaryPts(mesh.points().size(),false);
    int nCells = mesh.cells().size();
    //A cell is at index 1 if it is interior, and 0 if it is boundary.
    //Initially, everyone are interior. We will later on mark specific cells as boundary.
    nCellSubmeshIndexList.resize(nCells); nCellSubmeshIndexList.assign(nCells,1);
    //Iterate over boundary faces (faces with no neighbour) and mark their points as boundary.
    //Also, mark the face owner as a boundary cell.
    for(int nFaceIdx = mesh.faceNeighbour().size(); nFaceIdx < mesh.nFaces(); ++nFaceIdx) {
        if(!bValidFaceList[nFaceIdx]) continue;
        for(int nPointIdx : mesh.faces()[nFaceIdx]) isBoundaryPts[nPointIdx] = true;
        nCellSubmeshIndexList[mesh.faceOwner()[nFaceIdx]] = 0;
    }
 
    //Now we can iterate over the rest of the faces
    //and mark both the owner and the neighbour as boundary
    //if one of the points are on the boundary.
    int nInteriorCells = mesh.faceNeighbour().size();
    for(int nFaceIdx = 0; nFaceIdx < nInteriorCells; ++nFaceIdx) {
        for(int nPointIdx : mesh.faces()[nFaceIdx]) {
            if(isBoundaryPts[nPointIdx]) {
                nCellSubmeshIndexList[mesh.faceOwner()[nFaceIdx]] = 0;
                nCellSubmeshIndexList[mesh.faceNeighbour()[nFaceIdx]] = 0;
                break;
            }
        }
    }
    
    //nCellSubmeshIndexList is now initialized correctly.
}
Foam::point Foam::hpathRenumber::getCenter(const polyMesh& mesh) const {
    int n = mesh.nCells();
    double X = 0, Y = 0, Z = 0;
    for(Foam::point point : mesh.cellCentres()) X += point.x(), Y += point.y(), Z += point.z();
    return Foam::point(X/n,Y/n,Z/n);
}

double Foam::hpathRenumber::computeSquaredDistance(const point& P, const point& Q) const {
    double res = 0;
    for(int nCurDim = 0; nCurDim < P.dim; ++nCurDim) {
        res += (P[nCurDim]-Q[nCurDim])*(P[nCurDim]-Q[nCurDim]);
    }
    return res;
}


std::vector<int> Foam::hpathRenumber::getInteriorStartingCells(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex) const { 
    point pMeshCenter = getCenter(mesh);
    std::vector<int> nBoundaryCellsList;
    std::vector<int> bIsBoundaryList(mesh.nCells(), false);
    //boundary are those who have a neighbour from another submesh
    //let's iterate over faces and mark both of them
    //we will later have to make sure they are from our submesh
    for(int nFaceIdx = 0; nFaceIdx < int(mesh.faceNeighbour().size()); ++nFaceIdx) {
        int nCell = mesh.faceOwner()[nFaceIdx], nNei = mesh.faceNeighbour()[nFaceIdx];
        if(nCellSubmeshIndexList[nCell] != nCellSubmeshIndexList[nNei]) {
            bIsBoundaryList[nCell] = bIsBoundaryList[nNei] = true;
        }
    }
    for(int nCellIdx = 0; nCellIdx < int(mesh.nCells()); ++nCellIdx) {
        if(nCellSubmeshIndexList[nCellIdx]!=nSubmeshIndex) bIsBoundaryList[nCellIdx] = false;
        if(bIsBoundaryList[nCellIdx]) {
            nBoundaryCellsList.push_back(nCellIdx);
        }
    }

    // Find irregular cells (cells with only one regular neighbor)
    std::vector<int> nIrregularNeighsList;
    std::vector<bool> bIsIrregularCellNeigh(mesh.nCells(), false);

    int nTmpIrregularNeigh;
    for(int nBndCell : nBoundaryCellsList) {
        int nRegularNeiCnt = 0;
        for(int nFaceIdx : mesh.cells()[nBndCell]) {
            if(!bValidFaceList[nFaceIdx]) continue;
            int nCurNei = getNei(mesh, nBndCell, nFaceIdx);
            if(nCurNei < 0) continue;
            if(nCellSubmeshIndexList[nCurNei] != nSubmeshIndex) continue;
            nRegularNeiCnt++;
            nTmpIrregularNeigh = nCurNei;
        }
        if(nRegularNeiCnt<2) {
            nIrregularNeighsList.push_back(nTmpIrregularNeigh);
            bIsIrregularCellNeigh[nTmpIrregularNeigh] = true;
        }
    }

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Find cell next to and not part of boundary hpath
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
    std::vector<int> nTmpStartingCellsList;
    nTmpStartingCellsList.reserve(mesh.nCells());
    // Loop over boundary cells
	for(int nBndCell : nBoundaryCellsList) {
        if(nBndCell < 0) continue;
		// Condition 0: no neighbor can be an irregular cell's neighbor
        bool bValidStartingCell = true;
        for(int nFaceIdx : mesh.cells()[nBndCell]) {
            if(!bValidFaceList[nFaceIdx]) continue;
            int nCurNei = getNei(mesh, nBndCell, nFaceIdx);
            if(nCurNei < 0  || nCellSubmeshIndexList[nCurNei]!=nSubmeshIndex) continue;
            if(bIsIrregularCellNeigh[nCurNei]) {
                bValidStartingCell = false;
                break;
            }
        }
        if(!bValidStartingCell) continue;

		// Condition 1: the two neighbors must have at least two regular cells other than the current
		// cell
        int nValidNeiCnt = 0;
        for(int nFaceIdx : mesh.cells()[nBndCell]) {
            if(!bValidFaceList[nFaceIdx]) continue;
            int nCurNei = getNei(mesh, nBndCell, nFaceIdx);
            if(nCurNei < 0  || nCellSubmeshIndexList[nCurNei]!=nSubmeshIndex) continue;
            int nRegularNeisCnt = 0;
            for(int nNeiFaceIdx : mesh.cells()[nCurNei]) {
                if(!bValidFaceList[nNeiFaceIdx]) continue;
                int nNeiNei = getNei(mesh,nCurNei, nNeiFaceIdx);
                if(nNeiNei < 0 || nCellSubmeshIndexList[nNeiNei]!=nSubmeshIndex) continue;
                if(nNeiNei == nBndCell) continue;
                nRegularNeisCnt++;
            }
            if(nRegularNeisCnt > 1) nValidNeiCnt++;
        }
        if(nValidNeiCnt < 2) continue;
		// Condition 2: the two neighbors must have at least one interior cell (neighbour)
        nValidNeiCnt = 0;
        for(int nFaceIdx : mesh.cells()[nBndCell]) {
            int nCurNei = getNei(mesh, nBndCell, nFaceIdx);
            if(nCurNei < 0  || nCellSubmeshIndexList[nCurNei]!=nSubmeshIndex) continue;
            bool bHasInteriorCellNei = false;
            for(int nNeiFaceIdx : mesh.cells()[nCurNei]) {
                int nNeiNei = getNei(mesh, nCurNei, nNeiFaceIdx);
                if(nNeiNei < -1) continue;
                if(bIsBoundaryList[nNeiNei]) {
                    bHasInteriorCellNei = true;
                    break;
                }
            }
            if(bHasInteriorCellNei) ++nValidNeiCnt;
        }
        if(nValidNeiCnt >= 2) nTmpStartingCellsList.push_back(nBndCell);
    }

	// If no starting cells, take the farthest cell as a starting point
	if ( nTmpStartingCellsList.empty() ){
		double dMaxDist = computeSquaredDistance(mesh.cellCentres()[nBoundaryCellsList[0]], pMeshCenter);
        int nFarthestCell = nBoundaryCellsList[0];
		for(int nBndCell : nBoundaryCellsList) {
            if(nBndCell < 0) continue;
            double dCurDist = computeSquaredDistance(mesh.cellCentres()[nBndCell], pMeshCenter);
            if(dCurDist > dMaxDist) {
                dMaxDist = dCurDist;
                nFarthestCell = nBndCell;
            }
        }
        nTmpStartingCellsList.push_back(nFarthestCell);
	}
    
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Order starting cells randomly
	// 	- Hpaths beginning with starting cells close to each other often yield similar results
	// 	- make sure consecutive starting cells have a low probability of being next to each other
	// 	  thus finding an Hpath quicker (if one exists)
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Compute the angle of each starting cell with respect to the x-axis
    std::vector<std::vector<float>> fCellAnglesList(nTmpStartingCellsList.size(), std::vector<float>(2,0));
    for(int nCellIdx = 0; nCellIdx < int(nTmpStartingCellsList.size()); ++nCellIdx) {
        int nCell = nTmpStartingCellsList[nCellIdx];
        fCellAnglesList[nCellIdx][0] = atan2( mesh.cellCentres()[nCell].y()-pMeshCenter.y(), mesh.cellCentres()[nCell].x()-pMeshCenter.x());
        fCellAnglesList[nCellIdx][1] = nCellIdx;
    }
    sort(fCellAnglesList.begin(),fCellAnglesList.end());
    // std::random_shuffle(fCellAnglesList.begin(),fCellAnglesList.end());
	// std::random_shuffle(nTmpStartingCellsList.begin(), nTmpStartingCellsList.end());
    int nTmpStartingCellsCnt = nTmpStartingCellsList.size();
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Compute ratio of cells: assume the algorithm should go over 2M cells in the worst case scenario
	// e.g.: if mesh has 100K cells, we should have a max of 20 Hpaths attempts (feel free to change this number)
	// =======> this limits the time spent iterating on Hpaths. As the starting cells are spread out
	//          around the mesh's boundary, it should not affect significantly the best performance
	//          of the Hpath
	const int TRESHOLD = 2000000;
	int MAX_HPATH_ITER = TRESHOLD / mesh.nCells();

	int SELECT_RATIO   = MAX_HPATH_ITER > nTmpStartingCellsCnt ? 1 : nTmpStartingCellsCnt / MAX_HPATH_ITER;

	// Only select a handful of starting cells, spread out evenly along the boundary
	int nSelectedStartingCellsCnt = nTmpStartingCellsCnt / SELECT_RATIO + 1;
    std::vector<int> nFinalStartingCellsList;
    nFinalStartingCellsList.reserve(nSelectedStartingCellsCnt);

    for(int nListIdx = 0; nListIdx < nTmpStartingCellsCnt; nListIdx++) {
        if(nListIdx % SELECT_RATIO != 0) continue;
        //We could really just loop over the multiples of SELECT_RATIO. 
        //Replicating this for good luck of course.
        int nCellIdx = fCellAnglesList[nListIdx][1];
        nFinalStartingCellsList.push_back(nTmpStartingCellsList[nCellIdx]);
    }

	if(nFinalStartingCellsList.empty()){
		Info << "Could not find any good interior starting cells." << nl;
	}
	return nFinalStartingCellsList;
}


//In the interior algorithm, get (one) deadend cell neighbour (if it exists)
int  Foam::hpathRenumber::getInteriorDeadendNeighbour(const polyMesh& mesh, int nThisCell, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bWalkedCellsList, std::vector<bool>& bValidFaceList) const {
    int nNeiCell;
    bool bNeiDeadend = false;
    int nCurSubmesh = nCellSubmeshIndexList[nThisCell];
    //Count how many options do we have
	int nCellsLeft = 0;
    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        if(!bValidFaceList[nFaceIdx]) continue;
        nNeiCell = getNei(mesh, nThisCell, nFaceIdx);
        if(nNeiCell < 0 || nCellSubmeshIndexList[nNeiCell] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNeiCell]) continue;
        ++nCellsLeft;
    }
    //If no cells left, in particular no deadend nei
    //If one cell left, no choice anyway
    if(nCellsLeft <= 1) return -1;

	//If more than one path is available, check if one of them is a deadend

    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        if(!bValidFaceList[nFaceIdx]) continue;
        nNeiCell = getNei(mesh, nThisCell, nFaceIdx);
        if(nNeiCell < 0 || nCellSubmeshIndexList[nNeiCell] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNeiCell]) continue;
        bNeiDeadend = isNeighbourDeadend(mesh, nThisCell, nNeiCell, bWalkedCellsList, nCellSubmeshIndexList);
        if(bNeiDeadend) break;
    }
    return bNeiDeadend ? nNeiCell : -1;
}

//In the interior standard algo, get the next hpath cell
int  Foam::hpathRenumber::getInteriorNextHpathCell(const polyMesh& mesh, int nThisCell, std::vector<int>& nCellSubmeshIndexList, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, point& pMeshCenter) const {
    
    int nCurSubmesh = nCellSubmeshIndexList[nThisCell];
	int nNextCell = -1;
    double dDistance, dMaxDistance = 0;

	// Compute the number of available cells
    std::vector<int> nFreeCellsList;
    nFreeCellsList.reserve(mesh.cells()[nThisCell].size());
    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        int nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNei]) continue;
        nFreeCellsList.push_back(nNei);
    }

	// If only 1 cell is free, no further work required
    if(nFreeCellsList.size() == 1) return nFreeCellsList[0];
    
	// If no cells are available, we've reached the end of the Hpath
	if(nFreeCellsList.empty()) return -1;
	
	// If multiple free cells...
	for( int nFreeCell : nFreeCellsList ){
		// Skip deadend neighbor: will be added after Hpath is done
		if( nFreeCell == nDeadendNeighbour ) continue;
        dDistance = sqrt(computeSquaredDistance(mesh.cellCentres()[nFreeCell], pMeshCenter));
        
		if( dMaxDistance < dDistance ){
            dMaxDistance = dDistance;
            nNextCell = nFreeCell;
		}
	}
	return nNextCell;
}

//In the interior standard algo, get the next interior hpath cell
int  Foam::hpathRenumber::getInteriorNextHpathCellBnd(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, int nThisCell, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, bool& bUpdateHpathAlgo) const {
    //Currently implemented as nested loops. Could be changed to distance-bound BFS.
    int nNei;
    int nCurSubmesh = nCellSubmeshIndexList[nThisCell];
	// Priority 1: boundary cells
    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        //skip out of bounds neighs, walked cells, and deadend neighbour
        nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNei] || nNei == nDeadendNeighbour) continue;
        if(bIsCellBnd[nNei] && !bWalkedCellsList[nNei])
            return nNei;
    }
    
	// Priority 2: boundary's neighbors
    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNei] || nNei == nDeadendNeighbour) continue;
        if(bCellHasBndNei[nNei] && !bWalkedCellsList[nNei])
            return nNei;
    }

    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNei] || nNei == nDeadendNeighbour) continue;
        for(int nNeiFaceIdx : mesh.cells()[nNei]) {
            int nNeiNei = getNei(mesh, nNei, nNeiFaceIdx);
            if(nNeiNei < 0 || nCellSubmeshIndexList[nNeiNei] != nCurSubmesh) continue;
            if(bWalkedCellsList[nNeiNei]) continue;
            if(bCellHasBndNei[nNeiNei])
                return nNeiNei;
        }
    }

	// If we've reached this point, the boundary cells and their neighbors are probably all walked
    for(int nFaceIdx : mesh.cells()[nThisCell]) {
        nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) continue;
        if(bWalkedCellsList[nNei] || nNei == nDeadendNeighbour) continue;
        bUpdateHpathAlgo = true;
        return nNei;

    }
	return -1;
}

int  Foam::hpathRenumber::getInteriorNextHpathCellSpiral(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, int nThisCell, int nDeadendNeighbour, std::vector<bool>& bWalkedCellsList, std::vector<int>& nHpathIdxList, bool& bUpdateHpathAlgo) const {
    
    int nCurSubmesh = nCellSubmeshIndexList[nThisCell];
    std::vector<int> nBndCellsList;
    std::vector<int> nBndNeighsCellsList;
    std::vector<int> nNormalCellsList;
    std::vector<bool> bBadCells(mesh.cells()[nThisCell].size(),false);

    int nNei;

	// Classify all neighbors
    for(int nFaceIdxIdx  = 0; nFaceIdxIdx < mesh.cells()[nThisCell].size(); ++nFaceIdxIdx) {
        int nFaceIdx = mesh.cells()[nThisCell][nFaceIdxIdx];
        nNei = getNei(mesh, nThisCell, nFaceIdx);
        if(nNei < 0 || nCellSubmeshIndexList[nNei] != nCurSubmesh) {
            bBadCells[nFaceIdxIdx] = true;
            continue;
        }
        if(bWalkedCellsList[nNei] || nNei == nDeadendNeighbour) {
            bBadCells[nFaceIdxIdx] = true;
            continue;
        }
        if(bIsCellBnd[nNei]) {
            nBndCellsList.push_back(nNei);
            continue;
        } else if(bCellHasBndNei[nNei]) {
            nBndNeighsCellsList.push_back(nNei);
            continue;
        } else {
            nNormalCellsList.push_back(nNei);
        }
    }
	// Priority 1: boundary cells
	if (!nBndCellsList.empty()) return nBndCellsList[0]; 
	
	// Priority 2: boundary cell neighbors. For each candidate:
	// 		- give priority to cells with boundary cells neighbors which weren't walked
	// 		- if multiple candidates found, pick either one (for now, haven't found a case where this happens yet)
	// 		- if no candidates, pick cell with oldest walked boundary cell
	if (!nBndNeighsCellsList.empty()) {
		int nWalkedIndex = int(mesh.nCells())*2;  
		int nBestCandidate = -1;
        for(int nNei : nBndNeighsCellsList) {
            //Candidate cell
            for(int nFaceIdx : mesh.cells()[nNei]) {
                int nNeiNei = getNei(mesh, nNei, nFaceIdx);
                if(nNeiNei < 0 || nCellSubmeshIndexList[nNeiNei] != nCurSubmesh) continue;
                if(!bIsCellBnd[nNeiNei]) continue;
                if(nWalkedIndex > nHpathIdxList[nNeiNei]) {
                    nWalkedIndex = nHpathIdxList[nNeiNei];
                    nBestCandidate = nNei;
                }
            }
        }
		if(nBestCandidate  == -1) {
            Info << "Spiral algorithm error - boundary cell neighbors fails." << nl;
        } else return nBestCandidate;
	}
    if(nNormalCellsList.empty()) return -1;
    if(nNormalCellsList.size() == 1) return nNormalCellsList[0];

	// Priority 3: for each candidate, look for one that has a walked neighbor
	int nWalkedIndex = int(mesh.nCells())*2;
	int nBestCandidate = -1;
    for(int nNei : nNormalCellsList) {
        for(int nNeiFaceIdx : mesh.cells()[nNei]) {
            int nNeiNei = getNei(mesh, nNei, nNeiFaceIdx);
            if(nNeiNei < 0) continue;
            if(!bWalkedCellsList[nNeiNei]) continue;
            if(nNeiNei == nThisCell) continue;
            if(nWalkedIndex > nHpathIdxList[nNeiNei]) {
                nWalkedIndex =nHpathIdxList[nNeiNei];
                nBestCandidate = nNei;
            }
        }
    }
	if (nBestCandidate != -1) return nBestCandidate;

	// Priority 4: for each candidate, look for one that has a walked neighbor
	nWalkedIndex = int(mesh.nCells())*2;
	nBestCandidate = -1;
    for(int nNei : nNormalCellsList) {
        for(int nNeiFaceIdx : mesh.cells()[nNei]) {
            int nNeiNei = getNei(mesh, nNei, nNeiFaceIdx);
            if(nNeiNei < 0) continue;
            if(nNeiNei == nThisCell) continue;
            for(int nNeiNeiFaceIdx : mesh.cells()[nNeiNei]) {
                int nNei3 = getNei(mesh, nNeiNei, nNeiNeiFaceIdx);
                if(nNei3 < 0) continue;
                if(!bWalkedCellsList[nNei3]) continue;
				// Ignore cells that were very recently walked
                if(nHpathIdxList[nThisCell]-nHpathIdxList[nNei3] < 3) continue;
                
                if(nWalkedIndex > nHpathIdxList[nNei3]) {
                    nWalkedIndex = nHpathIdxList[nNei3];
                    nBestCandidate = nNei;
                }
            }
        }
    }
	if (nBestCandidate > -1)
		return nBestCandidate;

	// Priority 4(.2): for each candidate, look for one that has a walked neighbor
	nWalkedIndex = int(mesh.nCells())*2;
	nBestCandidate = -1;
    for(int nNei : nNormalCellsList) {
        for(int nNeiFaceIdx : mesh.cells()[nNei]) {
            int nNeiNei = getNei(mesh, nNei, nNeiFaceIdx);
            if(nNeiNei < 0) continue;
            for(int nNeiNeiFaceIdx : mesh.cells()[nNeiNei]) {
                int nNei3 = getNei(mesh, nNeiNei, nNeiNeiFaceIdx);
                if(nNei3 < 0) continue;
                if(!bWalkedCellsList[nNei3]) continue;
                if(nWalkedIndex > nHpathIdxList[nNei3]) {
                    nWalkedIndex = nHpathIdxList[nNei3];
                    nBestCandidate = nNei;
                }
            }
        }
    }
	if (nBestCandidate == -1) return nNormalCellsList[0];
	return nBestCandidate;
}
bool Foam::hpathRenumber::getInteriorNextHpathCellSwitch(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bIsCellBnd, std::vector<bool>& bCellHasBndNei, std::vector<bool>& bWalkedCellsList, point& pMeshCenter, int nThisCell, int nDeadendNeighbour, bool& bUpdateHpathAlgo, std::vector<int>& nHpathIdxList, HpathAlgoType eAlgoType, int& nNextCell) const {
    bool bNowSwitched = false;
    switch(eAlgoType) {
        case ALGO_STANDARD:
				if( bUpdateHpathAlgo ){
                    nNextCell = getInteriorNextHpathCell(mesh, nThisCell, nCellSubmeshIndexList, nDeadendNeighbour, bWalkedCellsList, pMeshCenter);
				}else{
                    nNextCell = getInteriorNextHpathCellBnd(mesh, nCellSubmeshIndexList, bIsCellBnd, bCellHasBndNei, nThisCell, nDeadendNeighbour, bWalkedCellsList, bUpdateHpathAlgo);
					if( bUpdateHpathAlgo ) bNowSwitched = true;
				}
            break;
        case ALGO_TRUE_SPIRAL:
            nNextCell = getInteriorNextHpathCellSpiral(mesh, nCellSubmeshIndexList, bIsCellBnd, bCellHasBndNei, nThisCell, nDeadendNeighbour, bWalkedCellsList, nHpathIdxList, bUpdateHpathAlgo);	
            break;
        default:
            Info << "Error: No such algorithm type exists("<< eAlgoType <<")." << nl;
            break;
    }
    return bNowSwitched;
}



void Foam::hpathRenumber::findInteriorHpathNoDeadend(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex, int nStartingCell, HpathAlgoType eAlgoType, int& nTotalHpathCells, std::vector<int>& nHpathList, std::vector<int>& nDeadendCellsList) const {
    // Find all boundary cells and their regular neighbors
    std::vector<bool> bIsCellBnd(mesh.nCells(), false);
    std::vector<bool> bCellHasBndNei(mesh.nCells(), false);
    Foam::point pMeshCenter = getCenter(mesh);
    int nTotalInteriorCells = std::count(nCellSubmeshIndexList.begin(),nCellSubmeshIndexList.end(), nSubmeshIndex);
    for(int nFaceIdx = 0; nFaceIdx < mesh.faceNeighbour().size(); ++nFaceIdx) {
        int nCell = mesh.faceOwner()[nFaceIdx], nNei = mesh.faceNeighbour()[nFaceIdx];
        if(nCellSubmeshIndexList[nCell]!=nCellSubmeshIndexList[nNei]) {
            if(nCellSubmeshIndexList[nCell]==nSubmeshIndex) bIsCellBnd[nCell] = true;
            if(nCellSubmeshIndexList[nNei]==nSubmeshIndex) bIsCellBnd[nCell] = true;
        }
    }
    for(int nFaceIdx = 0; nFaceIdx < mesh.faceNeighbour().size(); ++nFaceIdx) {
        int nCell = mesh.faceOwner()[nFaceIdx], nNei = mesh.faceNeighbour()[nFaceIdx];
        if(nCellSubmeshIndexList[nCell]!=nCellSubmeshIndexList[nNei]) continue;
        if(nCellSubmeshIndexList[nCell]!=nSubmeshIndex) continue;
        if(bIsCellBnd[nCell]) bCellHasBndNei[nNei] = true;
        if(bIsCellBnd[nNei]) bCellHasBndNei[nCell] = true;
    }

	// Walk Hamiltonian path
    int nTmpHpathIdx = 0;
    int nThisCell;
    int nNextCell = nStartingCell;
    bool bUpdateHpathAlgorithm = false; //This is some info we pass on to the
                                        //find next cell algorithm
                                        //and the algorithm updates it, too.
    std::vector<int> nTmpHpathList(mesh.nCells(),-1);
    std::vector<bool> bWalkedCellsList(mesh.nCells(), false);
    std::vector<int> nTmpDeadendList;
    nTmpDeadendList.reserve(mesh.nCells());

    int nNonHpathCellsCnt = 0;

	// Spiral algorithm
    std::vector<int> nHpathIdxList(mesh.nCells(), -1);
    std::vector<int> nHpathNei(mesh.nCells(), -1);
    std::vector<int> nHpathNeiNei(mesh.nCells(), -1);
	do{
        //Move to next cell and update corresponding arrays
        nThisCell = nNextCell;
        nHpathIdxList[nThisCell] = nTmpHpathIdx;
        nTmpHpathList[nTmpHpathIdx++] = nThisCell;
        bWalkedCellsList[nThisCell] = true;
		// Update hpath
		// If one of the neighboring cell is a deadend, split this_cell and deadend_cell in two
		// Only store information for now, the mesh and hpath will be updated at the end
		int nDeadendNeighbour = getInteriorDeadendNeighbour(mesh, nThisCell, nCellSubmeshIndexList, bWalkedCellsList, bValidFaceList);

		if( nDeadendNeighbour > -1 ){
            nNonHpathCellsCnt++;
            nTmpDeadendList.push_back(nDeadendNeighbour);
		}
        bool bNowSwitched = getInteriorNextHpathCellSwitch(mesh,nCellSubmeshIndexList,bIsCellBnd, bCellHasBndNei, bWalkedCellsList, pMeshCenter, nThisCell, nDeadendNeighbour, bUpdateHpathAlgorithm, nHpathIdxList, eAlgoType, nNextCell);
		if(bNowSwitched) {
            /*Until now we were looking for boundary cells, now we're looking for interior. */
        }
		//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
		// If path fails prematurely, most likely reason is long strip of deadend cells
		if( nNextCell == -1 && nTmpHpathIdx + nNonHpathCellsCnt < nTotalInteriorCells ){
			// Walk back Hpath out of deadend and change direction
            bool bFoundExit = false;
            int nWalkedBackCellsCnt = 1;
            std::vector<bool> bWalkedCellsRestartList = bWalkedCellsList;
            int nRestartNextCell = -1;

			// 1- Starting from last walked cell, walk back one cell
			// 2- set last walked cell as a deadend neighbor
			// 3- look for the next cell considering deadend neighbor
			// 4- If another path is found, take it
			//    else, go back to 1
			do{
				if( nTmpHpathIdx - nWalkedBackCellsCnt <= 0 ) break;
                int nDeadendNeighbour = nTmpHpathList[nTmpHpathIdx-nWalkedBackCellsCnt];
                int nRestartThisCell = nTmpHpathList[nTmpHpathIdx-nWalkedBackCellsCnt-1];

                if(nRestartThisCell == -1) break;

                // nRestartNextCell = getInteriorNextHpathCellSwitch(mesh, nCellSubmeshIndexList, bIsCellBnd, bCellHasBndNei, bWalkedCellsList, pMeshCenter, nRestartThisCell, nDeadendNeighbour, bUpdateHpathAlgorithm, nHpathIdxList, eAlgoType, nRestartNextCell);
                nRestartNextCell = getInteriorNextHpathCell(mesh, nRestartThisCell, nCellSubmeshIndexList, nDeadendNeighbour, bWalkedCellsList, pMeshCenter);

                if(nRestartNextCell != -1)
                    bFoundExit = true;
                else
                    ++nWalkedBackCellsCnt;
                
                if(nWalkedBackCellsCnt == 60 && eAlgoType == ALGO_STANDARD) break;
                if(nWalkedBackCellsCnt == 60 && eAlgoType == ALGO_TRUE_SPIRAL) break;

			}while( !bFoundExit );

			if( nRestartNextCell == -1 ){
				break;
			}
            for(int nCellIdx = 1; nCellIdx < nWalkedBackCellsCnt - 1; nCellIdx++) {
                int nNewDeadendCell = nTmpHpathList[nTmpHpathIdx-nCellIdx];
                bWalkedCellsList[nNewDeadendCell] = false;
                nTmpDeadendList.push_back(nNewDeadendCell);
                nHpathIdxList[nNewDeadendCell] = -1;
                nTmpHpathList[nTmpHpathIdx-nCellIdx] = -1;
                ++nNonHpathCellsCnt;
            }
            nNextCell = nRestartNextCell;
            nTmpHpathIdx -= nWalkedBackCellsCnt;
		}
        if(nNextCell != -1) {
            bWalkedCellsList[nNextCell] = true;
        }
	} while( nNextCell != -1 );

    nTmpHpathList.resize(nTmpHpathIdx);
    nHpathList = nTmpHpathList;
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// Set deadend cells:
	// 		- most of deadends were (hopefully) recorder along with the Hpath algo
	// 		- it's possible that large chunks of the mesh were left out of the Hpath ==> record
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    std::vector<bool> bIsDeadend(mesh.nCells(), true);
    nDeadendCellsList.clear();
    nDeadendCellsList.reserve(mesh.nCells());

	// Define the deadend cells as cells not walked by Hpath
    for(int nHpathIdx = 0; nHpathIdx < int(nHpathList.size()); ++nHpathIdx) 
        bIsDeadend[nHpathList[nHpathIdx]] = false;
	// Walk the Hpath, and record deadends along the way
    int nDeadendCnt = 0;
    std::vector<bool> bRecordedDeadends(mesh.nCells(), false);

    for(int nHpathCell : nHpathList) {
        for(int nFaceIdx : mesh.cells()[nHpathCell]) {
            int nCurNei = getNei(mesh, nHpathCell, nFaceIdx);
            if(nCurNei < 0 || nCellSubmeshIndexList[nCurNei] != nSubmeshIndex) continue;
            if(bIsDeadend[nCurNei] && !bRecordedDeadends[nCurNei]) {
                nDeadendCellsList.push_back(nCurNei);
                bRecordedDeadends[nCurNei] = true;
                nDeadendCnt++;
            }
        }
    }
    // Check if there are deadend cells not yet recorded (i.e. with no Hpath cells neighbor )
    for(int nCellIdx = 0; nCellIdx < mesh.nCells(); ++nCellIdx) {
        if(nCellSubmeshIndexList[nCellIdx] == nSubmeshIndex && bIsDeadend[nCellIdx] && !bRecordedDeadends[nCellIdx]) {
            nDeadendCellsList.push_back(nCellIdx);
            bRecordedDeadends[nCellIdx] = true;
            nDeadendCnt++;
        }
    }
    nDeadendCellsList.resize(nDeadendCnt);
    nTotalHpathCells = nTmpHpathIdx;
    if(nTotalInteriorCells == int(nHpathList.size() + nDeadendCellsList.size())) {
        return;
    }
    Info << "Error: hpath size ("<< int(nHpathList.size()) <<") + deadend cells size ("<< nDeadendCellsList.size() <<") != total number of interior cells(" << nTotalInteriorCells << ")." << nl;
    Info << "Shuffling interior" << nl;
    std::random_shuffle(nHpathList.begin(),nHpathList.end());
    std::random_shuffle(nDeadendCellsList.begin(),nDeadendCellsList.end());

}


bool Foam::hpathRenumber::getInteriorHpath(const polyMesh& mesh, std::vector<int>& nCellSubmeshIndexList, std::vector<bool>& bValidFaceList, int nSubmeshIndex, std::vector<int>& nResultHpathList, std::vector<int>& nResultDeadendCellsList) const {

    std::vector<bool> bWalkedCellsList(mesh.nCells(), false);
    std::vector<int> nStartingCellsList = getInteriorStartingCells(mesh, nCellSubmeshIndexList, bValidFaceList, nSubmeshIndex);
    std::vector<int> nDeadendCellsList;
    nDeadendCellsList.reserve(mesh.nCells());
	
    std::vector<int> nHpathList(mesh.nCells());

    //keep track of the best hpath to recover it later.
    int nTotalHpathCells;
    int nBestTotalHpathCells = -1;
    int nBestHpathTrial = -1;
    HpathAlgoType eBestHpathAlgo = ALGO_STANDARD;
    double dBestHpathPrct = -1;
    
    //count the number of interior cells to make sure
    //later that hpath+deadend=interior
    int nTotalInteriorCells = std::count(nCellSubmeshIndexList.begin(),nCellSubmeshIndexList.end(),nSubmeshIndex);
    for(int nAlgoType = 0; nAlgoType <= 1; nAlgoType++) {
        HpathAlgoType eAlgoType = HpathAlgoType(nAlgoType);   
        for(int nTrialIndex = 0; nTrialIndex < int(nStartingCellsList.size()); ++nTrialIndex) {
            findInteriorHpathNoDeadend(mesh, nCellSubmeshIndexList, bValidFaceList, nSubmeshIndex, nStartingCellsList[nTrialIndex],  eAlgoType ,nTotalHpathCells, nHpathList, nDeadendCellsList);
            
            if(nTotalInteriorCells != int(nHpathList.size() + nDeadendCellsList.size()))
                continue;
            
            if(nBestTotalHpathCells < nTotalHpathCells) {
                nBestTotalHpathCells = nTotalHpathCells;
                nBestHpathTrial = nTrialIndex;
                eBestHpathAlgo = eAlgoType;
                dBestHpathPrct = ((nTotalHpathCells) / double(nTotalInteriorCells)) * 100;
            }
        }
    }
    Info << "Best hpath percent is " << dBestHpathPrct << nl << nl;
    // Info << "best hpath: prct="<<dBestHpathPrct <<",algo="<<(eBestHpathAlgo == ALGO_STANDARD ? "ALGO_STANDARD" : "ALGO_TRUE_SPIRAL")<<",start="<<nStartingCellsList[nBestHpathTrial] << nl << nl;
	// Recover the hpath
    if(nBestHpathTrial < 0) return false;

    //Reconstruct best hpath
    findInteriorHpathNoDeadend(mesh, nCellSubmeshIndexList, bValidFaceList, nSubmeshIndex, nStartingCellsList[nBestHpathTrial],  eBestHpathAlgo ,nTotalHpathCells, nHpathList, nDeadendCellsList);
    

    nResultHpathList = nHpathList;
    nResultDeadendCellsList = nDeadendCellsList;
    return true;
}



// ************************************************************************* //
