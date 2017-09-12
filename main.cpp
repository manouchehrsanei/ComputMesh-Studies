
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"

#include "TPZVTKGeoMesh.h"

using namespace std;

TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el);


int main() 

{ 
    
    // ******************* (Create linear meshes: 1D) ***********************************************************
    
    REAL domain = 1.0;
    long num_el = 10;
    REAL size_el = domain/num_el;
    
    TPZGeoMesh *gmesh_OneDL = CreateOneDLGMesh(num_el, size_el); // function to create the 1D geometric mesh
    

return 0; 
}

// ************************************** Create 1D linear meshes ***************************************

TPZGeoMesh *CreateOneDLGMesh(long num_el, REAL size_el)
{
    TPZGeoMesh * gmesh_OneDL = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("geomesh OneDL"); // geometry name
    
    gmesh_OneDL->SetName(name);
    gmesh_OneDL->SetDimension(geometry_dim);
    
    long num_nodes = num_el + 1; // Number of the nodes
    gmesh_OneDL->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
    
    
    const int physical_id = 1; // Define id for material
    const int bc0 = -1; // Define id for left boundary condition
    const int bc1 = -2; // Define id for right boundary condition
    
    for (long i = 0 ; i < num_nodes; i++)
    {
        const REAL valElem = i * size_el;
        TPZVec <REAL> coord(3,0.);
        coord[0] = valElem;
        gmesh_OneDL->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
        gmesh_OneDL->NodeVec()[i].SetNodeId(i); // The id identification
    }
    
    // Creating linear element and  zero-dimensional boundary element
    TPZVec <long> Linear_topology(2); // Vector of the node index: One-dimensional element
    TPZVec <long> point_topology(1); // Vector of the node index: Zero-dimensional element
    long element_id = 0;
    
    
    // Elements
    
    for (long iel = 0; iel < num_el; iel++)
    {
        const long inod_l = iel;
        const long inod_r = iel + 1;
        Linear_topology[0] = inod_l;
        Linear_topology[1] = inod_r;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneDL);
        
    }
    
    element_id++;
    
    
    // Left boundary condition
    point_topology[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneDL);
    element_id++;
    
    
    // Right boundary condition
    point_topology[0] = num_nodes-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneDL);
    
    
    gmesh_OneDL->BuildConnectivity(); // Construct mesh neighbor connectivity
    
    
    std::ofstream outgmeshOneDL("geomesh_OneDL.txt");
    gmesh_OneDL->Print(outgmeshOneDL);
    
    std::ofstream vtkgmeshOneDL("geomesh_OneDL.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneDL, vtkgmeshOneDL);
    
    
    return gmesh_OneDL;
}
