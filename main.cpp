
#include <iostream>
#include <cmath>
#include <fstream>

#include "pzreal.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"

#include "pzgeopoint.h"
#include "TPZGeoLinear.h"

#include "TPZVTKGeoMesh.h"

#include "pzl2projection.h"

#include "pzinterpolationspace.h"

#include "pzgeotriangle.h"
#include "pzgeoquad.h"

#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "pzshapelinear.h"


// ************************************** (Create meshes: 1D) *****************************************************

TPZGeoMesh *CreateOneDGMesh(long num_el, REAL size_el);

TPZCompMesh * CreateComputationalOneDMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeOneD(TPZCompMesh * cmesh, int order);

// ********************************** (Create meshes: 2DTri) *****************************************************

TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly);

TPZCompMesh * CreateComputationalTwoDTriMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeTwoDTri(TPZCompMesh * cmesh, int order);

// -----------------------------------(Create meshes: 2DQuad) ---------------------------------------------

TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly);

TPZCompMesh * CreateComputationalTwoDQuadMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeTwoDQuad(TPZCompMesh * cmesh, int order);

// ********************************** (Create meshes: 3DTetra) *************************************************

TPZGeoMesh *CreateThreeDTetraGMesh(long nnodestetra, REAL Lx, REAL Ly, REAL Lz);

TPZCompMesh * CreateComputationalThreeDTetraMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeThreeDTetra(TPZCompMesh * cmesh, int order);

// -----------------------------------(Create meshes: 3DPyram) ---------------------------------------------

TPZGeoMesh *CreateThreeDPyramGMesh(long nnodestetra, REAL Lx, REAL Ly, REAL Lz);

TPZCompMesh * CreateComputationalThreeDPyramMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeThreeDPyram(TPZCompMesh * cmesh, int order);

// -----------------------------------(Create meshes: 3DPrism) ---------------------------------------------

TPZGeoMesh *CreateThreeDPrismGMesh(long nnodestetra, REAL Lx, REAL Ly, REAL Lz);

TPZCompMesh * CreateComputationalThreeDPrismMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeThreeDPrism(TPZCompMesh * cmesh, int order);

// -----------------------------------(Create meshes: 3DHexa) ---------------------------------------------

TPZGeoMesh *CreateThreeDHexaGMesh(long nnodeshexa, REAL Lx, REAL Ly, REAL Lz);

TPZCompMesh * CreateComputationalThreeDHexaMesh(TPZGeoMesh * geometry, int order);

void ComputationShapeThreeDHexa(TPZCompMesh * cmesh, int order);

// -----------------------------------(Create meshes: 1D) ---------------------------------------------


TPZCompMesh *CreateComputationOneDMesh(TPZGeoMesh * gmesh_OneD,int orderOneD);



// **************************************** (main) ***********************************************************


int main()

{
    
    
    REAL Lx = 1.0; // length of domain in x direction
    REAL Ly = 1.0; // length of domain in y direction
    REAL Lz = 1.0; // length of domain in z direction

    
    // ************************************** (Create meshes: 1D) *****************************************************
    
    REAL domain = 1.0;
    long num_el = 1;
    REAL size_el = domain/num_el;
    
    TPZGeoMesh * gmesh_OneD = CreateOneDGMesh(num_el, size_el); // function to create the 1D geometric mesh
//
//    int order1D = 4; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_OneD = CreateComputationalOneDMesh(gmesh_OneD,order1D);
//
//    ComputationShapeOneD(cmesh_OneD, order1D);
//    
//    
//    // ********************************** (Create meshes: 2DTri) *****************************************************
//    
//
//    long nnodes = 3; // Number of the nodes
//    
//    TPZGeoMesh *gmesh_TwoDTri = CreateTwoDTriGMesh(nnodes, Lx, Ly); // function to create the 2D geometric mesh
//    
//    int order2DTri = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_TwoDTri = CreateComputationalTwoDTriMesh(gmesh_TwoDTri,order2DTri);
//    
//    ComputationShapeTwoDTri(cmesh_TwoDTri, order2DTri);
//    
//    // ---------------------------------------------------------------------------------------------------------
//    
//    long nnodesqu = 4; // number of divition
//    
//    TPZGeoMesh *gmesh_TwoDQuad = CreateTwoDQuadGMesh(nnodesqu, Lx, Ly); // function to create the 2D geometric mesh
//    
//    int order2DQuad = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_TwoDQuad = CreateComputationalTwoDQuadMesh(gmesh_TwoDQuad,order2DQuad);
//    
//    ComputationShapeTwoDQuad(cmesh_TwoDQuad, order2DQuad);
//    
//    // ********************************** (Create meshes: 3DTetra) *************************************************
//
//    long nnodestetra = 4; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDTetra = CreateThreeDTetraGMesh(nnodestetra, Lx, Ly, Lz); // create the 3D geometric mesh
//    
//    int order3DTetra = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_ThreeDTetra = CreateComputationalThreeDTetraMesh(gmesh_ThreeDTetra,order3DTetra);
//    
//    ComputationShapeThreeDTetra(cmesh_ThreeDTetra, order3DTetra);
//    
//    // -----------------------------------(Create meshes: 3DPyram) ---------------------------------------------
//
//    long nnodespyram = 5; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDPyram = CreateThreeDPyramGMesh(nnodespyram, Lx, Ly, Lz); // create the 3D geometric mesh
//    
//    int order3DPyram = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_ThreeDPyram = CreateComputationalThreeDPyramMesh(gmesh_ThreeDPyram,order3DPyram);
//    
//    ComputationShapeThreeDPyram(cmesh_ThreeDPyram, order3DPyram);
//    
//    // -----------------------------------(Create meshes: 3DPrism) ---------------------------------------------
//
//    long nnodesprism = 6; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDPrism = CreateThreeDPrismGMesh(nnodesprism, Lx, Ly, Lz); // create the 3D geometric mesh
//    
//    int order3DPrism = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_ThreeDPrism = CreateComputationalThreeDPrismMesh(gmesh_ThreeDPrism,order3DPrism);
//    
//    ComputationShapeThreeDPrism(cmesh_ThreeDPrism, order3DPrism);
//    
//    // -----------------------------------(Create meshes: 3DHexa) ---------------------------------------------
//
//    long nnodeshexa = 8; // number of divition
//    
//    TPZGeoMesh *gmesh_ThreeDHexa = CreateThreeDHexaGMesh(nnodeshexa, Lx, Ly, Lz); // create the 3D geometric mesh
//    
//    int order3DHexa = 1; // order of the piecewise polynomial space
//    TPZCompMesh * cmesh_ThreeDHexa = CreateComputationalThreeDHexaMesh(gmesh_ThreeDHexa,order3DHexa);
//    
//    ComputationShapeThreeDHexa(cmesh_ThreeDHexa, order3DHexa);
//
    
    
    // ************************************** (Create meshes: 1D) *****************************************************

    
    int orderOneD = 1; // order of the piecewise polynomial space
    TPZCompMesh * ComputationMesh_OneD = CreateComputationOneDMesh(gmesh_OneD,orderOneD);
    
    
    
return 0;
    
}



// ************************************** Create 1D meshes ***************************************

TPZGeoMesh *CreateOneDGMesh(long num_el, REAL size_el)
{
    TPZGeoMesh * gmesh_OneD = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 1; // geometry dimension
    std::string name("Geomesh OneD"); // geometry name
    
    gmesh_OneD->SetName(name);
    gmesh_OneD->SetDimension(geometry_dim);
    
    long num_nodes = num_el + 1; // Number of the nodes
    gmesh_OneD->NodeVec().Resize(num_nodes); // Resize of the geometry mesh
    
    
    const int physical_id = 1; // Define id for material
    const int bc0 = -1; // Define id for left boundary condition
    const int bc1 = -2; // Define id for right boundary condition
    
    for (long i = 0 ; i < num_nodes; i++)
    {
        const REAL valElem = i * size_el;
        TPZVec <REAL> coord(3,0.);
        coord[0] = valElem;
        gmesh_OneD->NodeVec()[i].SetCoord(coord); // Set of cordinate on the vector
        gmesh_OneD->NodeVec()[i].SetNodeId(i); // The id identification
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
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear> (element_id, Linear_topology, physical_id, *gmesh_OneD);
        
    }
    
    element_id++;
    
    
    // Left boundary condition
    point_topology[0] = 0;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc0, *gmesh_OneD);
    element_id++;
    
    
    // Right boundary condition
    point_topology[0] = num_nodes-1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (element_id, point_topology, bc1, *gmesh_OneD);
    
    
    gmesh_OneD->BuildConnectivity(); // Construct mesh neighbor connectivity
    
    
    std::ofstream outgmeshOneD("Geomesh_OneD.txt");
    gmesh_OneD->Print(outgmeshOneD);
    
    std::ofstream vtkgmeshOneD("Geomesh_OneD.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_OneD, vtkgmeshOneD);
    
    
    return gmesh_OneD;
}


// ******************* (Create cmesh: 1D) ***********************************************************


TPZCompMesh * CreateComputationalOneDMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_OneD = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_OneD->SetReference(geometry);
    cmesh_OneD->InsertMaterialObject(material);
    cmesh_OneD->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_OneD->SetAllCreateFunctionsContinuous();
    cmesh_OneD->AutoBuild();
    
    return cmesh_OneD;
}

// ******************* (Create shape: 1D) ***********************************************************


void ComputationShapeOneD(TPZCompMesh * cmesh, int order){
    
    
    
    //    int ncel = cmesh->NElements();
    
    TPZCompEl * cel = cmesh->Element(0);
    cel->SetOrthogonalFunction(pzshape::TPZShapeLinear::Chebyshev);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}


// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// ********************************** (Create meshes: 2D) *****************************************************

TPZGeoMesh *CreateTwoDTriGMesh(long nnodes, REAL Lx, REAL Ly)

{
    TPZGeoMesh * gmesh_TwoDTri = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    
    std::string name("geomesh TwoDTri"); // geometry name
    gmesh_TwoDTri->SetName(name);
    gmesh_TwoDTri->SetDimension(geometry_dim);
    
    
    gmesh_TwoDTri->NodeVec().Resize(nnodes); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodes);
    
    
    TPZVec<long> Triangle_topology(3);
    TPZVec <long> Linear_topology(2);
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_left = -3; // define id for a material (border left)
    
        
        // triangle element
        {
            
            coord[0] = 0.0; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[0].SetNodeId(0);
            gmesh_TwoDTri->NodeVec()[0].SetCoord(coord);
            Triangle_topology[0] = 0; // index
            
            coord[0] = Lx; // x coordinate
            coord[1] = 0.0; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[1].SetNodeId(1);
            gmesh_TwoDTri->NodeVec()[1].SetCoord(coord);
            Triangle_topology[1] = 1;
            
            coord[0] = 0.0; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDTri->NodeVec()[2].SetNodeId(2);
            gmesh_TwoDTri->NodeVec()[2].SetCoord(coord);
            Triangle_topology[2] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (elementid,Triangle_topology,physical_id,*gmesh_TwoDTri); // create triangle element
            elementid++;
            
        }


        // bottom
            {
                
                Linear_topology[0] = 0;
                
                Linear_topology[1] = 1;
                
                new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDTri); // create boundary element; bottom
                elementid++;
            }
        
        // right
        {
            
            Linear_topology[0] = 1;
            
            Linear_topology[1] = 2;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDTri); // create boundary element; bottom
            elementid++;
        }
        
        // left
        {
            
            Linear_topology[0] = 2;
            
            Linear_topology[1] = 0;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDTri); // create boundary element; bottom
            elementid++;
        }
        
    // Build the mesh
    gmesh_TwoDTri->BuildConnectivity();
    
    std::ofstream outgmeshTwoDTri("geomesh_TwoDTri.txt");
    gmesh_TwoDTri->Print(outgmeshTwoDTri);
    
    std::ofstream vtkgmeshTwoDTri("geomesh_TwoDTri.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDTri, vtkgmeshTwoDTri);
    return gmesh_TwoDTri;
    
}


// ******************* (Create cmesh: 2D) ***********************************************************


TPZCompMesh * CreateComputationalTwoDTriMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_TwoDTri = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_TwoDTri->SetReference(geometry);
    cmesh_TwoDTri->InsertMaterialObject(material);
    cmesh_TwoDTri->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_TwoDTri->SetAllCreateFunctionsContinuous();
    cmesh_TwoDTri->AutoBuild();
    
    return cmesh_TwoDTri;
}


// ******************* (Create shape: 2D) ***********************************************************


void ComputationShapeTwoDTri(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;

    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}

// ---------------------------------------------------------------------------------------------------------


TPZGeoMesh *CreateTwoDQuadGMesh(long nnodesqu, REAL Lx, REAL Ly)
{
    TPZGeoMesh * gmesh_TwoDQuad = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 2; // geometry dimension
    
    std::string name("geomesh TwoDQuad"); // geometry name
    gmesh_TwoDQuad->SetName(name);
    gmesh_TwoDQuad->SetDimension(geometry_dim);
    
    
    gmesh_TwoDQuad->NodeVec().Resize(nnodesqu); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodesqu);
    
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec <long> Linear_topology(2);
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_top = -3; // define id for a material (border top)
    const int bc_left = -4; // define id for a material (border left)
    
    
        // 0th element
        {
            
            coord[0] = -Lx; // x coordinate
            coord[1] = -Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[0].SetNodeId(0);
            gmesh_TwoDQuad->NodeVec()[0].SetCoord(coord);
            Quadrilateral_topology[0] = 0; // index
            
            coord[0] = Lx; // x coordinate
            coord[1] = -Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[1].SetNodeId(1);
            gmesh_TwoDQuad->NodeVec()[1].SetCoord(coord);
            Quadrilateral_topology[1] = 1;
            
            coord[0] = Lx; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[2].SetNodeId(2);
            gmesh_TwoDQuad->NodeVec()[2].SetCoord(coord);
            Quadrilateral_topology[2] = 2;
            
            coord[0] = -Lx; // x coordinate
            coord[1] = Ly; // Y coordinate
            coord[2] = 0.0; // Z coordinate
            gmesh_TwoDQuad->NodeVec()[3].SetNodeId(3);
            gmesh_TwoDQuad->NodeVec()[3].SetCoord(coord);
            Quadrilateral_topology[3] = 3;
            
            new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (elementid,Quadrilateral_topology,physical_id,*gmesh_TwoDQuad); // create quadrilateral element
            elementid++;
            
        }

    // bottom
    {
        
        Linear_topology[0] = 0;
        
        Linear_topology[1] = 1;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_bottom,*gmesh_TwoDQuad); // create boundary element; bottom
        elementid++;
    }
    // right
    {
        
        Linear_topology[0] = 1;
        
        Linear_topology[1] = 2;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_right,*gmesh_TwoDQuad); // create boundary element; bottom
        elementid++;
    }
    // top
    {
        
        Linear_topology[0] = 2;
        
        Linear_topology[1] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_top,*gmesh_TwoDQuad); // create boundary element; bottom
        elementid++;
    }
    // left
    {
        
        Linear_topology[0] = 3;
        
        Linear_topology[1] = 0;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementid,Linear_topology,bc_left,*gmesh_TwoDQuad); // create boundary element; bottom
        elementid++;
    }

    // Build the mesh
    gmesh_TwoDQuad->BuildConnectivity();
    
    std::ofstream outgmeshTwoDQuad("geomesh_TwoDQuad.txt");
    gmesh_TwoDQuad->Print(outgmeshTwoDQuad);
    
    std::ofstream vtkgmeshTwoDQuad("geomesh_TwoDQuad.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_TwoDQuad, vtkgmeshTwoDQuad);
    
    
    return gmesh_TwoDQuad;
    
}


// ******************* (Create cmesh: 2D) ***********************************************************


TPZCompMesh * CreateComputationalTwoDQuadMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_TwoDQuad = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_TwoDQuad->SetReference(geometry);
    cmesh_TwoDQuad->InsertMaterialObject(material);
    cmesh_TwoDQuad->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_TwoDQuad->SetAllCreateFunctionsContinuous();
    cmesh_TwoDQuad->AutoBuild();
    
    return cmesh_TwoDQuad;
}


// ******************* (Create shape: 2D) ***********************************************************


void ComputationShapeTwoDQuad(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;
    
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}

// ********************************** (Create meshes: 3D Tetra) *****************************************************

TPZGeoMesh *CreateThreeDTetraGMesh(long nnodestetra, REAL Lx, REAL Ly, REAL Lz)

{
    TPZGeoMesh * gmesh_ThreeDTetra = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    
    std::string name("geomesh ThreeDTetra"); // geometry name
    gmesh_ThreeDTetra->SetName(name);
    gmesh_ThreeDTetra->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDTetra->NodeVec().Resize(nnodestetra); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodestetra);
    
    
    TPZVec<long> Tetrahedron_topology(4);
    TPZVec<long> Triangle_topology(3);

    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_back = -3; // define id for a material (border back)
    const int bc_left = -4; // define id for a material (border left)
    
    // 0th element
    
    {
        
        // 1st node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDTetra->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDTetra->NodeVec()[0].SetCoord(coord);
        Tetrahedron_topology[0] = 0;
        
        // 2nd node
        
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDTetra->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDTetra->NodeVec()[1].SetCoord(coord);
        Tetrahedron_topology[1] = 1;
        
        // 3rd node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDTetra->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDTetra->NodeVec()[2].SetCoord(coord);
        Tetrahedron_topology[2] = 2;
        
        // 4th node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDTetra->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDTetra->NodeVec()[3].SetCoord(coord);
        Tetrahedron_topology[3] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (elementid, Tetrahedron_topology, physical_id, *gmesh_ThreeDTetra);
        elementid++;

    }
    
    // bottom
    {
        
    Triangle_topology[0] = 0;
    Triangle_topology[1] = 1;
    Triangle_topology[2] = 2;

    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_bottom, *gmesh_ThreeDTetra);
         elementid++;
    }
    
    // right
    {
        
        Triangle_topology[0] = 1;
        Triangle_topology[1] = 2;
        Triangle_topology[2] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_right, *gmesh_ThreeDTetra);
        elementid++;
    }
    
    // back
    {
        
        Triangle_topology[0] = 0;
        Triangle_topology[1] = 2;
        Triangle_topology[2] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_back, *gmesh_ThreeDTetra);
        elementid++;
    }
    
    // left
    {
        
        Triangle_topology[0] = 1;
        Triangle_topology[1] = 0;
        Triangle_topology[2] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_left, *gmesh_ThreeDTetra);
        elementid++;
    }
    
    
    // Build the mesh
    gmesh_ThreeDTetra->BuildConnectivity();
    
    std::ofstream outgmeshThreeDTetra("geomesh_ThreeDTetra.txt");
    gmesh_ThreeDTetra->Print(outgmeshThreeDTetra);
    
    std::ofstream vtkgmeshThreeDTetra("geomesh_ThreeDTetra.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDTetra, vtkgmeshThreeDTetra);
    
    
    return gmesh_ThreeDTetra;
    
}


// ---------------------------------------------------------------------------------------------------------


TPZCompMesh * CreateComputationalThreeDTetraMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_ThreeDTetra = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_ThreeDTetra->SetReference(geometry);
    cmesh_ThreeDTetra->InsertMaterialObject(material);
    cmesh_ThreeDTetra->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_ThreeDTetra->SetAllCreateFunctionsContinuous();
    cmesh_ThreeDTetra->AutoBuild();
    
    return cmesh_ThreeDTetra;
}


// ---------------------------------------------------------------------------------------------------------


void ComputationShapeThreeDTetra(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;
    par_space[2] = 0.25;

    
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}


// ********************************** (Create meshes: 3D Pyram) *****************************************************

TPZGeoMesh *CreateThreeDPyramGMesh(long nnodespyram, REAL Lx, REAL Ly, REAL Lz)

{
    TPZGeoMesh * gmesh_ThreeDPyram = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    
    std::string name("geomesh ThreeDPyram"); // geometry name
    gmesh_ThreeDPyram->SetName(name);
    gmesh_ThreeDPyram->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDPyram->NodeVec().Resize(nnodespyram); // Resize of the geometry mesh
    TPZVec<TPZGeoNode> Node(nnodespyram);
    
    
    TPZVec<long> Pyramid_topology(5);
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec<long> Triangle_topology(3);
    
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_back = -3; // define id for a material (border back)
    const int bc_left = -4; // define id for a material (border left)
    const int bc_front = -5; // define id for a material (border front)

    
    // 0th element
    
    {
        
        // 1st node
        
        coord[0] = Lx; // x coordinate
        coord[1] = -Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPyram->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDPyram->NodeVec()[0].SetCoord(coord);
        Pyramid_topology[0] = 0;
        
        // 2nd node
        
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPyram->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDPyram->NodeVec()[1].SetCoord(coord);
        Pyramid_topology[1] = 1;
        
        // 3rd node
        
        coord[0] = -Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPyram->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDPyram->NodeVec()[2].SetCoord(coord);
        Pyramid_topology[2] = 2;
        
        // 4th node
        
        coord[0] = -Lx; // x coordinate
        coord[1] = -Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPyram->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDPyram->NodeVec()[3].SetCoord(coord);
        Pyramid_topology[3] = 3;
        
        // 5th node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDPyram->NodeVec()[4].SetNodeId(4);
        gmesh_ThreeDPyram->NodeVec()[4].SetCoord(coord);
        Pyramid_topology[4] = 4;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoPyramid> (elementid, Pyramid_topology, physical_id, *gmesh_ThreeDPyram);
        elementid++;
        
    }
    
    // bottom
    {
        
        Quadrilateral_topology[0] = 0;
        Quadrilateral_topology[1] = 1;
        Quadrilateral_topology[2] = 2;
        Quadrilateral_topology[3] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_bottom, *gmesh_ThreeDPyram);
        elementid++;
    }
    
    // right
    {
        
        Triangle_topology[0] = 0;
        Triangle_topology[1] = 1;
        Triangle_topology[2] = 4;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_right, *gmesh_ThreeDPyram);
        elementid++;
    }
    
    // back
    {
        
        Triangle_topology[0] = 1;
        Triangle_topology[1] = 2;
        Triangle_topology[2] = 4;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_back, *gmesh_ThreeDPyram);
        elementid++;
    }
    
    // left
    {
        
        Triangle_topology[0] = 2;
        Triangle_topology[1] = 3;
        Triangle_topology[2] = 4;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_left, *gmesh_ThreeDPyram);
        elementid++;
    }
    
    // front
    {
        
        Triangle_topology[0] = 0;
        Triangle_topology[1] = 4;
        Triangle_topology[2] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_front, *gmesh_ThreeDPyram);
        elementid++;
    }
    
    // Build the mesh
    gmesh_ThreeDPyram->BuildConnectivity();
    
    std::ofstream outgmeshThreeDPyram("geomesh_ThreeDPyram.txt");
    gmesh_ThreeDPyram->Print(outgmeshThreeDPyram);
    
    std::ofstream vtkgmeshThreeDPyram("geomesh_ThreeDPyram.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDPyram, vtkgmeshThreeDPyram);
    
    
    return gmesh_ThreeDPyram;
    
}


// ---------------------------------------------------------------------------------------------------------


TPZCompMesh * CreateComputationalThreeDPyramMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_ThreeDPyram = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_ThreeDPyram->SetReference(geometry);
    cmesh_ThreeDPyram->InsertMaterialObject(material);
    cmesh_ThreeDPyram->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_ThreeDPyram->SetAllCreateFunctionsContinuous();
    cmesh_ThreeDPyram->AutoBuild();
    
    return cmesh_ThreeDPyram;
}


// ---------------------------------------------------------------------------------------------------------


void ComputationShapeThreeDPyram(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;
    par_space[2] = 0.25;
    
    
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}


// ********************************** (Create meshes: 3D Prism) *****************************************************

TPZGeoMesh *CreateThreeDPrismGMesh(long nnodesprism, REAL Lx, REAL Ly, REAL Lz)

{
    TPZGeoMesh * gmesh_ThreeDPrism = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    
    std::string name("geomesh ThreeDPrism"); // geometry name
    gmesh_ThreeDPrism->SetName(name);
    gmesh_ThreeDPrism->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDPrism->NodeVec().Resize(nnodesprism); // Resize of the geometry
    TPZVec<TPZGeoNode> Node(nnodesprism);
    
    
    TPZVec<long> Prism_topology(6);
    TPZVec <long> Quadrilateral_topology(4);
    TPZVec<long> Triangle_topology(3);
    
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_back = -3; // define id for a material (border back)
    const int bc_left = -4; // define id for a material (border left)
    const int bc_top = -5; // define id for a material (border top)
    
    
    // 0th element
    
    {
        
        // 1st node
        
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDPrism->NodeVec()[0].SetCoord(coord);
        Prism_topology[0] = 0;
        
        // 2nd node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDPrism->NodeVec()[1].SetCoord(coord);
        Prism_topology[1] = 1;
        
        // 3rd node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDPrism->NodeVec()[2].SetCoord(coord);
        Prism_topology[2] = 2;
        
        // 4th node
        
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDPrism->NodeVec()[3].SetCoord(coord);
        Prism_topology[3] = 3;
        
        // 5th node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[4].SetNodeId(4);
        gmesh_ThreeDPrism->NodeVec()[4].SetCoord(coord);
        Prism_topology[4] = 4;
        
        // 6th node
        
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDPrism->NodeVec()[5].SetNodeId(5);
        gmesh_ThreeDPrism->NodeVec()[5].SetCoord(coord);
        Prism_topology[5] = 5;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoPrism> (elementid, Prism_topology, physical_id, *gmesh_ThreeDPrism);
        elementid++;
        
    }
    
    // bottom
    {
        
        Triangle_topology[0] = 0;
        Triangle_topology[1] = 1;
        Triangle_topology[2] = 2;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_bottom, *gmesh_ThreeDPrism);
        elementid++;
    }

    // right
    {
        
        Quadrilateral_topology[0] = 0;
        Quadrilateral_topology[1] = 1;
        Quadrilateral_topology[2] = 4;
        Quadrilateral_topology[3] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_right, *gmesh_ThreeDPrism);
        elementid++;
    }
    
    
    // back
    {
        
        Quadrilateral_topology[0] = 1;
        Quadrilateral_topology[1] = 4;
        Quadrilateral_topology[2] = 5;
        Quadrilateral_topology[3] = 2;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_back, *gmesh_ThreeDPrism);
        elementid++;
    }
    
    
    // left
    {
        
        Quadrilateral_topology[0] = 0;
        Quadrilateral_topology[1] = 3;
        Quadrilateral_topology[2] = 5;
        Quadrilateral_topology[3] = 2;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_left, *gmesh_ThreeDPrism);
        elementid++;
    }
    
    
    // top
    {
        
        Triangle_topology[0] = 3;
        Triangle_topology[1] = 4;
        Triangle_topology[2] = 5;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementid, Triangle_topology, bc_top, *gmesh_ThreeDPrism);
        elementid++;
    }
    
    // Build the mesh
    gmesh_ThreeDPrism->BuildConnectivity();
    
    std::ofstream outgmeshThreeDPrism("geomesh_ThreeDPrism.txt");
    gmesh_ThreeDPrism->Print(outgmeshThreeDPrism);
    
    std::ofstream vtkgmeshThreeDPrism("geomesh_ThreeDPrism.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDPrism, vtkgmeshThreeDPrism);
    
    
    return gmesh_ThreeDPrism;
    
}


// ---------------------------------------------------------------------------------------------------------


TPZCompMesh * CreateComputationalThreeDPrismMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * cmesh_ThreeDPrism = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    cmesh_ThreeDPrism->SetReference(geometry);
    cmesh_ThreeDPrism->InsertMaterialObject(material);
    cmesh_ThreeDPrism->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    cmesh_ThreeDPrism->SetAllCreateFunctionsContinuous();
    cmesh_ThreeDPrism->AutoBuild();
    
    return cmesh_ThreeDPrism;
}


// ---------------------------------------------------------------------------------------------------------


void ComputationShapeThreeDPrism(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;
    par_space[2] = 0.25;
    
    
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}


// ********************************** (Create meshes: 3D Hexa) *****************************************************

TPZGeoMesh *CreateThreeDHexaGMesh(long nnodeshexa, REAL Lx, REAL Ly, REAL Lz)

{
    TPZGeoMesh * gmesh_ThreeDHexa = new TPZGeoMesh; // Initilized of TPZGeoMesh class
    
    long geometry_dim = 3; // geometry dimension
    
    std::string name("geomesh ThreeDHexa"); // geometry name
    gmesh_ThreeDHexa->SetName(name);
    gmesh_ThreeDHexa->SetDimension(geometry_dim);
    
    
    gmesh_ThreeDHexa->NodeVec().Resize(nnodeshexa); // Resize of the geometry
    TPZVec<TPZGeoNode> Node(nnodeshexa);
    
    
    TPZVec<long> Hexahedron_topology(8);
    TPZVec <long> Quadrilateral_topology(4);
    
    
    TPZVec<REAL> coord(3,0.0);
    
    // Index of element
    
    long elementid = 0;
    int physical_id = 1;
    
    // Index of boundary element
    const int bc_bottom = -1; // define id for a material (border bottom)
    const int bc_right = -2; // define id for a material (border right)
    const int bc_back = -3; // define id for a material (border back)
    const int bc_left = -4; // define id for a material (border left)
    const int bc_top = -5; // define id for a material (border top)
    const int bc_front = -6; // define id for a material (border front)

    
    
    // 0th element
    
    {
        
        // 1st node
        
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[0].SetNodeId(0);
        gmesh_ThreeDHexa->NodeVec()[0].SetCoord(coord);
        Hexahedron_topology[0] = 0;
        
        // 2nd node
        
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[1].SetNodeId(1);
        gmesh_ThreeDHexa->NodeVec()[1].SetCoord(coord);
        Hexahedron_topology[1] = 1;
        
        // 3rd node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[2].SetNodeId(2);
        gmesh_ThreeDHexa->NodeVec()[2].SetCoord(coord);
        Hexahedron_topology[2] = 2;
        
        // 4th node
        
        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = 0.0; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[3].SetNodeId(3);
        gmesh_ThreeDHexa->NodeVec()[3].SetCoord(coord);
        Hexahedron_topology[3] = 3;
        
        // 5th node
        
        coord[0] = Lx; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[4].SetNodeId(4);
        gmesh_ThreeDHexa->NodeVec()[4].SetCoord(coord);
        Hexahedron_topology[4] = 4;
        
        // 6th node
        
        coord[0] = Lx; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[5].SetNodeId(5);
        gmesh_ThreeDHexa->NodeVec()[5].SetCoord(coord);
        Hexahedron_topology[5] = 5;
        
        // 7th node

        coord[0] = 0.0; // x coordinate
        coord[1] = Ly; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[6].SetNodeId(6);
        gmesh_ThreeDHexa->NodeVec()[6].SetCoord(coord);
        Hexahedron_topology[6] = 6;
        
        // 8th node

        coord[0] = 0.0; // x coordinate
        coord[1] = 0.0; // Y coordinate
        coord[2] = Lz; // Z coordinate
        gmesh_ThreeDHexa->NodeVec()[7].SetNodeId(7);
        gmesh_ThreeDHexa->NodeVec()[7].SetCoord(coord);
        Hexahedron_topology[7] = 7;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoCube> (elementid, Hexahedron_topology, physical_id, *gmesh_ThreeDHexa);
        elementid++;
        
    }
    
    // bottom
    {
        
        Quadrilateral_topology[0] = 0;
        Quadrilateral_topology[1] = 1;
        Quadrilateral_topology[2] = 2;
        Quadrilateral_topology[3] = 3;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_bottom, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    // right
    {
        
        Quadrilateral_topology[0] = 1;
        Quadrilateral_topology[1] = 2;
        Quadrilateral_topology[2] = 6;
        Quadrilateral_topology[3] = 5;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_right, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    
    // back
    {
        
        Quadrilateral_topology[0] = 2;
        Quadrilateral_topology[1] = 3;
        Quadrilateral_topology[2] = 7;
        Quadrilateral_topology[3] = 6;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_back, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    
    // left
    {
        
        Quadrilateral_topology[0] = 3;
        Quadrilateral_topology[1] = 0;
        Quadrilateral_topology[2] = 4;
        Quadrilateral_topology[3] = 7;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_left, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    // front
    {
        
        
        Quadrilateral_topology[0] = 0;
        Quadrilateral_topology[1] = 1;
        Quadrilateral_topology[2] = 5;
        Quadrilateral_topology[3] = 4;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_front, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    // top
    {
        
        Quadrilateral_topology[0] = 4;
        Quadrilateral_topology[1] = 5;
        Quadrilateral_topology[2] = 6;
        Quadrilateral_topology[3] = 7;
        
        new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementid, Quadrilateral_topology, bc_top, *gmesh_ThreeDHexa);
        elementid++;
    }
    
    // Build the mesh
    gmesh_ThreeDHexa->BuildConnectivity();
    
    std::ofstream outgmeshThreeDHexa("geomesh_ThreeDHexa.txt");
    gmesh_ThreeDHexa->Print(outgmeshThreeDHexa);
    
    std::ofstream vtkgmeshThreeDHexa("geomesh_ThreeDHexa.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh_ThreeDHexa, vtkgmeshThreeDHexa);
    
    
    return gmesh_ThreeDHexa;
    
}


// ---------------------------------------------------------------------------------------------------------


TPZCompMesh * CreateComputationalThreeDHexaMesh(TPZGeoMesh * geometry, int order){
    
    TPZCompMesh * gmesh_ThreeDHexa = new TPZCompMesh;
    
    int material_id = 1;
    int dim = geometry->Dimension();
    int nstate = 1;
    TPZVec<STATE> sol;
    TPZL2Projection * material = new TPZL2Projection(material_id, dim, nstate, sol);
    
    gmesh_ThreeDHexa->SetReference(geometry);
    gmesh_ThreeDHexa->InsertMaterialObject(material);
    gmesh_ThreeDHexa->SetDefaultOrder(order);
    // States the functions related to the paper http://www.sciencedirect.com/science/article/pii/S0045782509000255
    gmesh_ThreeDHexa->SetAllCreateFunctionsContinuous();
    gmesh_ThreeDHexa->AutoBuild();
    
    return gmesh_ThreeDHexa;
}


// ---------------------------------------------------------------------------------------------------------


void ComputationShapeThreeDHexa(TPZCompMesh * cmesh, int order){
    
    
    TPZCompEl * cel = cmesh->Element(0);
    TPZInterpolationSpace * inteporlated_el = dynamic_cast< TPZInterpolationSpace * >(cel);
    
    int n_shapes    = inteporlated_el->NShapeF();
    int dim         = inteporlated_el->Dimension();
    
    TPZVec<REAL> par_space(dim);
    TPZFMatrix<REAL> phi(n_shapes,1);
    TPZFMatrix<REAL> dphidxi(dim,n_shapes);
    par_space[0] = 0.25;
    par_space[1] = 0.25;
    par_space[2] = 0.25;
    
    
    
    inteporlated_el->Shape(par_space, phi, dphidxi);
    phi.Print(std::cout);
    dphidxi.Print(std::cout);
}
// ---------------------------------------------------------------------------------------------------------



TPZCompMesh * CreateComputationOneDMesh(TPZGeoMesh * gmesh_OneD,int orderOneD){
    
    TPZCompMesh * ComputationMesh_OneD = new TPZCompMesh;
    
    std::string name("Commesh OneDMesh"); // Computation mesh name
    ComputationMesh_OneD->SetName(name);
    

    return ComputationMesh_OneD;

}

