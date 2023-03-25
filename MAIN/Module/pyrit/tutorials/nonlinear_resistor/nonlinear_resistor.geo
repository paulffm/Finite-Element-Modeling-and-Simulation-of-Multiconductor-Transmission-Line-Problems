// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 1
//
//  Variables, elementary entities (points, curves, surfaces), physical groups
//  (points, curves, surfaces)
//
// -----------------------------------------------------------------------------

// The simplest construction in Gmsh's scripting language is the
// `affectation'. The following command defines a new variable `lc':

r_i = 0.01;
r_m = 0.02;
r_o = 0.03;
height = 0.01;
lc =  (r_o - r_i) / 20;

Point(1) = {r_i,0,0,lc};
Point(2) = {r_m,0,0,lc};
Point(3) = {r_o,0,0,lc};
Point(4) = {r_o,height,0,lc};
Point(5) = {r_m,height,0,lc};
Point(6) = {r_i,height,0,lc};

Line(1) = {1,2};
Line(2) = {2,5};
Line(3) = {5,6};
Line(4) = {6,1};
Line(5) = {2,3};
Line(6) = {3,4};
Line(7) = {4,5};

Curve Loop(1) = {1,2,3,4};
Curve Loop(2) = {5,6,7,-2};

Plane Surface(1) = {1};
Plane Surface(2) = {2};


Physical Surface("NONLINEAR_MATERIAL") = {1,2};

Physical Curve("GROUND") = {4};
Physical Curve("VOLTAGE") = {6};




