// cable
Mesh.MshFileVersion = 2.2;

// Parameters
r1 = 0.002;
r2 = 0.0035;

SetFactory("OpenCASCADE");

// Wire
Circle(1) = {0, 0, 0, r1, 0, 2*Pi};

Curve Loop(1) = {1};
Plane Surface(1) = {1};
Physical Surface("WIRE") = {1};

// Shell
Circle(2) = {0, 0, 0, r2, 0, 2*Pi};
Curve Loop(2) = {2};
Plane Surface(2) = {2,1};
Physical Surface("SHELL") = {2};

// Boundary condition
Physical Line("GND") = {2};
