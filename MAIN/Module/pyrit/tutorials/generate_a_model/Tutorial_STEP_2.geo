Merge "Tutorial_STEP.stp";
//+
SetFactory("OpenCASCADE");
Delete {Surface{1};}
Delete {Line{1:4};}
Line Loop(15) = {5,6};
Line Loop(16) = {7,8};
Plane Surface(4) = {15,16};
Physical Surface('rotor') = {4};
Physical Surface('steel') = {2};
Physical Surface('stator') = {3};
Physical Line('external_temperature')={9,10};