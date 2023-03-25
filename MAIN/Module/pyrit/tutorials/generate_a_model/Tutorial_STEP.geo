SetFactory("OpenCASCADE");
v() = ShapeFromFile("Tutorial_STEP.stp");
BooleanFragments{ Surface{v()}; Delete; }{}

Physical Surface('rotor') = {1};
Physical Surface('steel') = {2};
Physical Surface('stator') = {3};
Physical Line('external_temperature') = {5,6};