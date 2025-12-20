// Gmsh project created on Sat Dec 13 22:44:03 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {2, 0, 0, 1.0};
//+
Point(2) = {10, 0, 0, 1.0};
//+
Point(3) = {10, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 0.5, 0, 1.0};
//+
Point(6) = {2, 0.5, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 1};
//+
Curve Loop(1) = {3, 4, 5, 6, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("BOTTOM", 7) = {1};
//+
Physical Curve("RIGHT", 8) = {2};
//+
Physical Curve("TOP", 9) = {3};
//+
Physical Curve("LEFT", 10) = {4};
//+
Physical Curve("STEP_HORIZONTAL", 11) = {5};
//+
Physical Curve("STEP_VERTICAL", 12) = {6};
//+
Physical Surface("FLUID", 13) = {1};
