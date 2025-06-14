// Gmsh project created on Sat Jun 14 11:34:30 2025
//+
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve("DOWN", 5) = {1};
Physical Curve("RIGHT", 6) = {2};
Physical Curve("TOP", 7) = {3};
Physical Curve("LEFT", 8) = {4};
//+
Physical Surface("FLUID", 9) = {1};
