

L = 2.2;
H = 0.41;
c_x = 0.2;
c_y = 0.2;
r = 0.05;
gdim = 2;
ps   = 0.025;
psc  = 0.025;

res_min = r/4;

Point(1) = {0,0,0,ps};
Point(2) = {L,0,0,ps};
Point(3) = {L,H,0,ps};
Point(4) = {0,H,0,ps};

Point(5) = {c_x, c_y, 0, psc};
Point(6) = {c_x, c_y+r,0, psc};
Point(7) = {c_x, c_y-r,0, psc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

//+
Circle(5) = {6, 5, 7};
//+
Circle(6) = {7, 5, 6};


//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Curve Loop(2) = {6, 5};
//+
Plane Surface(1) = {1, 2};
//+
Physical Curve("WALLS", 1) = {3, 1};
//+
Physical Curve("INFLOW", 2) = {4};
//+
Physical Curve("OUTFLOW", 3) = {2};
//+
Physical Curve("OBSTACLE", 4) = {6, 5};
//+
Physical Surface("FLUID", 5) = {1};



//+
Field[1] = Distance;
//+

Field[1].CurvesList = {5, 6};


//+
Field[2] = Threshold;
//+
Field[2].InField = 1;

//+
Field[2].DistMax = 2*H;
//+
Field[2].DistMin = r;
//+
Field[2].SizeMax = 0.25*H;
//+
Field[2].SizeMin = res_min;



//+
Field[3] = Min;
//+
Field[3].FieldsList = {2};
//+
Background Field = 3;
//+

