lc = .5;
Point(1) = {0, 0, 0, lc};
Point(2) = {0, 1, 0, lc};
Point(3) = {8, 1, 0, lc};
Point(4) = {8, 0, 0, lc};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {2, 3, 4, 1};
Plane Surface(1) = {1};
Transfinite Curve {2, 4} = 8 Using Progression 1;
Transfinite Curve {1, 3} = 3 Using Progression 1;
Transfinite Surface {1} = {1,2,3,4};
Recombine Surface {1};

Physical Surface(1) = {1};
Physical Curve(1) = {1};
Physical Curve(2) = {3};
Physical Curve(3) = {2,4};

Mesh 2;
Mesh.MshFileVersion = 2.2;
Mesh.ElementOrder = 4;
