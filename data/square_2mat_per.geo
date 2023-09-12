//https://bthierry.pages.math.cnrs.fr/tutorial/gmsh/occ_basics/
//http://onelab.info/pipermail/gmsh/2017/011277.html
//https://fossies.org/linux/gmsh/demos/boolean/periodic.geo
//https://scicomp.stackexchange.com/questions/30651/gmsh-for-3d-volume-with-inclusions

level=0;
cl=.25;

Point(1) = { 0  ,  0, level, cl} ; 
Point(2) = { 0.5,  0, level, cl} ;
Point(3) = { 1  ,  0, level, cl} ;
Point(4) = { 1  ,  1, level, cl} ;
Point(5) = { 0.5,  1, level, cl} ;
Point(6) = { 0  ,  1, level, cl} ;

Line(1) = {1,2}; 
Line(2) = {2,3};
Line(3) = {3,4}; 
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {2,5};
Line Loop(1) = {  1, 7, 5, 6 } ;
Plane Surface(1) = {1} ;
Line Loop(2) = {  2, 3, 4, -7} ;
Plane Surface(2) = {2} ;

Periodic Line {5} = {-1} Translate {0,1,0};
Periodic Line {4} = {-2} Translate {0,1,0};
Periodic Line {3} = {-6} Translate {1,0,0};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};
Physical Line(5) = {5};
Physical Line(6) = {6};
Physical Line(7) = {7};

Physical Surface(1) = {1};
Physical Surface(2) = {2};

For pt In {1:22}
    Physical Point(pt)  = {pt}; 
EndFor

//Transfinite Surface {1,2};
//Recombine Surface {1,2};


Mesh 2;
Mesh.MshFileVersion = 2.2;
