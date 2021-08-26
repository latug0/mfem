lc=0.05;
N = 9;
M = 5;
N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

N = N*2-1;
M = M*2-1;

ray=1.;
cx=10.;
cy=10.;
Point(1)={-cx,-cy,0,lc};

//morceaux d’ellipse
Point(5)={0,0,0,lc};
Point(6)={-ray,0,0,lc};
Point(9)={-ray,0,0,lc};
Point(8)={0,-ray,0,lc};

Point(10)={0,-cy,0,lc};
Point(11)={-cx,0,0,lc};

Point(14)={-ray*(Sqrt(2)/2),-ray*(Sqrt(2)/2),0,lc};
Ellipse(17)={8,5,8,14};
Ellipse(18)={14,5,14,9};

Line(29)={1,11};
Line(30)={11,9};
Line(32)={14,1};
Line(33)={8,10};
Line(34)={10,1};
//generationdescontours
Line Loop(3)={29,30,-18,32};
Line Loop(4)={-32,-17,33,34};

//generationdelasurface
Plane Surface(1)={3};
Plane Surface(2)={4};
Physical Surface(1)={1,2};
Physical Curve(1)={34};	//bas
Physical Curve(2)={33};	//haut
Physical Curve(3)={17,18};	//demicercle
Physical Curve(4)={30};	//droite
Physical Curve(5)={29};	//gauche

//maillage quadrangle
Transfinite Curve {17,18,29,34} = M;
Transfinite Curve {30,32,33} =  N;
Transfinite Surface {1,2};
Recombine Surface {1,2};

