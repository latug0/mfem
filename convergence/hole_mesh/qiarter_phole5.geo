lc=0.05;
N=110;
M=N+2;
Point(1)={0,0,0,lc};

ray=.25;
cx=.5;
cy=.5;

//morceaux d’ellipse
Point(5)={cx,cy,0,lc};
Point(6)={cx-ray,cy,0,lc};
Point(8)={cx-ray,cy,0,lc};
Point(9)={cx,cy-ray,0,lc};

Point(10)={0,0.5,0,lc};
Point(11)={0.5,0,0,lc};

Point(14)={cx-ray*(Sqrt(2)/2),cy-ray*(Sqrt(2)/2),0,lc};
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
Physical Curve(1)={34};	//gauche
Physical Curve(2)={33};	//haut
Physical Curve(3)={17,18};	//demicercle
Physical Curve(4)={30};	//droite
Physical Curve(5)={29};	//bas

//maillage quadrangle
Transfinite Curve {17,18,29,34} = M;
Transfinite Curve {30,32,33} =  N;
Transfinite Surface {1,2};
Recombine Surface {1,2};

