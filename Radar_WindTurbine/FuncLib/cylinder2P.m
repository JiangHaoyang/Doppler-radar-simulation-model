function [X,Y,Z]=cylinder2P(P1,P2,h,r,N)

Cntr = (P1+P2)/2;  % ellipsoid center
Lc = norm(P2-P1);

% the axis defined by: P1+V*[0:Lc]
V = (P1-P2)/Lc;   %normalized cylinder's axis-vector;
U = rand(1,3);     %linear independent vector
U = V-U/(U*V');    %orthogonal vector to V
U = U/sqrt(U*U');  %orthonormal vector to V
W = cross(V,U);    %vector orthonormal to V and U
W = W/sqrt(W*W');  %orthonormal vector to V and U 

% generate the ellipsoid at (0,0,0)
[Xc,Yc,Zc] = cylinder(r,N);
Zc=(Zc-0.5)*h;

A = kron(U',Xc);
B = kron(W',Yc);
C = kron(V',Zc);
TMP = A+B+C;

X = TMP([1,2],:)+Cntr(1);
Y = TMP([3,4],:)+Cntr(2);
Z = TMP([5,6],:)+Cntr(3);

end