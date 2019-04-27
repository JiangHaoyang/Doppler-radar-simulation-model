function [X,Y,Z] = ellipsoid2P1A(P1,P2,rotate,a,b,c,N)

Cntr = (P1+P2)/2;  % ellipsoid center
% generate the ellipsoid at (0,0,0)
[Xc,Yc,Zc] = ellipsoid(0,0,0,a,b,c,N);

% roll for a certain angle
U=[1,0,0];
V=[0,cos(rotate),sin(rotate)];
W=[0,-sin(rotate),cos(rotate)];
A = kron(U',Xc);
B = kron(V',Yc);
C = kron(W',Zc);
TMP = A+B+C;
nt = size(TMP,2);
X = TMP(1:nt,:);
Y = TMP(nt+1:2*nt,:);
Z = TMP(2*nt+1:end,:);

% new base vector
Lc = norm(P2-P1);
U = (P1-P2)/Lc;   %normalized
W = [0,0,1];
V = cross(W,U);    %vector orthonormal to V and U
V = V/sqrt(V*V');  %normalized
W = cross(U,V);%vector orthonormal to W and U
W = W/sqrt(W*W');  %normalized

A = kron(U',X);
B = kron(V',Y);
C = kron(W',Z);
TMP = A+B+C;
nt = size(TMP,2);

X = TMP(1:nt,:)+Cntr(1);
Y = TMP(nt+1:2*nt,:)+Cntr(2);
Z = TMP(2*nt+1:end,:)+Cntr(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%