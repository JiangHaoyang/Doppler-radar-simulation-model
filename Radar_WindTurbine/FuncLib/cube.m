function [ X,Y,Z ] = cube( center,a,b,c )
%CUBE Summary of this function goes here
%   Detailed explanation goes here
x=[0,0,0,0,0;
   a/2,a/2,-a/2,-a/2,a/2;
   a/2,a/2,-a/2,-a/2,a/2;
   0,0,0,0,0];
y=[0,0,0,0,0;
   -b/2,b/2,b/2,-b/2,-b/2;
   -b/2,b/2,b/2,-b/2,-b/2;
   0,0,0,0,0];
z=[c/2,c/2,c/2,c/2,c/2;
   c/2,c/2,c/2,c/2,c/2;
   -c/2,-c/2,-c/2,-c/2,-c/2;
   -c/2,-c/2,-c/2,-c/2,-c/2];
X=x+center(1);
Y=y+center(2);
Z=z+center(3);

end

