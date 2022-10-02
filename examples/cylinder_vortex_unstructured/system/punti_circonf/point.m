#!/usr/bin/octave -q

np = 360;
r=   75;

x=zeros(np,1);
y=zeros(np,1);
z=zeros(np,1);
theta=zeros(np,1);
%
   for i=1:np
     theta(i) = (2*pi/(np-1))*(i-1);
     x(i) = r*cos(pi - theta(i)) ;
     y(i) = r*sin(pi - theta(i)) ;
   end


dlmwrite('file_cerchio.dat', [x,y,z] ,  ' ');

