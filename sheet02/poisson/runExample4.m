% runExample4 
% example to illustrate the use of L2projection
clear all
f = @(x) x(:,1) .* x(:,2);
coordinates = [0,0;1,0;1,1;0,1;1/2,1/2];
elements = [1,2,5;2,3,5;3,4,5;4,1,5];
x = L2projection(coordinates,elements,f);
trisurf(elements,coordinates(:,1),coordinates(:,2),...
         x,'edgecolor','k','facecolor','interp')
title('L^2 projection of x*y')  
view(6,26)



