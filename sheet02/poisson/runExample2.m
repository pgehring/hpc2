% runExample2 run an example to illustrate 
% the use of initmesh
clear all
g = [2  2  2  2  2  2 
     0  2  2  1  1  0 
     2  2  1  1  0  0 
     0  0  1  1  2  2 
     0  1  1  2  2  0 
     1  1  1  1  1  1    
     0  0  0  0  0  0];
[p,e,t]=initmesh(g,'hmax',0.5);
elements = t(1:3,:)';
coordinates = p';
trisurf(elements,coordinates(:,1),coordinates(:,2), ...
    0*coordinates(:,2),'edgecolor','k','facecolor','none')
view(2), axis equal, axis off