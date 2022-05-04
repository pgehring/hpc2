% runExample4 
% example to illustrate the use of L2projection
clear all
%*** input data
f = @(x)    ones(size(x,1),1);
g = @(x)    ones(size(x,1),1);
uD = @(x)   zeros(size(x,1),1);
%*** build geometry
geom = [2  2  2  2  2  2 
        0  2  2  1  1  0 
        2  2  1  1  0  0 
        0  0  1  1  2  2 
        0  1  1  2  2  0 
        1  1  1  1  1  1    
        0  0  0  0  0  0];
[p,e,t]=initmesh(geom,'hmax',0.1);
elements = t(1:3,:)'; 
coordinates = p';
%*** select edges on boundary
boundary = { [3,4] , [1,2,5,6] };      
for k=1:length(boundary)
  idx = find(ismember(e(5,:),boundary{k}));
  bdry{k} = e(1:2,idx)';
  idx = find( e(6,idx) == 0 );
  bdry{k}(idx,[2,1]) = bdry{1}(idx,:); 
end
dirichlet = bdry{1};neumann = bdry{1};
%*** solve Poisson's equation
x = solvePoisson(coordinates,elements, ...
                    dirichlet,neumann,f,g,uD);
%*** plot solution
trisurf(elements,coordinates(:,1),coordinates(:,2),...
         x,'edgecolor','k','facecolor','interp')
title('Solution of Poisson''s equation')  
view(100,10)
