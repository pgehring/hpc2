% runExample3 run an example to illustrate 
% the use of initmesh
clear all
g = [2.0000   2.0000    1.0000    1.0000    1.0000
          0        0   -1.0000    0.0000    0.0000
     1.0000        0    0.0000    1.0000   -1.0000
          0   1.0000   -0.0000   -1.0000    1.0000
          0        0   -1.0000         0   -0.0000
          0        0    1.0000    1.0000    1.0000
     1.0000   1.0000         0         0         0
          0        0         0         0         0
          0        0         0         0         0
          0        0    1.0000    1.0000    1.0000];
[p,e,t]=initmesh(g,'hmax',0.3);
elements = t(1:3,:)';
coordinates = p';
%*** plot triangular mesh
trisurf(elements,coordinates(:,1),coordinates(:,2), ...
                 0*coordinates(:,2),'edgecolor','k','facecolor','none')

%*** select edges on boundary
%    indicated by there original geometric edges in g 
boundary = { [1,2] , [3,4,5] };      
for k=1:length(boundary)
  idx = find(ismember(e(5,:),boundary{k}));
  bdry{k} = e(1:2,idx)';
  idx = find( e(6,idx) == 0 );
  bdry{k}(idx,[2,1]) = bdry{1}(idx,:); 
end
%*** plot boundary of mesh
hold on
plot(reshape(coordinates(bdry{1},1),[],2)', ...
     reshape(coordinates(bdry{1},2),[],2)','r-','linewidth',3)
plot(reshape(coordinates(bdry{2},1),[],2)', ...
     reshape(coordinates(bdry{2},2),[],2)','g:','linewidth',3)
hold off
view(2), axis equal, axis off