g = 'lshapeg';
[p,e,t]=initmesh(g,'hmax',0.5);
elements = t(1:3,:)';
coordinates = p';
trimesh(elements,coordinates(:,1),coordinates(:,2), ...
                 0*coordinates(:,2),'edgecolor','k')
view(2), axis equal, axis off