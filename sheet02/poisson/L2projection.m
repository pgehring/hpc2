function x = L2projection(coordinates,elements,f)
%*** Assembly of mass matrix
M = sparse(size(coordinates,1),size(coordinates,1));
for k = 1:size(elements,1)
    nodes = elements(k,:);
    area = det([1 1 1 ; coordinates(nodes,:)'])/2;
    M_k = area/12 * [2,1,1;1,2,1;1,1,2];
    M(nodes,nodes) = M(nodes,nodes) + M_k;
end
%*** Assembly of right-hand side
b = zeros(size(coordinates,1),1);
for k = 1:size(elements,1)
    nodes = elements(k,:);
    area = det([1 1 1 ; coordinates(nodes,:)'])/2;
    f_k = f(coordinates(nodes,:));
    b(nodes) = b(nodes) + area/3 * f_k;
end
%*** Computation of L2-projection 
x = M\b;
