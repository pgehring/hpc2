function x = solvePoisson(coordinates,elements, ...
                                       dirichlet,neumann,f,g,uD)
%*** Assembly of stiffness matrix
A = sparse(size(coordinates,1),size(coordinates,1));
for i = 1:size(elements,1)
    nodes = elements(i,:);
    B = [1 1 1 ; coordinates(nodes,:)'];
    grad = B \ [0 0 ; 1 0 ; 0 1];
    A(nodes,nodes) = A(nodes,nodes) + det(B)*grad*grad'/2;
end
%*** Prescribe values at Dirichlet nodes
dirichlet = unique(dirichlet);
x = zeros(size(coordinates,1),1);
x(dirichlet) = uD(coordinates(dirichlet,:));
%*** Assembly of right-hand side
b = -A*x;
for i = 1:size(elements,1)
    nodes = elements(i,:);
    sT = [1 1 1]*coordinates(nodes,:)/3;
    b(nodes) = b(nodes) + det([1 1 1 ; coordinates(nodes,:)'])*f(sT)/6;
end
for i = 1:size(neumann,1)
    nodes = neumann(i,:);
    mE = [1 1]*coordinates(nodes,:)/2;
    b(nodes) = b(nodes) + norm([1 -1]*coordinates(nodes,:))*g(mE)/2;
end
%*** Computation of P1-FEM approximation
freenodes = setdiff(1:size(coordinates,1), dirichlet);
x(freenodes) = A(freenodes,freenodes)\b(freenodes);

