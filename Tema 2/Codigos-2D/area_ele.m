function [A] = area_ele(nodes)
    if (size(nodes,1) == 3) % Triangular element
        A = 0.5*det([1 nodes(1,1) nodes(1,2);
                     1 nodes(2,1) nodes(2,2);
                     1 nodes(3,1) nodes(3,2);]);
    else % Quadrangular element
        A1 = 0.5*det([1 nodes(1,1) nodes(1,2);
                     1 nodes(2,1) nodes(2,2);
                     1 nodes(3,1) nodes(3,2);]);
        A2 = 0.5*det([1 nodes(1,1) nodes(1,2);
                     1 nodes(3,1) nodes(3,2);
                     1 nodes(4,1) nodes(4,2);]);
        A = A1+A2; 
    end
end