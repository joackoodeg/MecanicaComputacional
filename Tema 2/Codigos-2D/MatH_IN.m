function [H] = MatH_IN(xnode,rho,cp)
    if (size(xnode,1) == 3) % elemento triangular
        A = 0.5*det([1 xnode(1,1) xnode(1,2);
                     1 xnode(2,1) xnode(2,2);
                     1 xnode(3,1) xnode(3,2);]);
        H = rho*cp*A/12*[2 1 1;1 2 1;1 1 2];
    else
        N = @(s,t)[0.25*(1-s)*(1-t), 0.25*(1+s)*(1-t), 0.25*(1+s)*(1+t), 0.25*(1-s)*(1+t)];
        DN = @(s,t)[(-1+t)/4,( 1-t)/4,( 1+t)/4,(-1-t)/4;
                    (-1+s)/4,(-1-s)/4,( 1+s)/4,( 1-s)/4]; 
        % cuatro puntos de Gauss con peso w=1
        p = sqrt(3)/3;
        pospg = [-p,p];
        H = zeros(4,4);
        for i=1:2
            for j=1:2
                Nnum = N(pospg(i),pospg(j));
                DNnum = DN(pospg(i),pospg(j));
                J = DNnum*xnode;
                H = H + Nnum'*Nnum*det(J);
            end
        end
        H = rho*cp*H;
    end

end