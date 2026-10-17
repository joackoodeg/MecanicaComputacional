function [H] = MatH_simb(xnode,rho,cp)
    if (size(xnode,1) == 3) % elemento triangular
        A = 0.5*det([1 xnode(1,1) xnode(1,2);
                     1 xnode(2,1) xnode(2,2);
                     1 xnode(3,1) xnode(3,2);]);
        H = rho*cp*A/12*[2 1 1;1 2 1;1 1 2];
    else
        syms e n
        N1 = 0.25*(1-e)*(1-n); N2 = 0.25*(1+e)*(1-n);
        N3 = 0.25*(1+e)*(1+n); N4 = 0.25*(1-e)*(1+n);
        N = [N1 N2 N3 N4];
        % x = xnode(1,1)*N1+xnode(2,1)*N2+xnode(3,1)*N3+xnode(4,1)*N4;
        x = N*xnode(:,1);
        % y = xnode(1,2)*N1+xnode(2,2)*N2+xnode(3,2)*N3+xnode(4,2)*N4;
        y = N*xnode(:,2);
        J = [diff(x,e)  diff(y,e); diff(x,n)  diff(y,n)];
        detJ = det(J);
        integ = N'*N*detJ;
        H = rho*cp*double(int(int(integ,e,-1,1),n,-1,1));
    end
end