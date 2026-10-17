function [K] = MatCondNat(xnode,kx,ky)
    k = [kx 0;0 ky];
    if (size(xnode,1) == 3) % elemento triangular
        J = [xnode(2,1)-xnode(1,1)  xnode(2,2)-xnode(1,2);
             xnode(3,1)-xnode(1,1)  xnode(3,2)-xnode(1,2)];
        DN = [-1 1 0;-1 0 1];
        V = inv(J)*DN;
        A = 0.5;
        K = (V'*k*V)*det(J)*A;
    else
        syms e n
        N1 = 0.25*(1-e)*(1-n); N2 = 0.25*(1+e)*(1-n);
        N3 = 0.25*(1+e)*(1+n); N4 = 0.25*(1-e)*(1+n);
        N = [N1 N2 N3 N4];
        x = xnode(1,1)*N1+xnode(2,1)*N2+xnode(3,1)*N3+xnode(4,1)*N4;
        y = xnode(1,2)*N1+xnode(2,2)*N2+xnode(3,2)*N3+xnode(4,2)*N4;
        J = [diff(x,e)  diff(y,e); diff(x,n)  diff(y,n)];
        DNde = diff(N,e); DNdn = diff(N,n); DN = [DNde; DNdn];
        V = inv(J)*DN;
        integ = V'*k*V*det(J);
        K = double(int(int(integ,e,-1,1),n,-1,1)); 
end