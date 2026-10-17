function evaluarPunto(xnode,conn,T,punto,k)
    px = punto(1); py = punto(2);
    Nele = size(conn,1);
    for e = 1:Nele
        nod = conn(e,:);
        xe = xnode(nod,:);
        [Nfun,b,c,A] = triShapeFull(xe);
        Nv = Nfun(px,py);
        % Coordenadas baricentricas: dentro del triangulo si todas >= -tol
        if all(Nv >= -1e-9)
            Te = T(nod);
            Tval = Nv'*Te;
            qx = -k*(b'*Te)/(2*A);
            qy = -k*(c'*Te)/(2*A);
            fprintf('\nPunto (%.2f,%.2f) pertenece al elemento con nodos [%d %d %d]\n',px,py,nod(1),nod(2),nod(3));
            fprintf('  N = [%.4f %.4f %.4f]\n',Nv(1),Nv(2),Nv(3));
            fprintf('  T(%.2f,%.2f) = %.4f\n',px,py,Tval);
            fprintf('  qx = %.4f   qy = %.4f\n',qx,qy);
        end
    end
end
