function [K,F,Ke,Fe,T] = resolverMalla(xnode,conn,k,Q,cN,dirNodes,dirVal)
    Nnod = size(xnode,1);
    Nele = size(conn,1);
    K = zeros(Nnod,Nnod);
    F = zeros(Nnod,1);
    Ke = cell(Nele,1); Fe = cell(Nele,1);

    for e = 1:Nele
        nod = conn(e,:);
        xe  = xnode(nod,:);
        Ke{e} = MatCondCar(xe,k,k);
        Fe{e} = genF(xe,Q);
        K = assembly(K,Ke{e},nod);
        F = assembly(F,Fe{e},nod);
    end

    % Condicion de Neumann
    % neuCond.m calcula F(lado) -= q*l/2 (q = flujo saliente, convencion
    % de Fourier). La solucion oficial del parcial usa q positivo como
    % flujo ENTRANTE, por lo que se invierte el signo de q al llamar
    % a la funcion del curso (sin modificar neuCond.m).
    cN(:,3) = -cN(:,3);
    F = neuCond(F,cN,xnode);

    % Condicion de Dirichlet (reemplazo de fila)
    for i = 1:length(dirNodes)
        n = dirNodes(i);
        K(n,:) = 0; K(n,n) = 1;
        F(n) = dirVal(i);
    end

    T = K\F;
end
