function C = quadShapeBilineal(xe)
    % Funciones de forma bilineales "fisicas" N_i(x,y) = a_i+b_i*x+c_i*y+d_i*x*y
    % para un cuadrilatero de 4 nodos (orden CCW), resolviendo el sistema
    % lineal que impone N_i(nodo_j) = delta_ij.
    % Devuelve C (4x4): columna i = [a_i;b_i;c_i;d_i]
    x = xe(:,1); y = xe(:,2);
    Amat = [ones(4,1) x y x.*y];
    C = Amat\eye(4);
end
