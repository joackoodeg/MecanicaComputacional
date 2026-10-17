function [Nfun,b,c,A] = triShapeCST(xe)
    % Funciones de forma lineales completas (CST), con b,c (gradientes)
    % y A (area con signo, requiere orden CCW de nodos)
    x = xe(:,1); y = xe(:,2);
    A = 0.5*det([1 x(1) y(1); 1 x(2) y(2); 1 x(3) y(3)]);
    a = [x(2)*y(3)-x(3)*y(2); x(3)*y(1)-x(1)*y(3); x(1)*y(2)-x(2)*y(1)];
    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    Nfun = @(px,py) (a + b*px + c*py)/(2*A);
end
