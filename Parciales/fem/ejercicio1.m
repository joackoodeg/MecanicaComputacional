 % Inciso 3: dos elementos lineales, un paso de Backward Euler con dt = 0.1
%   rho*cp*T_t = k*T'' + G,  k = 1, G = 100, rho*cp = 1
%   T(0) = 0 (Dirichlet);  T'(1) = -T(1)  ->  Robin h = 1, Tinf = 0
%   T(x,0) = 10 x (1-x)^2
model.k = 1;  model.c = 0;  model.rho = 1;  model.cp = 1;
model.G = 100;
model.tipoFF = 1;                      % lineales

xnode = [0 0.5 1];                     % 2 elementos lineales
cb = [1 0 -1;                          % izquierda: Dirichlet T = 0
      3 1  0];                         % derecha:   Robin h = 1, Tinf = 0

T0 = 10*xnode.*(1 - xnode).^2;         % condicion inicial en los nodos: 0  1.25  0
dt = 0.1;
et = [2, 1, 1e-12, dt];                % [Backward Euler, 1 paso, tol, dt]

T = elemFinitos(xnode, model, cb, et, T0);
disp('T en t = 0 y t = 0.1 (columnas):'); disp(T);
printf('Primer paso: T = [%g  %.4f  %.4f]\n', T(:,2));
