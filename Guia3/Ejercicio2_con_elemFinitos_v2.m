% Ejercicio 2 (Guia 3) resuelto con elemFinitos.m
%   k T'' + G = 0 en [-1, 1], k = 2, G = 10 x^3
%   T(-1) = 0 (Dirichlet);  k T' + h (T - Tinf) = 0 en x = 1 (Robin, h = 5, Tinf = 50)
%   Elemento 1 lineal [-1, 0], elemento 2 cuadratico [0, 1]

model.k   = 2;
model.c   = 0;
model.rho = 1;  model.cp = 1;        % no se usan en estacionario
model.G   = @(x) 10*x.^3;
model.tipoFF = [1 2];                % lineal + cuadratico (requiere la version mixta)

xnode = [-1, 0, 0.5, 1];             % 0.5 es el nodo central del cuadratico
cb = [1, 0, -1;                      % izquierda: Dirichlet T = 0
      3, 5, 50];                     % derecha:   Robin h = 5, Tinf = 50

T = elemFinitos(xnode, model, cb, 0);
disp('c) Temperaturas nodales [x, T]:'); disp([xnode', T]);

% d) Interpolacion con las funciones de forma del elemento que contiene a x
% x = -0.25 cae en el elemento 1 (lineal, nodos 1-2): psi = 0.5
psi = 0.5;   N = [(1-psi)/2, (1+psi)/2];
fprintf('d) T(-0.25) = %.4f\n', N * T([1 2]));
% x = 0.25 cae en el elemento 2 (cuadratico, nodos 2-3-4): psi = -0.5
psi = -0.5;  N = [psi*(psi-1)/2, 1-psi^2, psi*(psi+1)/2];
fprintf('   T( 0.25) = %.4f\n', N * T([2 3 4]));
