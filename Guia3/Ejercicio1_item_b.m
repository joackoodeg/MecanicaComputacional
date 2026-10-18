% Ejercicio 1 (Guia 3), item b - prueba de elemFinitos.m
%
%   k*T'' - c*T + G = 0,  x en [L1,L2] = [0,2]
%   T(0) = 100 (Dirichlet),  q(2) = 0 (Neumann)
%   cp = 0, K = 1, c = 1, G = 0
%
%   Solucion analitica:  T(x) = 100*exp(-x).*(exp(2x)+exp(4)) / (1+exp(4))
%
% Este caso agrega, respecto al item a: condicion de Neumann y termino
% de reaccion (c != 0).

graphics_toolkit('gnuplot');
set(0, 'defaultfigurevisible', 'off');

Tanalitica = @(x) 100*exp(-x).*(exp(2*x) + exp(4)) / (1 + exp(4));

model.k   = 1;
model.c   = 1;
model.rho = 1;  model.cp = 0;   % no se usan en estacionario
model.G   = 0;

cb = [1, 100, -1;   % izquierda (x=0): Dirichlet T=100
      2, 0,   -1];  % derecha   (x=2): Neumann q=0

et = 0;             % estado estacionario

Nelem = 10;

% --- Funciones de forma lineales ---
model.tipoFF = 1;
xnode_lin = linspace(0, 2, Nelem+1);
T_lin = elemFinitos(xnode_lin, model, cb, et);

% --- Funciones de forma cuadraticas ---
model.tipoFF = 2;
xnode_cuad = linspace(0, 2, 2*Nelem+1);
T_cuad = elemFinitos(xnode_cuad, model, cb, et);

% --- Comparacion con la analitica ---
xfino = linspace(0, 2, 200);

fh = figure('visible', 'off');
graphics_toolkit(fh, 'gnuplot');
plot(xfino, Tanalitica(xfino), 'k-', 'linewidth', 1.5); hold on;
plot(xnode_lin,  T_lin,  'bo-', 'linewidth', 1.2);
plot(xnode_cuad, T_cuad, 'rs-', 'linewidth', 1.2);
legend('Analitica', 'FEM lineal', 'FEM cuadratica', 'location', 'northeast');
xlabel('x'); ylabel('T(x)');
title('Ejercicio 1.b: FEM vs analitica (Neumann + reaccion)');
grid on;

err_lin  = norm(T_lin  - Tanalitica(xnode_lin)(:), Inf);
err_cuad = norm(T_cuad - Tanalitica(xnode_cuad)(:), Inf);

fprintf('\nError maximo (lineal,    %d elem): %.3e\n', Nelem, err_lin);
fprintf('Error maximo (cuadratica, %d elem): %.3e\n', Nelem, err_cuad);

print(fh, 'Ejercicio1_item_b.png', '-dpng', '-r120');
