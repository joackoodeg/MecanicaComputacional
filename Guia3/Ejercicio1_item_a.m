% Ejercicio 1 (Guia 3), item a - prueba de elemFinitos.m
%
%   k*T'' - c*T + G = 0,  x en [L1,L2] = [0,1]
%   T(0) = 10 (Dirichlet),  T(1) = 50 (Dirichlet)
%   cp = 0, K = 2, c = 0, G = 100
%
%   Solucion analitica:  T(x) = -25 x^2 + 65 x + 10

graphics_toolkit('gnuplot');
set(0, 'defaultfigurevisible', 'off');

Tanalitica = @(x) -25*x.^2 + 65*x + 10;

model.k   = 2;
model.c   = 0;
model.rho = 1;  model.cp = 0;   % no se usan en estacionario
model.G   = 100;

cb = [1, 10, -1;    % izquierda (x=0): Dirichlet T=10
      1, 50, -1];   % derecha   (x=1): Dirichlet T=50

et = 0;             % estado estacionario

Nelem = 5;

% --- Funciones de forma lineales ---
model.tipoFF = 1;
xnode_lin = linspace(0, 1, Nelem+1);
T_lin = elemFinitos(xnode_lin, model, cb, et);

% --- Funciones de forma cuadraticas ---
model.tipoFF = 2;
xnode_cuad = linspace(0, 1, 2*Nelem+1);
T_cuad = elemFinitos(xnode_cuad, model, cb, et);

% --- Comparacion con la analitica ---
xfino = linspace(0, 1, 200);

fh = figure('visible', 'off');
graphics_toolkit(fh, 'gnuplot');
plot(xfino, Tanalitica(xfino), 'k-', 'linewidth', 1.5); hold on;
plot(xnode_lin,  T_lin,  'bo-', 'linewidth', 1.2);
plot(xnode_cuad, T_cuad, 'rs-', 'linewidth', 1.2);
legend('Analitica', 'FEM lineal', 'FEM cuadratica', 'location', 'best');
xlabel('x'); ylabel('T(x)');
title('Ejercicio 1.a: FEM vs analitica');
grid on;

err_lin  = norm(T_lin  - Tanalitica(xnode_lin)(:), Inf);
err_cuad = norm(T_cuad - Tanalitica(xnode_cuad)(:), Inf);

fprintf('\nError maximo (lineal,   %d elem): %.3e\n', Nelem, err_lin);
fprintf('Error maximo (cuadratica, %d elem): %.3e\n', Nelem, err_cuad);

print(fh, 'Ejercicio1_item_a.png', '-dpng', '-r120');
