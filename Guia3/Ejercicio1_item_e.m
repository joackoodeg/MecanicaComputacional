% Ejercicio 1 (Guia 3), item e - prueba de elemFinitos.m
%
%   rho*cp*dT/dt = k*T'' - c*T + G,  x en [L1,L2] = [5,10]
%   Robin h=2, Tinf=100 en x=5;  T(10) = 50 (Dirichlet)
%   cp = 1, K = 2, c = 0, G = x^3
%
%   Solucion analitica (estado estacionario):
%     T(x) = -x^5/40 + 1225*x/3 - 4600/3
%
% Este caso agrega, respecto a los items a/b: condicion de Robin, fuente
% G(x) variable (no constante), y rho*cp != 0 -> permite probar los tres
% esquemas temporales (Forward Euler, Backward Euler, Crank-Nicolson)
% convergiendo al mismo estacionario.

graphics_toolkit('gnuplot');
set(0, 'defaultfigurevisible', 'off');

Tanalitica = @(x) -x.^5/40 + 1225*x/3 - 4600/3;

model.k   = 2;
model.c   = 0;
model.rho = 1;  model.cp = 1;
model.G   = @(x) x.^3;

cb = [3, 2, 100;   % izquierda (x=5):  Robin h=2, Tinf=100
      1, 50, -1];  % derecha   (x=10): Dirichlet T=50

Nelem = 10;
L1 = 5; L2 = 10;

% ===================== 1) Chequeo espacial (estacionario) =====================
model.tipoFF = 1;
xnode_lin = linspace(L1, L2, Nelem+1);
T_lin = elemFinitos(xnode_lin, model, cb, 0);

model.tipoFF = 2;
xnode_cuad = linspace(L1, L2, 2*Nelem+1);
T_cuad = elemFinitos(xnode_cuad, model, cb, 0);

err_lin  = norm(T_lin  - Tanalitica(xnode_lin)(:), Inf);
err_cuad = norm(T_cuad - Tanalitica(xnode_cuad)(:), Inf);
fprintf('--- Estacionario (et=0) ---\n');
fprintf('Error maximo (lineal,    %d elem): %.3e\n', Nelem, err_lin);
fprintf('Error maximo (cuadratica, %d elem): %.3e\n', Nelem, err_cuad);

% ===================== 2) Chequeo temporal (3 esquemas, lineal) ===============
model.tipoFF = 1;
dts    = [0.01, 0.05, 0.05];   % Forward Euler necesita dt mas chico (estabilidad)
nPasos = [3000, 2000, 2000];
tol    = 1e-8;

nombres = {'Forward Euler', 'Backward Euler', 'Crank-Nicolson'};
colores = {'r', 'g', 'm'};
T_final = zeros(Nelem+1, 3);

fprintf('\n--- Esquemas temporales (tol=%.1e) ---\n', tol);
for s = 1:3
  dt = dts(s);
  et = [s, nPasos(s), tol, dt];
  Thist = elemFinitos(xnode_lin, model, cb, et);
  T_final(:,s) = Thist(:,end);
  err = norm(T_final(:,s) - Tanalitica(xnode_lin)(:), Inf);
  fprintf('%-16s dt=%.3g  pasos usados: %4d   error max vs analitica: %.3e\n', ...
          nombres{s}, dt, size(Thist,2)-1, err);
endfor

% ===================== Graficos =====================
xfino = linspace(L1, L2, 200);

fh = figure('visible', 'off');
graphics_toolkit(fh, 'gnuplot');
plot(xfino, Tanalitica(xfino), 'k-', 'linewidth', 1.5); hold on;
plot(xnode_lin,  T_lin,  'bo-', 'linewidth', 1.2);
plot(xnode_cuad, T_cuad, 'cs-', 'linewidth', 1.2);
for s = 1:3
  plot(xnode_lin, T_final(:,s), [colores{s} '+--'], 'linewidth', 1.0);
endfor
legend(['Analitica', 'FEM lineal (et=0)', 'FEM cuadratica (et=0)', nombres], ...
       'location', 'northeast');
xlabel('x'); ylabel('T(x)');
title('Ejercicio 1.e: FEM vs analitica (Robin + G(x) + esquemas temporales)');
grid on;

print(fh, 'Ejercicio1_item_e.png', '-dpng', '-r120');
