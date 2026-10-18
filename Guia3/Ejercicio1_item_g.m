% Ejercicio 1 (Guia 3), item g - prueba de elemFinitos.m
%
%   rho*cp*dT/dt = k*T'' - c*T + G,  x en [L1,L2] = [0,1]
%   T(0) = 50 (Dirichlet);  q(1) = 5 (Neumann, flujo saliente no nulo)
%   cp = 1, K = 2, c = -2 (reaccion NEGATIVA), G = 0
%
%   Solucion analitica (estado estacionario):
%     T(x) = 73.2433*sin(x) + 50*cos(x)
%
% Este caso agrega, respecto a los anteriores: condicion de Neumann con
% flujo no nulo, y coeficiente de reaccion c < 0 (caso limite interesante
% para la estabilidad del esquema explicito Forward Euler).

graphics_toolkit('gnuplot');
set(0, 'defaultfigurevisible', 'off');

Tanalitica = @(x) 73.2433*sin(x) + 50*cos(x);

model.k   = 2;
model.c   = -2;
model.rho = 1;  model.cp = 1;
model.G   = 0;

cb = [1, 50, -1;   % izquierda (x=0): Dirichlet T=50
      2, 5,  -1];  % derecha   (x=1): Neumann q=5

Nelem = 10;
L1 = 0; L2 = 1;

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
% Malla fina (dx=0.1) con k=2, rho*cp=1 -> estabilidad explicita requiere
% dt <~ dx^2/(2*k/(rho*cp)) = 0.01/4 = 0.0025. Usamos dt=0.002.
model.tipoFF = 1;
dt     = 0.002;
nPasos = 4000;
tol    = 1e-8;

nombres = {'Forward Euler', 'Backward Euler', 'Crank-Nicolson'};
colores = {'r', 'g', 'm'};
T_final = zeros(Nelem+1, 3);

fprintf('\n--- Esquemas temporales (dt=%.3g, tol=%.1e) ---\n', dt, tol);
for s = 1:3
  et = [s, nPasos, tol, dt];
  Thist = elemFinitos(xnode_lin, model, cb, et);
  T_final(:,s) = Thist(:,end);
  err = norm(T_final(:,s) - Tanalitica(xnode_lin)(:), Inf);
  estado = 'OK';
  if any(isnan(T_final(:,s))) || any(isinf(T_final(:,s))) || err > 10
    estado = 'DIVERGIO';
  endif
  fprintf('%-16s pasos usados: %4d   error max vs analitica: %.3e   [%s]\n', ...
          nombres{s}, size(Thist,2)-1, err, estado);
endfor

% ===================== Graficos =====================
xfino = linspace(L1, L2, 200);

fh = figure('visible', 'off');
graphics_toolkit(fh, 'gnuplot');
plot(xfino, Tanalitica(xfino), 'k-', 'linewidth', 1.5); hold on;
plot(xnode_lin,  T_lin,  'bo-', 'linewidth', 1.2);
plot(xnode_cuad, T_cuad, 'cs-', 'linewidth', 1.2);
for s = 1:3
  if all(isfinite(T_final(:,s)))
    plot(xnode_lin, T_final(:,s), [colores{s} '+--'], 'linewidth', 1.0);
  endif
endfor
legend(['Analitica', 'FEM lineal (et=0)', 'FEM cuadratica (et=0)', nombres], ...
       'location', 'northeast');
xlabel('x'); ylabel('T(x)');
title('Ejercicio 1.g: FEM vs analitica (Neumann q!=0 + reaccion c<0)');
grid on;

print(fh, 'Ejercicio1_item_g.png', '-dpng', '-r120');
