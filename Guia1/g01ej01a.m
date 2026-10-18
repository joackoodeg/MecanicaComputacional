# Guia 1 - Ejercicio 1 - Item a
# rho*cp*dT/dt = k*T'' - c*T + G   en [0,1]
# T(0)=10, T(1)=50, rho*cp=0, k=2, c=0, G=100
# Solucion analitica: T(x) = -25x^2 + 65x + 10
#
# rho*cp = 0  -> solo estado estacionario (sin esquemas temporales)
# Ambas CB Dirichlet -> no aplica el analisis de orden de las CB

clear all; close all; clc;

L1 = 0; L2 = 1;
T_ex = @(x) -25*x.^2 + 65*x + 10;

function [T, err] = resolver_a(xnode, T_ex)
  N = length(xnode);
  model.rho = 0;
  model.cp  = 0;
  model.k   = 2;
  model.c   = 0;
  model.G   = 100*ones(1,N);
  cb = [1 10 -1; 1 50 -1];   # Dirichlet en ambos extremos
  et = [0 0 0 0];            # estado estacionario
  T = tp_dif_finitas(xnode, model, cb, et);
  err = max(abs(T(:) - T_ex(xnode(:))));
endfunction

# --- 1) Resolucion con malla uniforme y comparacion con la analitica ---
N = 11;
xnode = linspace(L1, L2, N);
[T, err] = resolver_a(xnode, T_ex);

figure(1); clf; hold on;
xfino = linspace(L1, L2, 200);
plot(xfino, T_ex(xfino), 'r-', 'linewidth', 1.5);
plot(xnode, T, 'bx', 'markersize', 8);
legend('Analitica', 'Dif. finitas', 'location', 'northwest');
xlabel('x'); ylabel('T'); title('Ej 1a - estacionario, N=11'); grid on;

printf('\nItem a - N=%d uniforme: max|T - T_ex| = %.3e\n\n', N, err);

# --- 2) Malla uniforme vs no uniforme ---
Ns = [5 10 20 40 80];
figure(2);  # tp_dif_finitas plotea internamente; lo mandamos a una figura aparte
printf('%5s %12s %12s %12s %12s\n', 'N', 'h_unif', 'err_unif', 'h_nounif', 'err_nounif');
for N = Ns
  xu = linspace(L1, L2, N);                     # uniforme
  xn = L1 + (L2-L1)*linspace(0, 1, N).^2;       # no uniforme (refinada a la izquierda)
  [~, eu] = resolver_a(xu, T_ex);
  [~, en] = resolver_a(xn, T_ex);
  printf('%5d %12.4e %12.4e %12.4e %12.4e\n', N, max(diff(xu)), eu, max(diff(xn)), en);
endfor
close(2);
# El error es de redondeo en todos los casos: la solucion exacta es cuadratica
# y el stencil de 3 puntos (uniforme o no) es exacto para polinomios de grado 2.
