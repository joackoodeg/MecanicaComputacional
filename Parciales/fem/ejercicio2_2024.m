% EJERCICIO 2 - FEM 1D transitorio
%   rho*cp*T_t = k*T'' + G(x),  x en [2, 4],  rho*cp = 1, k = 2, G = 100*(4 - x)
%   T(2) = 10 (Dirichlet);  -k dT/dx = 1 en x = 4 (Neumann, q = 1 sale)
%   Malla uniforme dx = 0.5 -> nodos 2, 2.5, 3, 3.5, 4
% Usa elemFinitos.m (version con T0 como 5to argumento)
clear; clc; close all;

model.k = 2;  model.c = 0;  model.rho = 1;  model.cp = 1;
model.G = @(x) 100*(4 - x);
model.tipoFF = 1;                        % lineales

xnode = 2:0.5:4;
cb = [1 10 -1;                           % izquierda: Dirichlet T = 10
      2  1 -1];                          % derecha:   Neumann q = 1 (sale)
h = 0.5;  N = length(xnode);

%% a) dt critico del explicito (con las mismas matrices que arma elemFinitos)
K = zeros(N);  C = zeros(N);
for e = 1:N-1
  n = [e e+1];
  K(n,n) += (model.k/h)*[1 -1; -1 1];
  C(n,n) += (model.rho*model.cp*h/6)*[2 1; 1 2];
end
lib = 2:N;                               % el nodo 1 es Dirichlet
lmax = max(eig(C(lib,lib) \ K(lib,lib)));
dtc  = 2/lmax;
printf('a) lambda_max = %.4f  ->  dt_crit = 2/lambda_max = %.5f\n', lmax, dtc);
printf('   estimacion rho*cp*h^2/(6k) = %.5f\n', model.rho*model.cp*h^2/(6*model.k));

%% b) Un paso implicito (Backward Euler) con dt = 0.1 desde el T^n dado
Tn = [10 39.848 52.149 54.173 52.852];
dt = 0.1;
T  = elemFinitos(xnode, model, cb, [2, 1, 1e-12, dt], Tn);
Tn1 = T(:,2);
disp('b) T^{n+1} [x, T]:'); disp([xnode', Tn1]);

%% c) Interpolacion en x = 3.2 (elemento 3: nodos 3 y 4, x = 3 y 3.5)
xp = 3.2;  e = 3;  xa = xnode(e);  xb = xnode(e+1);
Na = (xb - xp)/h;  Nb = (xp - xa)/h;
Tp = Na*Tn1(e) + Nb*Tn1(e+1);
printf('c) T(3.2) = %.4f*%.4f + %.4f*%.4f = %.4f\n', Na, Tn1(e), Nb, Tn1(e+1), Tp);

%% d) Flujo en x = 3.2: q = -k dT/dx (constante en el elemento lineal)
dTdx = (Tn1(e+1) - Tn1(e))/h;
q = -model.k*dTdx;
if q > 0, dir = 'hacia +x'; else, dir = 'hacia -x'; end
printf('d) dT/dx = %.4f  ->  q = -k dT/dx = %.4f  (%s)\n', dTdx, q, dir);
