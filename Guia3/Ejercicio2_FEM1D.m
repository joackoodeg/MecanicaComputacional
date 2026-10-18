% EJERCICIO 2 - FEM 1D - Ecuacion del calor
%   k T'' + G = 0 en [-1, 1],  k = 2,  G(x) = 10 x^3
%   T(-1) = 0                                   (Dirichlet)
%   k dT/dx + h (T - Tinf) = 0 en x = 1          (Robin, h = 5, Tinf = 50)
%   Elemento 1 [-1, 0]: lineal.  Elemento 2 [0, 1]: cuadratico.
%
% Nodos globales: 1 (x=-1), 2 (x=0), 3 (x=0.5, central del cuadratico), 4 (x=1)
% Script autocontenido: no usa paquetes extra.

clear; clc;

k = 2;
G = @(x) 10*x.^3;
h_rob = 5;  T_inf = 50;  T_izq = 0;

xnode = [-1, 0, 0.5, 1];
icone = {[1 2], [2 3 4]};
grado = [1, 2];

% Elemento master [-1, 1]: funciones de forma
N_lin  = @(psi) [(1 - psi)/2, (1 + psi)/2];
N_cuad = @(psi) [psi.*(psi - 1)/2, (1 - psi).*(1 + psi), psi.*(psi + 1)/2];

% Gauss de 3 puntos: exacto hasta grado 5 (10x^3 por una cuadratica es grado 5)
pg = [-sqrt(15)/5, 0, sqrt(15)/5];
wg = [5/9, 8/9, 5/9];

%% a) Matrices y vectores elementales
K = zeros(4);  F = zeros(4, 1);
for e = 1:2
  nn = icone{e};
  xe = xnode(nn);
  he = xe(end) - xe(1);
  J  = he/2;
  xc = (xe(1) + xe(end))/2;
  if grado(e) == 1
    Ke = (k/he) * [1 -1; -1 1];
    Nf = N_lin;
  else
    Ke = (k/(3*he)) * [7 -8 1; -8 16 -8; 1 -8 7];
    Nf = N_cuad;
  end
  Fe = zeros(length(nn), 1);
  for g = 1:3
    Fe = Fe + wg(g) * G(xc + pg(g)*J) * Nf(pg(g))' * J;
  end
  fprintf('\na) Elemento %d (nodos %s, h = %g):\n', e, mat2str(nn), he);
  disp('   K^e ='); disp(Ke);
  disp('   F^e ='); disp(Fe);
  K(nn, nn) = K(nn, nn) + Ke;
  F(nn)     = F(nn) + Fe;
end

%% b) Matriz y vector globales
disp('b) K global (antes de condiciones de borde):'); disp(K);
disp('   F global (antes de condiciones de borde):'); disp(F);

% Robin en el nodo 4: K44 += h, F4 += h*Tinf
K(4, 4) = K(4, 4) + h_rob;
F(4)    = F(4) + h_rob * T_inf;
K_sinD = K;  F_sinD = F;                 % copia para el flujo reactivo

% Dirichlet en el nodo 1: se reemplaza la fila 1
K(1, :) = 0;  K(1, 1) = 1;  F(1) = T_izq;

disp('   K global (con condiciones de borde):'); disp(K);
disp('   F global (con condiciones de borde):'); disp(F);

%% c) Temperaturas
T = K \ F;
disp('c) Temperaturas nodales [x, T]:'); disp([xnode', T]);

%% d) Temperaturas en x = -0.25 (elemento 1) y x = 0.25 (elemento 2)
% x -> psi: psi = (2x - x1 - x2) / (x2 - x1)
psi1 = (2*(-0.25) - (-1) - 0) / 1;       % = 0.5
psi2 = (2*( 0.25) -   0  - 1) / 1;       % = -0.5
T_m025 = N_lin(psi1)  * T([1 2]);
T_p025 = N_cuad(psi2) * T([2 3 4]);
fprintf('d) T(-0.25) = %.4f   (N = %s)\n', T_m025, mat2str(N_lin(psi1), 4));
fprintf('   T( 0.25) = %.4f   (N = %s)\n', T_p025, mat2str(N_cuad(psi2), 4));

%% Extra: flujo reactivo y balance
Q_izq = K_sinD(1, :) * T - F_sinD(1);    % calor que ENTRA por x = -1
Q_rob = h_rob * (T_inf - T(4));          % calor que ENTRA por x = 1
fprintf('\nExtra: calor que entra por x = -1: %g  (negativo = sale)\n', Q_izq);
fprintf('       calor que entra por x =  1: %g\n', Q_rob);
fprintf('       fuente total: %g  -> balance: %g (debe ser 0)\n', sum(F_sinD) - h_rob*T_inf, ...
        Q_izq + Q_rob + (sum(F_sinD) - h_rob*T_inf));
