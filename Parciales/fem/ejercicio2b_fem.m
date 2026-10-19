% EJERCICIO 2 (barra de 5 elementos) - FEM 1D
%   rho*cp*T_t = k*T'' + G,  x en [0, 0.5] m,  5 elementos lineales de h = 0.1 m
%   k = 25, rho*cp = 1, G = 3000 W/m solo en e5 (x entre 0.4 y 0.5)
%   dT/dx = 0.5 en x = 0  -> Neumann, q saliente = +k*dT/dx = 12.5
%   T = 100 en x = 0.5     -> Dirichlet
% Todo en metros: h = 0.1 m, k = 25 W/m C -> k/h = 250; G = 30 W/cm = 3000 W/m.
% Usa elemFinitos.m (version con T0 como 5to argumento)
clear; clc; close all;

model.k = 25;  model.c = 0;  model.rho = 1;  model.cp = 1;
model.G = @(x) 3000*(x > 0.4);              % fuente solo en el elemento 5
model.tipoFF = 1;
xnode = 0:0.1:0.5;  h = 0.1;
cb = [2  12.5 -1;                        % izquierda: Neumann, q = k*T'(0) = 12.5 (sale)
      1 100   -1];                       % derecha:   Dirichlet T = 100

%% a) Estacionario
T = elemFinitos(xnode, model, cb, 0);
disp('a) T nodal [x, T]:'); disp([xnode', T]);

%% b) x = 0.25 m cae en e3 (nodos 3 y 4, x = 0.2 y 0.3)
N3 = (0.3 - 0.25)/h;  N4 = (0.25 - 0.2)/h;
T25 = N3*T(3) + N4*T(4);
q25 = -model.k*(T(4) - T(3))/h;
printf('b) T(0.25) = %.3f,  q = -k dT/dx = %.3f\n', T25, q25);

%% c) Un paso implicito con dt = 0.05 a partir de T2, T3, T4
dt = 0.05;
% Fila 3 a mano: (h/(6dt))[1 4 1] + (k/h)[-1 2 -1] sobre (T2,T3,T4)^{n+1}
a  = h/(6*dt);  b = model.k/h;
Tn = [63.5 67.5 74.5];
rhs = a*(Tn(1) + 4*Tn(2) + Tn(3));       % F3 = 0 (no hay fuente en e2 ni e3)
% Con vecinos a n+1 aproximados por sus valores a n (unico dato disponible):
T3_fila = (rhs - (a - b)*(Tn(1) + Tn(3))) / (4*a + 2*b);
printf('c) Solo con la fila 3 (vecinos congelados): T3^{n+1} = %.3f\n', T3_fila);

% Sistema completo: necesita T1^n y T5^n (no dados). Se prueban valores razonables:
for T15 = [59.5 85; 63.5 74.5; 61.5 81.5]'
  T0 = [T15(1) 63.5 67.5 74.5 T15(2) 100];
  Tt = elemFinitos(xnode, model, cb, [2, 1, 1e-12, dt], T0);
  printf('   Sistema completo con T1^n = %.1f, T5^n = %.1f:  T3^{n+1} = %.3f\n', ...
         T15(1), T15(2), Tt(3,2));
end
