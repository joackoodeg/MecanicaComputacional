% Practica - Problema 2: placa [-4,4]x[-3,3], MDF 2D, estacionario
%   k = 3,  G = 100*x  (fuente NO uniforme, varia linealmente con x), c = 0
%   Izquierdo (x=-4): T = 20              (Dirichlet)
%   Derecho   (x=4):  Robin h = 2, Tinf = 10
%   Superior  (y=3):  q = 0               (Neumann, aislado)
%   Inferior  (y=-3): q = 0               (Neumann, aislado)
%   Malla dada: dx = dy = 2  ->  5x4 nodos (20 nodos, 12 elementos)
clear all; close all; clc;

%% Datos del modelo
x1 = -4; x2 = 4;
y1 = -3; y2 = 3;

k_val = 3;
c_val = 0;      % sin reaccion

%% ========================================================================
%% a) Malla 2D del enunciado (dx=dy=2)
%% ========================================================================
nx = 5;   % cantidad de nodos en x (-4,-2,0,2,4)
ny = 4;   % cantidad de nodos en y (-3,-1,1,3)

xcoord = linspace(x1, x2, nx);
ycoord = linspace(y1, y2, ny);
h = xcoord(2) - xcoord(1);

nnodes = nx*ny;
xnode = zeros(nnodes, 2);

node = @(i,j) (j-1)*nx + i;   % numeracion de nodos: recorre x y luego y

for j = 1:ny
    for i = 1:nx
        xnode(node(i,j), :) = [xcoord(i), ycoord(j)];
    end
end

%% Conectividad: elementos cuadrangulares (P1 inf-izq, P2 inf-der, P3 sup-der, P4 sup-izq)
nelem = (nx-1)*(ny-1);
icone = zeros(nelem, 4);

e = 1;
for j = 1:ny-1
    for i = 1:nx-1
        P1 = node(i,   j);
        P2 = node(i+1, j);
        P3 = node(i+1, j+1);
        P4 = node(i,   j+1);
        icone(e,:) = [P1 P2 P3 P4];
        e = e+1;
    end
end

%% Condiciones de borde
% Esquinas izquierdas (-4,-3) y (-4,3): Dirichlet manda sobre el aislamiento
%   (Neumann) de arriba/abajo. Esquinas derechas (4,-3) y (4,3): Robin +
%   Neumann se suman (ninguna de las dos manda sobre la otra).
DIR = []; NEU = []; ROB = [];
for j = 1:ny
    DIR = [DIR; node(1,j), 20];                     % izquierdo (incluye ambas esquinas)
end
for j = 1:ny
    ROB = [ROB; node(nx,j), 2, 10, 2];               % derecho: h=2, Tinf=10, direccion 2 = Este
end
for i = 2:nx
    NEU = [NEU; node(i,ny), 0, 3];                   % superior: q=0, direccion 3 = Norte
    NEU = [NEU; node(i,1),  0, 1];                   % inferior: q=0, direccion 1 = Sur
end

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = 100 * xnode(:,1);      % fuente G(x,y) = 100*x, no uniforme
model.ts = 'X';                  % [0] explicito, [1] implicito, [X] estacionario

%% Resolver
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% a) Temperatura en (0,0)
% (0,0) no es nodo: x=0 coincide con la columna i=3, pero y=0 cae justo
% entre y=-1 (j=2) e y=1 (j=3) -> como es el punto medio, alcanza con
% promediar esos dos nodos (interpolacion lineal en y).
T00_2D = (PHI(node(3,2)) + PHI(node(3,3))) / 2;
printf('\na) T(0,0) [malla 2D, %dx%d nodos] = %.4f\n', nx, ny, T00_2D);

%% Graficos
fdm2d_graph_mesh(full(PHI), Q, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(full(PHI), Q, xnode, icone, 0, 1);   % Flujo de calor (vectorial)

%% ========================================================================
%% b) Reduccion a 1D por simetria (21 nodos)
%% ========================================================================
% G = 100*x no depende de y, y las condiciones de borde superior e inferior
% son ambas de aislamiento (q=0), tambien independientes de y. Por lo tanto
% T no depende de y y el problema es equivalente a uno unidimensional a lo
% largo del eje x, con las mismas condiciones Dirichlet/Robin en los extremos.
N1 = 21;
x1d = linspace(x1, x2, N1);
dx1d = x1d(2) - x1d(1);

m.rho = 1; m.cp = 1;            % no se usan en el esquema estacionario
m.k = k_val;
m.c = c_val;
m.G = 100 * x1d;                % fuente G(x) = 100*x en cada nodo

cb = [1 20 -1;     % izquierdo: tipo 1 = Dirichlet, T = 20
      3 2  10];    % derecho:   tipo 3 = Robin, h = 2, Tinf = 10

T1d = difFinitas(x1d, m, cb, [0 0 0 0], zeros(1,N1));

i0 = find(abs(x1d) < 1e-9, 1);         % nodo x = 0 (nodo 11 de 21)
T00_1D = T1d(i0);
printf('\nb) T(0,0) [modelo 1D, %d nodos] = %.4f  (coincide con el 2D por simetria)\n', N1, T00_1D);

%% ========================================================================
%% c) Gradiente de temperatura en (4,0), aproximacion de 2do orden
%% ========================================================================
% (4,0) es el ultimo nodo de la malla 1D (borde derecho, x = x1d(21)): no hay
% nodo a su derecha, por lo que se usa la formula de diferencias finitas de
% 2do orden NO centrada (3 puntos hacia atras):
%   dT/dx|_N = (3*T_N - 4*T_(N-1) + T_(N-2)) / (2*dx)
dTdx_4 = (3*T1d(N1) - 4*T1d(N1-1) + T1d(N1-2)) / (2*dx1d);
printf('\nc) dT/dx(4,0) = %.4f   (dT/dy = 0 por simetria -> grad T = [%.4f, 0])\n', dTdx_4, dTdx_4);
printf('   q(4,0) = -k*dT/dx = %.4f\n', -k_val*dTdx_4);
