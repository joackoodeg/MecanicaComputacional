% Practica - Transferencia de calor estacionaria con fuente G = 100*x
% Dominio: x en [-4, 4], y en [-3, 3]
% Datos: k = 3, G = 100*x, c = 0 (sin reaccion), estacionario
%   Superior (y=3):   aislado q = 0           (Neumann, dir 3 = Norte)
%   Inferior (y=-3):  aislado q = 0           (Neumann, dir 1 = Sur)
%   Izquierdo (x=-4): T = 20                  (Dirichlet)
%   Derecho (x=4):    Robin h = 2, Tinf = 10  (Robin, dir 2 = Este)
%   dx = dy = 2
clear all; close all; clc;

%% ========================================================================
%% INCISO a) Resolucion 2D con la malla propuesta (dx = dy = 2)
%% ========================================================================

%% Datos del modelo
x1 = -4; x2 = 4;
y1 = -3; y2 = 3;

k_val = 3;
c_val = 0;      % sin reaccion

%% Generacion de la malla (uniforme, coincide con xnode e icone del enunciado)
nx = 5;   % cantidad de nodos en x (-4, -2, 0, 2, 4) -> dx = 2
ny = 4;   % cantidad de nodos en y (-3, -1, 1, 3)    -> dy = 2

xcoord = linspace(x1, x2, nx);
ycoord = linspace(y1, y2, ny);
dx = xcoord(2) - xcoord(1);
dy = ycoord(2) - ycoord(1);

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
% Esquinas izquierdas: (-4,-3) [nodo 1] y (-4,3) [nodo 16] -> T = 20 (Dirichlet manda)
% Esquinas derechas:   (4,-3)  [nodo 5] y (4,3)  [nodo 20] -> Robin der + Neumann abajo/arriba (se suman)
DIR = []; NEU = []; ROB = [];

for j = 1:ny
    DIR = [DIR; node(1,j), 20];                      % izquierdo: T = 20 (incluye ambas esquinas)
end
for i = 2:nx
    NEU = [NEU; node(i,1), 0, 1];                    % inferior: q = 0, direccion 1 = Sur
end
for i = 2:nx
    NEU = [NEU; node(i,ny), 0, 3];                   % superior: q = 0, direccion 3 = Norte
end
for j = 1:ny
    ROB = [ROB; node(nx,j), 2, 10, 2];               % derecho: h = 2, Tinf = 10, direccion 2 = Este
end

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = 100 * xnode(:,1);     % fuente G = 100*x
model.ts = -1;                  % [0] explicito, [1] implicito, [X] estacionario

%% Resolver 2D
[Tf, Qf] = fdm2d(xnode, icone, DIR, NEU, ROB, model);
% Mirando la tabla xnode: el punto (0,0) esta justo en el medio
% entre el nodo 8 (0,-1) y el nodo 13 (0,1), asi que promediamos:
T00_2D = (Tf(8) + Tf(13)) / 2;
printf('\na) Temperatura en (0,0) [Malla 2D]:\n');
printf('   Nodo 8  (0,-1) = %.4f\n', Tf(8));
printf('   Nodo 13 (0, 1) = %.4f\n', Tf(13));
printf('   T(0,0)         = %.4f\n', T00_2D);

%% Graficos 2D
fdm2d_graph_mesh(full(Tf), Qf, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(full(Tf), Qf, xnode, icone, 0, 1);   % Flujo de calor (vectorial)


%% ========================================================================
%% INCISO b) Resolucion 1D con 21 nodos totales usando difFinitas.m
%% ========================================================================
N1d = 21;
x1d = linspace(x1, x2, N1d);    % 21 nodos -> dx = 8/20 = 0.4
dx1d = x1d(2) - x1d(1);

m1d.rho = 1; m1d.cp = 1; m1d.k = k_val; m1d.c = c_val;
m1d.G   = 100 * x1d;            % G(x) = 100*x en cada nodo

cb1d = [1, 20, -1;              % izquierda (x=-4): tipo 1 = Dirichlet, T = 20
        3,  2, 10];             % derecha   (x=4):  tipo 3 = Robin, h = 2, Tinf = 10

figure;
T1d = difFinitas(x1d, m1d, cb1d, [0 0 0 0], zeros(1, N1d));
% Entre 21 nodos de -4 a 4, el punto central (x = 0) es el nodo 11:
T00_1D = T1d(11);
printf('\nb) Temperatura en (0,0) [Modelo 1D con 21 nodos]: %.4f\n', T00_1D);


%% ========================================================================
%% INCISO c) Gradiente de temperatura en (4,0) de segundo orden
%% ========================================================================
% 1) Por Taylor hacia atras de 3 puntos (orden dx^2):
grad_x_3pt = (3*T1d(N1d) - 4*T1d(N1d-1) + T1d(N1d-2)) / (2*dx1d);

% 2) Por diferencia centrada con el nodo ficticio T_fic (o de la CB Robin):
T_fic = 2*T1d(N1d) - T1d(N1d-1) - (dx1d^2 / k_val)*m1d.G(N1d);
grad_x_fic = (T_fic - T1d(N1d-1)) / (2*dx1d);

printf('\nc) Gradiente de temperatura en (4,0) [2do orden, O(dx^2)]:\n');
printf('   - Con 3 puntos backward:  dT/dx = %.4f  -> grad(T) = [%.4f, 0]\n', grad_x_3pt, grad_x_3pt);
printf('   - Con nodo ficticio/Robin: dT/dx = %.4f  -> grad(T) = [%.4f, 0]\n', grad_x_fic, grad_x_fic);