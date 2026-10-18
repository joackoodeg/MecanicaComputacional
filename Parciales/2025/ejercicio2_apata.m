% Parcial 15/09/2025 - Problema 2: placa 1x1, MDF 2D, esquema explicito
%   k = 10, G = 500 (uniforme), c = 0 (sin reaccion), rho*cp = 1
%   Superior (y=1): T = 150        (Dirichlet)
%   Inferior (y=0): Robin h = 50, Tinf = 20
%   Izquierdo (x=0): T = 100       (Dirichlet)
%   Derecho (x=1): aislado         (Neumann q = 0)
%   dx = dy = 1/6, T inicial = 125, tolerancia 1e-5
clear all; close all; clc;

%% Datos del modelo
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;

k_val = 10;
c_val = 0;      % sin reaccion
G_val = 500;    % fuente uniforme

%% Generacion de la malla (uniforme)
nx = 7;   % cantidad de nodos en x
ny = 7;   % cantidad de nodos en y

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
% Esquinas: (0,1) y (1,1) -> 150 (superior, Dirichlet manda);
%           (0,0) -> 100 (izquierdo, Dirichlet manda sobre Robin);
%           (1,0) -> Robin abajo + Neumann a la derecha (se suman, no se pisan)
DIR = []; NEU = []; ROB = [];
for i = 1:nx
    DIR = [DIR; node(i,ny), 150];                  % superior (incluye ambas esquinas)
end
for j = 1:ny-1
    DIR = [DIR; node(1,j), 100];                    % izquierdo (sin la esquina superior)
end
for i = 2:nx
    ROB = [ROB; node(i,1), 50, 20, 1];               % inferior: h, Tinf, direccion 1 = Sur
end
for j = 1:ny-1
    NEU = [NEU; node(nx,j), 0, 2];                   % derecho: q = 0, direccion 2 = Este
end

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = G_val * ones(nnodes,1);
model.rho = 1;
model.cp  = 1;
model.ts    = 0;          % [0] explicito, [1] implicito, [X] estacionario
model.maxit = 1e6;
model.tol   = 1e-5;

% Condicion inicial: 125 uniforme, con los nodos Dirichlet ya en su valor
% (si no, el codigo los lleva a su valor de a poco y falsea el tiempo).
model.PHI_n = 125*ones(nnodes,1);
model.PHI_n(DIR(:,1)) = DIR(:,2);

%% Resolver
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);
dt = fdm2d_explicit_delta_t(xnode, model);
nsteps = columns(PHI) - 1;
Tf = PHI(:,end);
Qf = Q(:, end-1:end);                               % flujo del ultimo paso [Qx Qy]

%% 1) Tiempo para alcanzar el estacionario
printf('\n1) dt = %g, pasos = %d, t estacionario = %.4f s\n', dt, nsteps, nsteps*dt);

%% 2) Temperaturas
% (0.5,0.5) y (1,0) son nodos; (0.25,0.25) NO (cae entre nodos) -> interpolacion bilineal
TP = @(x,y) Tf(node(round(x/h)+1, round(y/h)+1));
printf('2) T(0.5,0.5) = %.4f\n', TP(0.5,0.5));
printf('   T(1,0)     = %.4f\n', TP(1,0));
Tgrid = reshape(Tf, nx, ny);                        % Tgrid(i,j)
T025 = interp2(xcoord, ycoord, Tgrid', 0.25, 0.25, 'linear');
printf('   T(0.25,0.25) = %.4f  (interpolacion bilineal)\n', T025);

%% 3) Flujo en el centro
P = node(round(0.5/h)+1, round(0.5/h)+1);
printf('3) q(0.5,0.5) = [%.4f, %.4f], |q| = %.4f\n', Qf(P,1), Qf(P,2), norm(Qf(P,:)));

%% Graficos
fdm2d_graph_mesh(full(Tf), Qf, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(full(Tf), Qf, xnode, icone, 0, 1);   % Flujo de calor (vectorial)
