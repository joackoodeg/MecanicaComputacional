clear all; close all; clc;

%% Datos del modelo (item a de la tabla)
x1 = 0; x2 = 4;
y1 = 0; y2 = 3;

k_val = 1;
c_val = 0;

%% Generacion de la malla (uniforme)
nx = 31;   % cantidad de nodos en x
ny = 41;   % cantidad de nodos en y

xcoord = linspace(x1, x2, nx);
ycoord = linspace(y1, y2, ny);

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
DIR = [];
NEU = [];
for j = 1:ny
    for i = 1:nx
        if j==ny                  % borde superior: T=50
            DIR = [DIR; node(i,j), 50];
        elseif i==1                % borde izquierdo: T=50
            DIR = [DIR; node(i,j), 50];
        elseif j==1                 % borde inferior: q=0 (Neumann, direccion 1=Sur)
            NEU = [NEU; node(i,j), 0, 1];
        elseif i==nx                  % borde derecho: q=0 (Neumann, direccion 2=Este)
            NEU = [NEU; node(i,j), 0, 2];
        end
    end
end
ROB = [];   % No hay condiciones Robin en el item (c)

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = 2*xnode(:,1).^2 + 3*xnode(:,2).^3;
model.ts = 2;   % distinto de 0 y de 1 => fdm2d_solve resuelve el caso estacionario

%% Resolucion (unicamente con las funciones fdm2d_*.m de la catedra)
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% Resultado en el centro de la placa
xc = 2.0; yc = 1.5;
[~, idx] = min((xnode(:,1)-xc).^2 + (xnode(:,2)-yc).^2);
printf("Temperatura en el centro (%.2f, %.2f): T = %.4f\n", xnode(idx,1), xnode(idx,2), PHI(idx));

%% Graficos
fdm2d_graph_mesh(PHI, Q, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(PHI, Q, xnode, icone, 0, 1);   % Flujo de calor (vectorial)
