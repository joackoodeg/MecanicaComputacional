% ejercicio_1.m
%
% Guia de Trabajos Practicos N2 - Metodo de Diferencias Finitas 2D
% Ejercicio 1 - item (a)
%
% Placa cuadrada [0,1]x[0,1] con T=20 impuesta en los 4 bordes.
% Datos del modelo (tabla del enunciado): rho*cp=0, k=2, c=0, G=300.
%
% Como rho*cp=0, el termino temporal de la ecuacion se anula y el problema
% queda definido exclusivamente por su version estacionaria:
%   0 = k*(d2T/dx2 + d2T/dy2) - c*T + G
% Por esta razon no se aplican los esquemas Forward/Backward Euler en este
% item (dividen por rho*cp en fdm2d_explicit_delta_t.m y fdm2d_implicit.m,
% lo cual no tiene sentido con rho*cp=0): se resuelve con el modo
% estacionario de fdm2d_solve.m.
%
% El script usa exclusivamente las funciones fdm2d_*.m de la catedra;
% solo arma la malla y las condiciones de borde (xnode, icone, DIR, NEU,
% ROB, model) que esas funciones reciben como entrada.

clear all; close all; clc;

%% Datos del modelo (item a de la tabla)
x1 = 0; x2 = 1;
y1 = 0; y2 = 1;

k_val = 2;
c_val = 0;
G_val = 300;
Tborde = 20;

%% Generacion de la malla (uniforme)
nx = 21;   % cantidad de nodos en x
ny = 21;   % cantidad de nodos en y

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

%% Condiciones de borde: Dirichlet T=20 en los 4 lados de la placa
DIR = [];
for j = 1:ny
    for i = 1:nx
        if (i==1 || i==nx || j==1 || j==ny)
            DIR = [DIR; node(i,j), Tborde];
        end
    end
end

NEU = [];   % No hay condiciones Neumann en el item (a)
ROB = [];   % No hay condiciones Robin en el item (a)

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = G_val * ones(nnodes,1);
model.ts = 2;   % distinto de 0 y de 1 => fdm2d_solve resuelve el caso estacionario

%% Resolucion (unicamente con las funciones fdm2d_*.m de la catedra)
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% Resultado en el centro de la placa
xc = 0.5; yc = 0.5;
[~, idx] = min((xnode(:,1)-xc).^2 + (xnode(:,2)-yc).^2);
printf("Temperatura en el centro (%.2f, %.2f): T = %.4f\n", xnode(idx,1), xnode(idx,2), PHI(idx));

%% Graficos
fdm2d_graph_mesh(PHI, Q, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(PHI, Q, xnode, icone, 0, 1);   % Flujo de calor (vectorial)
