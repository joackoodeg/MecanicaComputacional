clear all; close all; clc;

%% Datos del modelo (item a de la tabla)
x1 = 0; x2 = 3;
y1 = 0; y2 = 5;

k_val = 2;
c_val = 0;
G_val = 0;

%% Generacion de la malla (uniforme)
nx = 31;   % cantidad de nodos en x
ny = 31;   % cantidad de nodos en y

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
ROB = [];
for j = 1:ny
    for i = 1:nx
        if j==1                       % borde inferior
            ROB = [ROB; node(i,j), 2, 100, 1];
        elseif j==ny                   % borde superior:
            DIR = [DIR; node(i,j), 50];
        elseif i==1                     % borde izquierdo
            DIR = [DIR; node(i,j), 0];
        elseif i==nx                     % borde derecho: q=10 (Neumann, direccion 2=Este)
            DIR = [DIR; node(i,j), 25];
        end
    end
end

%% Resultado en el centro de la placa (indice del nodo mas cercano al centro)
xc = (x1+x2)/2;
yc = (y1+y2)/2;
[~, idx] = min((xnode(:,1)-xc).^2 + (xnode(:,2)-yc).^2);

%% Datos del modelo (vectores por nodo + esquema de resolucion)
model.nnodes = nnodes;
model.k = k_val * ones(nnodes,1);
model.c = c_val * ones(nnodes,1);
model.G = 10*xnode(:,1).^3;

%% Datos para los esquemas temporales (necesarios porque rho*cp=1 en este item)
model.rho = 1;
model.cp  = 1;          % rho*cp = 1, segun la tabla
model.PHI_n = zeros(nnodes,1);   % condicion inicial T(x,0)=0
model.tol = 1e-5;
model.maxit = 5000;

%% Resolucion: estacionario (el que ya tenias)
model.ts = 2;
[PHI_est, Q_est] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% Resolucion: Forward Euler (explicito)
model.ts = 0;
[PHI_exp, Q_exp] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% Resolucion: Backward Euler (implicito)
model.ts = 1;
model.dt = 0.05;
[PHI_imp, Q_imp] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

%% Comparacion en el centro de la placa
printf("Estacionario  : T = %.4f\n", PHI_est(idx));
printf("Forward Euler : T = %.4f (ultima iteracion)\n", PHI_exp(idx,end));
printf("Backward Euler: T = %.4f (ultima iteracion)\n", PHI_imp(idx,end));

%% Graficos (evolucion temporal del esquema explicito)
fdm2d_graph_mesh(PHI_exp, Q_exp, xnode, icone, 0, 0);   % Temperatura (escalar)
fdm2d_graph_mesh(PHI_exp, Q_exp, xnode, icone, 0, 1);   % Flujo de calor (vectorial)
