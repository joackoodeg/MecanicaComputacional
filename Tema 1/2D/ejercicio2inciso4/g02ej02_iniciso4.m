close all; clear all; more off;

dx=1/10;
dy=1/10;

% Inicializar la matriz de nodos
xnode = [];

for j = 0:10  % Iterar sobre y
    for i = 0:10  % Iterar sobre x
        xnode = [xnode; i*dx, j*dy];
    end
end

plot(xnode(:,1),xnode(:,2),'o');

NEU = [];
ROB = [];
DIR = [
  1,0;
  2,0;
  3,0;
  4,0;
  5,0;
  6,0;
  7,0;
  8,0;
  9,0;
  10,0;
  11,0;
  12,0;
  22,0;
  23,0;
  33,0;
  34,0;
  44,0;
  45,0;
  55,0;
  56,0;
  66,0;
  67,0;
  77,0;
  78,0;
  88,0;
  89,0;
  99,0;
  100,0;
  110,0;
  111,0;
  121,0;
  ];

nx = 11;  % Número de nodos en x
ny = 11;  % Número de nodos en y
% Inicializar la matriz icone
icone = [];

% Crear los elementos
for j = 1:(ny-1)
    for i = 1:(nx-1)
        % Definir los nodos de cada cuadrado (elemento)
        n1 = i + (j-1)*nx;       % Nodo inferior izquierdo
        n2 = n1 + 1;             % Nodo inferior derecho
        n3 = n1 + nx;            % Nodo superior izquierdo
        n4 = n3 + 1;             % Nodo superior derecho
        % Añadir el elemento a la matriz icone
        icone = [icone; n1, n2, n4, n3];
    end
end

##icone


####INCISO4######
##% Configuración inicial
model.nnodes = size(xnode,1);
model.k = ones(1,length(xnode));
model.c = zeros(1,length(xnode));
model.G = zeros(1,length(xnode));
model.rho = 1.0;
model.cp = 1.0;
model.dt = 1/400; #-> despejo de la ecuacion del inciso deltat y da 0.0025
model.maxit = 1000; #Num de iteraciones = t/deltat -> 0.025 / 0.0025 = 10

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 1;
model.tol = 1.000000e-5;
model.PHI_n = sin(xnode(:,1) .* pi) .* sin(xnode(:,2) .* (2*pi));
[PHI,Q,dt,cantIteraciones] = fdm2d(xnode, icone, DIR, NEU, ROB, model);
##PHI: Contendrá los valores de la temperatura en cada nodo después de 𝑡=0.025
mode = 0;
graph = 0;
##fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);
disp('el dt del esquema explicito es: ')
dt
disp('cantidad de iteraciones')
cantIteraciones
disp('El estado estacionario se alcanza en (iteraciones * dt) segundos: ')
cantIteraciones*dt


