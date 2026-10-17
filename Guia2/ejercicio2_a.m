close all; clear all; more off;

xnode = [
  0, 0;
  0, 2;
  0, 4;
  0, 6;
  1, 0;
  1, 2;
  1, 4;
  1, 6;
  2, 0;
  2, 2;
  3, 0;
  3, 2;
];

icone = [
       1,      5,      6,      2;
       2,      6,      7,      3;
       3,      7,      8,      4;
       5,      9,     10,      6;
       9,     11,     12,     10;
];

DIR = [
       1, 20.0000000000000000;
       2, 20.0000000000000000;
       3, 20.0000000000000000;
       4, 20.0000000000000000;
       8,  0.0000000000000000;
];

NEU = [
       7, 0.0000000000000000, 2;
       8, 0.0000000000000000, 2;
      10, 0.0000000000000000, 3;
      11, 0.0000000000000000, 2;
      12, 0.0000000000000000, 2;
      12, 0.0000000000000000, 3;
];

ROB = [
       1, 2.0000000000000000, 50.0000000000000000, 1;
       5, 2.0000000000000000, 50.0000000000000000, 1;
       9, 2.0000000000000000, 50.0000000000000000, 1;
      11, 2.0000000000000000, 50.0000000000000000, 1;
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);

model.k = [
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
    1.0000000000000000;
];

model.c = [
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
    0.0000000000000000;
];

model.G = [
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
    100.0000000000000000;
];

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 2;

% Parámetros para esquemas temporales
model.rho = 1.0000000000000000;
model.cp = 1.0000000000000000;
model.maxit =            1;
model.tol = 1.000000e-05;

% Condición inicial
model.PHI_n = mean(DIR(:,2))*ones(model.nnodes,1);

disp('Iniciando el método numérico...');

% Llamada principal al Método de Diferencias Finitas
[PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

disp('Finalizada la ejecución del método numérico.');

disp('---------------------------------------------------------------');
disp('Iniciando el post-procesamiento...');

% mode ---> modo de visualización:
%           [0] 2D - Con malla
%           [1] 3D - Con malla
%           [2] 2D - Sin malla
%           [3] 3D - Sin malla
% graph --> tipo de gráfica:
%           [0] Temperatura (escalar)
%           [1] Flujo de Calor (vectorial)
%           [2] Flujo de Calor eje-x (escalar)
%           [3] Flujo de Calor eje-y (escalar)
%           [4] Magnitud de Flujo de Calor (escalar)
mode = 0;
graph = 0;
fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);

disp('Finalizado el post-procesamiento.');

