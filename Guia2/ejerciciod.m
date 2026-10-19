% ejercicio 2d) transitorio, la fuente G=100 se apaga a los 2 s
% (basado en el 2b, malla de bajo refinamiento)
close all; clear all; more off;

xnode = [
  0.0000000000000000, 0.0000000000000000;
  0.0000000000000000, 0.1666666666666667;
  0.0000000000000000, 0.3333333333333333;
  0.0000000000000000, 0.5000000000000000;
  0.0000000000000000, 0.6666666666666666;
  0.0000000000000000, 0.8333333333333333;
  0.0000000000000000, 1.0000000000000000;
  0.1666666666666667, 0.0000000000000000;
  0.1666666666666667, 0.1666666666666667;
  0.1666666666666667, 0.3333333333333333;
  0.1666666666666667, 0.5000000000000000;
  0.1666666666666667, 0.6666666666666666;
  0.1666666666666667, 0.8333333333333333;
  0.1666666666666667, 1.0000000000000000;
  0.3333333333333333, 0.0000000000000000;
  0.3333333333333333, 0.1666666666666667;
  0.3333333333333333, 0.3333333333333333;
  0.3333333333333333, 0.5000000000000000;
  0.3333333333333333, 0.6666666666666666;
  0.3333333333333333, 0.8333333333333333;
  0.3333333333333333, 1.0000000000000000;
  0.5000000000000000, 0.0000000000000000;
  0.5000000000000000, 0.1666666666666667;
  0.5000000000000000, 0.3333333333333333;
  0.6666666666666666, 0.0000000000000000;
  0.6666666666666666, 0.1666666666666667;
  0.6666666666666666, 0.3333333333333333;
  0.8333333333333333, 0.0000000000000000;
  0.8333333333333333, 0.1666666666666667;
  0.8333333333333333, 0.3333333333333333;
  1.0000000000000000, 0.0000000000000000;
  1.0000000000000000, 0.1666666666666667;
  1.0000000000000000, 0.3333333333333333;
];

% escalo los valores a las dimensiones reales (la malla viene en [0,1]x[0,1])
xnode(:,1) = xnode(:,1)*3;
xnode(:,2) = xnode(:,2)*6;

icone = [
       1,      8,      9,      2;
       2,      9,     10,      3;
       3,     10,     11,      4;
       4,     11,     12,      5;
       5,     12,     13,      6;
       6,     13,     14,      7;
       8,     15,     16,      9;
       9,     16,     17,     10;
      10,     17,     18,     11;
      11,     18,     19,     12;
      12,     19,     20,     13;
      13,     20,     21,     14;
      15,     22,     23,     16;
      16,     23,     24,     17;
      22,     25,     26,     23;
      23,     26,     27,     24;
      25,     28,     29,     26;
      26,     29,     30,     27;
      28,     31,     32,     29;
      29,     32,     33,     30;
];

% Nodo 7 (esquina (0,6)): T=20 (izq.) y T=0 (sup.) -> promedio = 10
DIR = [
       1, 20.0;
       2, 20.0;
       3, 20.0;
       4, 20.0;
       5, 20.0;
       6, 20.0;
       7, 10.0;
      14, 0.0;
      21, 0.0;
];

NEU = [
      18, 0.0, 2;
      19, 0.0, 2;
      20, 0.0, 2;
      21, 0.0, 2;
      24, 0.0, 3;
      27, 0.0, 3;
      30, 0.0, 3;
      31, 0.0, 2;
      32, 0.0, 2;
      33, 0.0, 2;
      33, 0.0, 3;
];

ROB = [
       1, 2.0, 50.0, 1;
       8, 2.0, 50.0, 1;
      15, 2.0, 50.0, 1;
      22, 2.0, 50.0, 1;
      25, 2.0, 50.0, 1;
      28, 2.0, 50.0, 1;
      31, 2.0, 50.0, 1;
];

disp('---------------------------------------------------------------');
disp('Inicializando modelo de datos...');

model.nnodes = size(xnode,1);

model.k = ones(model.nnodes,1);      % k = 1
model.c = zeros(model.nnodes,1);     % sin reaccion

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 1;

% Parámetros para esquemas temporales
model.rho = 1.0;
model.cp  = 1.0;
model.dt  = 0.05;
t_off     = 2.0;                     % instante en que se apaga la fuente

% Condición inicial
model.PHI_n = zeros(model.nnodes,1);

disp('Iniciando el método numérico...');

% ---------------- Fase 1: fuente encendida, 0 <= t <= 2 s ----------------
model.G     = 100*ones(model.nnodes,1);
model.maxit = round(t_off/model.dt);   % 40 pasos
model.tol   = 0;                       % que no corte antes de t = 2 s
[PHI1,Q1] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

% ---------------- Fase 2: fuente apagada, t > 2 s ----------------
model.G     = zeros(model.nnodes,1);
model.PHI_n = PHI1(:,end);             % arranca desde el estado en t = 2 s
model.maxit = 100;                    % tope de pasos
model.tol   = 1.0e-5;                  % corta al llegar al estacionario
[PHI2,Q2] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

% Se unen ambas fases (la 1ra columna de PHI2 / las 2 primeras de Q2
% repiten el ultimo estado de la fase 1)
PHI = [PHI1, PHI2(:,2:end)];
Q   = [Q1,   Q2(:,3:end)];

printf('Columnas fase 1: %d (esperado %d)\n', size(PHI1,2), round(t_off/model.dt)+1);
printf('Pasos fase 2: %d  -> tiempo final = %.2f s\n', size(PHI2,2)-1, ...
       t_off + (size(PHI2,2)-1)*model.dt);

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