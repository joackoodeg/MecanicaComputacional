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

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 1;
model.dt = 0.05;

model.rho = 1.0000000000000000;
model.cp  = 1.0000000000000000;

% Vecindad de la malla (no cambia entre etapas, se calcula una sola vez)
[neighb] = fdm2d_neighbors(icone);

% ------------------------------------------------------------------
% Etapa 1: 0 <= t < 2s, fuente G=100 uniforme (igual que el item b)
% ------------------------------------------------------------------
disp('Etapa 1: 0 a 2s, con fuente G=100...');

t_apagado = 2;                        % instante en que se apaga la fuente
n1 = round(t_apagado/model.dt);       % cantidad de pasos hasta t=2s (40)

model.G = 100*ones(model.nnodes,1);
model.PHI_n = zeros(model.nnodes,1);  % condicion inicial T(x,0)=0
model.maxit = n1;
model.tol = 0;                        % para que no corte antes de los 40 pasos

[K,F] = fdm2d_initialize(model.nnodes);
[K,F] = fdm2d_gen_system(K,F,xnode,neighb,model.k,model.c,model.G);
[F] = fdm2d_neumann(F,xnode,neighb,NEU);
[K,F] = fdm2d_robin(K,F,xnode,neighb,ROB);
[K,F] = fdm2d_dirichlet(K,F,DIR);

[PHI1,Q1] = fdm2d_implicit(K,F,xnode,neighb,model,model.dt);

% ------------------------------------------------------------------
% Etapa 2: t >= 2s, se apaga la fuente (G=0), hasta estado estacionario
% ------------------------------------------------------------------
disp('Etapa 2: desde t=2s, sin fuente, hasta estado estacionario...');

model.G = zeros(model.nnodes,1);
model.PHI_n = PHI1(:,end);            % continua desde el ultimo estado de la etapa 1
model.maxit = 20000;                  % cota de seguridad, corta por tolerancia
model.tol = 1.000000e-05;

[K,F] = fdm2d_initialize(model.nnodes);
[K,F] = fdm2d_gen_system(K,F,xnode,neighb,model.k,model.c,model.G);
[F] = fdm2d_neumann(F,xnode,neighb,NEU);
[K,F] = fdm2d_robin(K,F,xnode,neighb,ROB);
[K,F] = fdm2d_dirichlet(K,F,DIR);

[PHI2,Q2] = fdm2d_implicit(K,F,xnode,neighb,model,model.dt);

disp('Finalizada la ejecución del método numérico.');

% ------------------------------------------------------------------
% Union de ambas etapas (se descarta la primera columna de la etapa 2,
% que duplica el ultimo estado de la etapa 1)
% ------------------------------------------------------------------
PHI = [PHI1, PHI2(:,2:end)];
Q   = [Q1,   Q2(:,3:end)];

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
