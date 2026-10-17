close all; clear all; more off;

dx = 0.0025;
dy = -0.004;
xnode = [
  0, 0; #Nodo 1
  dx, 0; #Nodo 2
  dx*2, 0; #Nodo 3
  dx*3, 0; #Nodo 4
  dx*4 , 0; #Nodo 5

  0, dy; #Nodo 6
  dx , dy; #Nodo 7
  dx*2,dy; #Nodo 8
  dx*3,dy; #Nodo 9
  dx*4,dy; #Nodo 10

  dx*2,dy*2; #Nodo 11
  dx*3,dy*2; #Nodo 12
  dx*4,dy*2; #Nodo 13

  0, dy*3; #Nodo 14
  dx , dy*3; #Nodo 15
  dx*2,dy*3; #Nodo 16
  dx*3,dy*3; #Nodo 17
  dx*4,dy*3; #Nodo 18

  0, dy*4; #Nodo 19
  dx , dy*4; #Nodo 20
  dx*2,dy*4; #Nodo 21
  dx*3,dy*4; #Nodo 22
  dx*4,dy*4; #Nodo 23

  0, dy*5; #Nodo 24
  dx , dy*5; #Nodo 25
  dx*2,dy*5; #Nodo 26
  dx*3,dy*5; #Nodo 27
  dx*4,dy*5 #Nodo 28
];



plot(xnode(:,1),xnode(:,2),'o');

%     (1) Flujo en dirección eje-y, sentido negativo (S – South - Sur)
%     (2) Flujo en dirección eje-x, sentido positivo (E – East - Este)
%     (3) Flujo en dirección eje-y, sentido positivo (N – North - Norte)
%     (4) Flujo en dirección eje-x, sentido negativo (W – West – Oeste)


NEU = [
  1, 0, 3;
  1, 0, 4;
  2, 0, 3;
  3, 0, 3;
  4, 0, 3;
  5, 0, 3;
  5, 0, 2;
  6, 0, 4;
  10, 0, 2;
  13, 0, 2;
  14, 0, 4;
  18, 0, 2;
  19, 0, 4;
  23, 0, 2;
  24, 0, 4;
  24, -10^5, 1;
  25, -10^5, 1;
  26, -10^5, 1;
  27, -10^5, 1;
  28, 0, 2;
  28, -10^5, 1;
];


ROB = [
    6, 5000, 15, 1;
    7, 5000, 15, 1;
    11, 5000, 15, 4; #2 o 4 el flujo ? con 2 no anda, pq 4?
    15, 5000, 15, 3;
    14, 5000, 15, 3;
  ];

DIR = [];

icone = [
  24,25,20,19;
  25,26,21,20;
  26,27,22,21;
  27,28,23,22;

  19,20,15,14;
  20,21,16,15;
  21,22,17,16;
  22,23,18,17;

  16,17,12,11;
  17,18,13,12;

  11,12,9,8;
  12,13,10,9;

  6,7,2,1;
  7,8,3,2;
  8,9,4,3;
  9,10,5,4
];

##RESOLUCION##

model.nnodes = size(xnode,1);
model.k = 190 * ones(1,length(xnode));
model.c = zeros(1,length(xnode));
model.G = zeros(1,length(xnode));

% Esquema Temporal: [0] Explícito, [1] Implícito, [X] Estacionario
model.ts = 2;

% Parámetros para esquemas temporales
model.rho = 0.0000000000000000;
model.cp = 0.0000000000000000;
model.maxit =            1;
model.tol = 1.000000e-05;

% Condición inicial
model.PHI_n = zeros(model.nnodes,1);
disp('Iniciando el método numérico...');

% Llamada principal al Método de Diferencias Finitas
[PHI,Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);

mode = 0;
graph = 0;
fdm2d_graph_mesh(full(PHI),Q,xnode,icone,mode,graph);;

disp('El flujo de calor en nodo 9:')
qx=Q(9,1)
qy=Q(9,2)
q9=sqrt(qx^2+qy^2);
disp(q9)

disp('El flujo de calor en nodo 17')
qx=Q(17,1)
qy=Q(17,2)
q17=sqrt(qx^2+qy^2);
disp(q17)

disp('El flujo de calor en nodo 20')
qx=Q(20,1)
qy=Q(20,2)
q20=sqrt(qx^2+qy^2);
disp(q20)
