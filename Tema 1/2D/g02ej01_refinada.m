close all; clear all; more off;

% Refinar la malla
dx = 0.00125; % La mitad de dx original
dy = -0.002;  % La mitad de dy original

xnode_refinada = [
  0, 0; # Nodo 1
  dx, 0; # Nodo 2
  dx*2, 0; # Nodo 3
  dx*3, 0; # Nodo 4
  dx*4, 0; # Nodo 5
  dx*5, 0; # Nodo 6
  dx*6, 0; # Nodo 7
  dx*7, 0; # Nodo 8
  dx*8, 0; # Nodo 9

  0, dy; # Nodo 10
  dx, dy; # Nodo 11
  dx*2, dy; # Nodo 12
  dx*3, dy; # Nodo 13
  dx*4, dy; # Nodo 14
  dx*5, dy; # Nodo 15
  dx*6, dy; # Nodo 16
  dx*7, dy; # Nodo 17
  dx*8, dy; # Nodo 18

  0, dy*2; #Nodo 19
  dx, dy*2; #Nodo 20
  dx*2, dy*2; #Nodo 21
  dx*3, dy*2; #Nodo 22
  dx*4, dy*2; #Nodo 23
  dx*5, dy*2; #Nodo 24
  dx*6, dy*2; #Nodo 25
  dx*7, dy*2; #Nodo 26
  dx*8, dy*2; #Nodo 27

  dx*4, dy*3; #Nodo 28
  dx*5, dy*3; #Nodo 29
  dx*6, dy*3; #Nodo 30
  dx*7, dy*3; #Nodo 31
  dx*8, dy*3; #Nodo 32

  dx*4, dy*4; #Nodo 33
  dx*5, dy*4; #Nodo 34
  dx*6, dy*4; #Nodo 35
  dx*7, dy*4; #Nodo 36
  dx*8, dy*4; #Nodo 37

  dx*4, dy*5; #Nodo 38
  dx*5, dy*5; #Nodo 39
  dx*6, dy*5; #Nodo 40
  dx*7, dy*5; #Nodo 41
  dx*8, dy*5; #Nodo 42

   0, dy*6;  #Nodo 43
  dx, dy*6; #Nodo 44
  dx*2, dy*6; #Nodo 45
  dx*3, dy*6; #Nodo 46
  dx*4, dy*6; #Nodo 47
  dx*5, dy*6; #Nodo 48
  dx*6, dy*6; #Nodo 49
  dx*7, dy*6; #Nodo 50
  dx*8, dy*6; #Nodo 51

  0, dy*7;  #Nodo 52
  dx, dy*7; #Nodo 53
  dx*2, dy*7; #Nodo 54
  dx*3, dy*7; #Nodo 55
  dx*4, dy*7; #Nodo 56
  dx*5, dy*7; #Nodo 57
  dx*6, dy*7; #Nodo 58
  dx*7, dy*7; #Nodo 59
  dx*8, dy*7; #Nodo 60

  0, dy*8; #Nodo 61
  dx, dy*8; #Nodo 62
  dx*2, dy*8; #Nodo 63
  dx*3, dy*8; #Nodo 64
  dx*4, dy*8; #Nodo 65
  dx*5, dy*8;#Nodo 66
  dx*6, dy*8;#Nodo 67
  dx*7, dy*8;#Nodo 68
  dx*8, dy*8;#Nodo 69

  0, dy*9; #Nodo 70
  dx, dy*9; #Nodo 71
  dx*2, dy*9; #Nodo 72
  dx*3, dy*9; #Nodo 73
  dx*4, dy*9; #Nodo 74
  dx*5, dy*9; #Nodo 75
  dx*6, dy*9; #Nodo 76
  dx*7, dy*9; #Nodo 77
  dx*8, dy*9; #Nodo 78

  0, dy*10; #Nodo 79
  dx, dy*10; #Nodo 80
  dx*2, dy*10;#Nodo 81
  dx*3, dy*10;#Nodo 82
  dx*4, dy*10;#Nodo 83
  dx*5, dy*10;#Nodo 84
  dx*6, dy*10;#Nodo 85
  dx*7, dy*10;#Nodo 86
  dx*8, dy*10;#Nodo 87

];

% Plot de los nuevos nodos
plot(xnode_refinada(:,1), xnode_refinada(:,2), 'o');

icone_refinada = [

79,80,71,70;
80,81,72,71;
81,82,73,72;
82,83,74,73;
83,84,75,74;
84,85,76,75;
85,86,77,76;
86,87,78,77;

70,71,62,61;
71,72,63,62;
72,73,64,63;
73,74,65,64;
74,75,66,65;
75,76,67,66;
76,77,68,67;
77,78,69,68;

61,62,53,52;
62,63,54,53;
63,64,55,54;
64,65,56,55;
65,66,57,56;
66,67,58,57;
67,68,59,58;
68,69,60,59;

52,53,44,43;
53,54,45,44;
54,55,46,45;
55,56,47,46;
56,57,48,47;
57,58,49,48;
58,59,50,49;
59,60,51,50;

47, 48, 39, 38;
48, 49, 40, 39;
49 ,50, 41, 40;
50, 51, 42, 41;

38,39,34,33;
39,40,35,34;
40,41,36,35;
41,42,37,36;

33,34,29,28;
34,35,30,29;
35,36,31,30;
36,37,32,31;

28,29,24,23;
29,30,25,24;
30,31,26,25;
31,32,27,26;

19,20,11,10;
20,21,12,11;
21,22,13,12;
22,23,14,13;
23,24,15,14;
24,25,16,15;
25,26,17,16;
26,27,18,17;

10, 11, 2, 1;
11, 12, 3, 2;
12, 13, 4, 3;
13, 14, 5, 4;
14, 15, 6, 5;
15, 16, 7, 6;
16, 17, 8, 7;
17, 18, 9, 8;
];

NEU = [

];

ROB = [
];

DIR = [];



% Configuración del modelo con malla refinada
model.nnodes = size(xnode_refinada, 1);
model.k = 190 * ones(1, length(xnode_refinada));
model.c = zeros(1, length(xnode_refinada));
model.G = zeros(1, length(xnode_refinada));
model.ts = 2;
model.rho = 0.0;
model.cp = 0;
model.maxit = 1;
model.tol = 1.000000e-05;
model.PHI_n = zeros(model.nnodes, 1);

% Llamada principal al Método de Diferencias Finitas
disp('Iniciando el método numérico con malla refinada');
[PHI_refinada, Q_refinada] = fdm2d(xnode_refinada, icone_refinada, DIR, NEU, ROB, model);

% Graficar resultados
mode = 0;
graph = 0;
fdm2d_graph_mesh(full(PHI_refinada), Q_refinada, xnode_refinada, icone_refinada, mode, graph);

% Mostrar resultados de flujo refinado en los mismos nodos
disp('El flujo de calor en nodo 9:');
qx = Q_refinada(9, 1);
qy = Q_refinada(9, 2);
q9 = sqrt(qx^2 + qy^2);
disp(q9);

disp('El flujo de calor en nodo 17:');
qx = Q_refinada(17, 1);
qy = Q_refinada(17, 2);
q17 = sqrt(qx^2 + qy^2);
disp(q17);

disp('El flujo de calor en nodo 20:');
qx = Q_refinada(20, 1);
qy = Q_refinada(20, 2);
q20 = sqrt(qx^2 + qy^2);
disp(q20);

