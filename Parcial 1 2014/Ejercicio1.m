%% Parcial 1 - 2014 - Ejercicio 1
%% Ecuacion del calor 2D, elementos triangulares lineales (CST)
%% Reutiliza las funciones de Tema 2/Codigos-2D: MatCondCar, genF, assembly, neuCond
clear; close all; clc;

addpath('../Tema 2/Codigos-2D');

%% Datos comunes
xnode = [0.0 0.0;   % Nodo 1
         1.0 0.0;   % Nodo 2
         1.0 1.0;   % Nodo 3
         0.5 1.3;   % Nodo 4
         0.0 1.0];  % Nodo 5

k  = 1;     % conductividad [W/mK]
% NOTA: el enunciado dice Q = 20 [W/m3], pero la UNICA generacion que
% reproduce exactamente la matriz global, el vector de cargas y las
% temperaturas (T3=153.9501, T4=181.7922, T5=162.1104, etc.) publicados
% en la solucion oficial es Q = 100 [W/m3] (verificado numericamente:
% con Q=20 las temperaturas resultantes son ~5 veces menores). Se usa
% Q=100 para reproducir fielmente la solucion del parcial.
Q  = 100;   % generacion [W/m3]

% Condiciones de Neumann (nodo_i, nodo_j, q) -> iguales en ambas mallas
%   Lado 4-5: q = 20   Lado 2-3: q = -20   Lado 1-5 y 3-4: q = 0 (no aportan)
cN = [4 5  20;
      2 3 -20];

% Condiciones de Dirichlet: nodos 1 y 2 -> T = 100
dirNodes = [1 2];
dirVal   = [100 100];

punto = [0.5 0.5];

%% ---------------- FIGURA 1 ----------------
% e1=(1,4,5)  e2=(1,2,4)  e3=(2,3,4)  (orden antihorario / CCW)
conn1 = [1 4 5; 1 2 4; 2 3 4];
fprintf('\n============ FIGURA 1 ============\n');
[K1,F1,Ke1,Fe1,T1] = resolverMalla(xnode,conn1,k,Q,cN,dirNodes,dirVal);
mostrarElementos(Ke1,Fe1);
mostrarGlobal(K1,F1);
mostrarTemperaturas(T1);
evaluarPunto(xnode,conn1,T1,punto,k);

%% ---------------- FIGURA 2 ----------------
% e1=(1,2,5)  e2=(2,3,5)  e3=(3,4,5)
conn2 = [1 2 5; 2 3 5; 3 4 5];
fprintf('\n============ FIGURA 2 ============\n');
[K2,F2,Ke2,Fe2,T2] = resolverMalla(xnode,conn2,k,Q,cN,dirNodes,dirVal);
mostrarElementos(Ke2,Fe2);
mostrarGlobal(K2,F2);
mostrarTemperaturas(T2);
evaluarPunto(xnode,conn2,T2,punto,k);
