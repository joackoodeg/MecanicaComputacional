% Recuperatorio Parcial 1 - 27/11/2025
%   Diferencias Finitas - Ecuacion del calor, estado estacionario, con fuente
%   Dominio rectangular Lx=0.50m x Ly=0.20m, centrado en el origen
%   k = 210 [W/m.C],  Q = 500 [W] (potencia total generada en la placa)
%   Izquierda (x=-Lx/2): Neumann q = 0   [W/m2]  (aislado)
%   Derecha   (x=+Lx/2): Neumann q = 50  [W/m2]  (sale calor)
%   Superior  (y=+Ly/2): Robin h = 10 [W/m2.C], Tinf = 20 [C]
%   Inferior  (y=-Ly/2): Robin h = 1  [W/m2.C], Tinf = 20 [C]
%   No hay condiciones Dirichlet.
%   Punto A = (0,0) (centro del dominio).  Punto B = (0.15,-0.05)
clear all; close all; clc;

%% Datos del modelo
Lx = 0.50; Ly = 0.20;
x1 = -Lx/2; x2 = Lx/2;
y1 = -Ly/2; y2 = Ly/2;

k_val = 210;                 % [W/m.C]
c_val = 0;                    % sin termino reactivo

% El enunciado da la potencia TOTAL generada, Q = 500 [W]. Para usarla en la
% ecuacion -k*Laplaciano(T) = G hace falta la fuente POR UNIDAD DE VOLUMEN
% (espesor unitario -> por unidad de area): G = Q / (Lx*Ly).
Q_tot = 500;
G_val = Q_tot / (Lx*Ly);      % [W/m^3] = 5000

qL = 0;    dirL = 4;           % izquierda: Neumann, saliente en -x (W)
qR = 50;   dirR = 2;           % derecha:   Neumann, saliente en +x (E)
hTop = 10; TinfTop = 20; dirTop = 3;   % superior: Robin, saliente en +y (N)
hBot = 1;  TinfBot = 20; dirBot = 1;   % inferior: Robin, saliente en -y (S)

xA = 0; yA = 0;
xB = 0.15; yB = -0.05;

%% ========================================================================
%% 1) Existencia y unicidad de solucion sin condiciones Dirichlet
%% ========================================================================
disp('1) Existencia/unicidad sin Dirichlet');
disp(['   Si TODO el contorno fuera Neumann, la matriz K (que discretiza ', ...
      '-k*Laplaciano) resulta singular: un vector constante T=cte queda en ', ...
      'su nucleo (sumar una constante a toda la solucion no cambia ningun ', ...
      'flujo). En ese caso solo hay solucion si se cumple una condicion de ', ...
      'compatibilidad (balance global de energia) y, de existir, queda ', ...
      'definida a menos de una constante aditiva.']);
disp(['   Aqui, en cambio, los bordes superior e inferior son de tipo Robin ', ...
      '(convectivos): la condicion liga el flujo a (T - Tinf), es decir, al ', ...
      'VALOR de T y no solo a su derivada. Sumar una constante c a toda la ', ...
      'solucion ya NO deja invariante esa ecuacion (h*(T+c-Tinf) =/= h*(T-Tinf) ', ...
      'si c =/= 0), por lo que el vector constante deja de pertenecer al nucleo ', ...
      'de K. En el ensamble, fdm2d_robin suma 2*h/d > 0 en la diagonal de cada ', ...
      'nodo Robin, lo que vuelve a K no singular (diagonal estrictamente dominante ', ...
      'en esos nodos).']);
disp(['   Conclusion: SI es posible llegar a una solucion unica, sin necesidad ', ...
      'de ninguna condicion Dirichlet ni de verificar un balance global de ', ...
      'energia previo, porque la condicion Robin fija el nivel absoluto de T.']);

%% ========================================================================
%% 2-6) Discretizacion mas simple con A como nodo: malla 3x3
%% ========================================================================
% La forma mas simple de tener a A=(0,0) como nodo INTERIOR (con los 4
% vecinos reales: uno en cada borde) es una malla de 3x3 nodos, 2x2
% elementos: dx = Lx/2 = 0.25, dy = Ly/2 = 0.10. Asi A queda exactamente en
% el centro, rodeado por un nodo de cada borde (izq., der., sup., inf.).
nx = 3; ny = 3;
xcoord = linspace(x1,x2,nx);
ycoord = linspace(y1,y2,ny);
dx = xcoord(2) - xcoord(1);     % 0.25
dy = ycoord(2) - ycoord(1);     % 0.10

node = @(i,j) (j-1)*nx + i;

nnodes = nx*ny;
xnode = zeros(nnodes,2);
for j = 1:ny
    for i = 1:nx
        xnode(node(i,j),:) = [xcoord(i), ycoord(j)];
    end
end

nelem = (nx-1)*(ny-1);
icone = zeros(nelem,4);
e = 1;
for j = 1:ny-1
    for i = 1:nx-1
        icone(e,:) = [node(i,j) node(i+1,j) node(i+1,j+1) node(i,j+1)];
        e = e+1;
    end
end

%% Condiciones de borde (sin Dirichlet; en las esquinas se suman las dos)
DIR = [];
NEU = []; ROB = [];
for j = 1:ny
    NEU = [NEU; node(1,j),  qL, dirL];   % izquierda
    NEU = [NEU; node(nx,j), qR, dirR];   % derecha
end
for i = 1:nx
    ROB = [ROB; node(i,1),  hBot, TinfBot, dirBot];  % inferior
    ROB = [ROB; node(i,ny), hTop, TinfTop, dirTop];  % superior
end

A = node(2,2);    % nodo central = punto A

printf('\n2) Malla propuesta: %dx%d nodos (dx=%.3f m, dy=%.3f m), A = nodo %d\n', nx, ny, dx, dy, A);
printf('3) Grados de libertad: no hay nodos Dirichlet -> los %d nodos son incognitas (DOF = %d)\n', nnodes, nnodes);

%% Ensamble manual (para poder mostrar K y F) y resolucion
[K,F] = fdm2d_initialize(nnodes);
[neighb] = fdm2d_neighbors(icone);
[K,F] = fdm2d_gen_system(K,F,xnode,neighb, k_val*ones(nnodes,1), c_val*ones(nnodes,1), G_val*ones(nnodes,1));
[F]   = fdm2d_neumann(F,xnode,neighb,NEU);
[K,F] = fdm2d_robin(K,F,xnode,neighb,ROB);
% (no se aplica fdm2d_dirichlet: no hay nodos Dirichlet)

PHI = K \ F;

printf('\n4) Temperaturas en los nodos de la malla 3x3 [C] (Tgrid(i,j), i=x, j=y):\n');
Tgrid = reshape(full(PHI), nx, ny);
disp(Tgrid);
printf('   T(A) = T(0,0) = %.4f C\n', full(PHI(A)));

printf('\n5) Matriz global K y vector F del sistema (9x9):\n');
disp(full(K));
disp(full(F));

%% 6) Stencil del nodo A, paso a paso
% Balance en diferencias finitas centradas para un nodo INTERIOR (A tiene
% los 4 vecinos reales, a distancia dx en x y dy en y a cada lado):
%
%   -k*(T_E - 2*T_A + T_W)/dx^2 - k*(T_N - 2*T_A + T_S)/dy^2 = G_A
%
% Reordenando (forma K*T = F que entrega fdm2d_gen_system, con ax=cx=1/dx^2,
% bx=-2/dx^2, ay=cy=1/dy^2, by=-2/dy^2 al ser dx, dy simetricos a ambos lados):
%
%   (2k/dx^2 + 2k/dy^2)*T_A - (k/dx^2)*(T_E+T_W) - (k/dy^2)*(T_N+T_S) = G_A
%
E = node(3,2); W = node(1,2); N = node(2,3); S = node(2,1);
kdx2 = k_val/dx^2; kdy2 = k_val/dy^2;
printf('\n6) Stencil del nodo A (forma final simplificada, dx=%.2f, dy=%.2f):\n', dx, dy);
printf('   (2k/dx^2 + 2k/dy^2)*T_A - (k/dx^2)*(T_E+T_W) - (k/dy^2)*(T_N+T_S) = G_A\n');
printf('   k/dx^2 = %.1f ,  k/dy^2 = %.1f\n', kdx2, kdy2);
printf('   %.1f*T_A - %.1f*(T_E+T_W) - %.1f*(T_N+T_S) = %.1f\n', 2*kdx2+2*kdy2, kdx2, kdy2, G_val);
printf('   (coincide con la fila de K del nodo A: K(A,A)=%.1f, K(A,E)=K(A,W)=%.1f, K(A,N)=K(A,S)=%.1f, F(A)=%.1f)\n', ...
        full(K(A,A)), full(K(A,E)), full(K(A,N)), full(F(A)));

%% ========================================================================
%% 7) Temperatura y flujo en B=(0.15,-0.05), que NO es nodo de esta malla
%% ========================================================================
% B cae dentro del elemento de esquinas (0,-0.1)-(0.25,-0.1)-(0.25,0)-(0,0).
P00 = node(2,1); P10 = node(3,1); P01 = node(2,2); P11 = node(3,2);
T00 = full(PHI(P00)); T10 = full(PHI(P10)); T01 = full(PHI(P01)); T11 = full(PHI(P11));

xi  = (xB - xnode(P00,1)) / dx;      % 0.6
eta = (yB - xnode(P00,2)) / dy;      % 0.5

TB_coarse = T00*(1-xi)*(1-eta) + T10*xi*(1-eta) + T01*(1-xi)*eta + T11*xi*eta;

% Gradiente de la misma interpolacion bilineal (lineal por franjas) evaluado en B
dTdx_B = ((T10-T00)*(1-eta) + (T11-T01)*eta) / dx;
dTdy_B = ((T01-T00)*(1-xi)  + (T11-T10)*xi ) / dy;
qB_coarse = -k_val*[dTdx_B, dTdy_B];

printf('\n7) B=(0.15,-0.05) no es nodo -> interpolacion bilineal dentro del elemento que lo contiene:\n');
printf('   T(B) = %.4f C\n', TB_coarse);
printf('   q(B) = [%.4f, %.4f] W/m^2 ,  |q(B)| = %.4f W/m^2\n', qB_coarse(1), qB_coarse(2), norm(qB_coarse));

%% ========================================================================
%% 8) Refinar la malla para que B sea un nodo, y comparar
%% ========================================================================
% Con paso uniforme h=0.05 tanto A=(0,0) como B=(0.15,-0.05) caen justo
% sobre nodos de la malla (h divide exactamente a 0.25, 0.10, 0.15 y 0.05).
h2 = 0.05;
nx2 = round(Lx/h2) + 1;   % 11
ny2 = round(Ly/h2) + 1;   % 5

xcoord2 = linspace(x1,x2,nx2);
ycoord2 = linspace(y1,y2,ny2);
dx2 = xcoord2(2)-xcoord2(1);
dy2 = ycoord2(2)-ycoord2(1);

node2 = @(i,j) (j-1)*nx2 + i;

nnodes2 = nx2*ny2;
xnode2 = zeros(nnodes2,2);
for j = 1:ny2
    for i = 1:nx2
        xnode2(node2(i,j),:) = [xcoord2(i), ycoord2(j)];
    end
end

nelem2 = (nx2-1)*(ny2-1);
icone2 = zeros(nelem2,4);
e = 1;
for j = 1:ny2-1
    for i = 1:nx2-1
        icone2(e,:) = [node2(i,j) node2(i+1,j) node2(i+1,j+1) node2(i,j+1)];
        e = e+1;
    end
end

NEU2 = []; ROB2 = [];
for j = 1:ny2
    NEU2 = [NEU2; node2(1,j),   qL, dirL];
    NEU2 = [NEU2; node2(nx2,j), qR, dirR];
end
for i = 1:nx2
    ROB2 = [ROB2; node2(i,1),   hBot, TinfBot, dirBot];
    ROB2 = [ROB2; node2(i,ny2), hTop, TinfTop, dirTop];
end

[K2,F2] = fdm2d_initialize(nnodes2);
[neighb2] = fdm2d_neighbors(icone2);
[K2,F2] = fdm2d_gen_system(K2,F2,xnode2,neighb2, k_val*ones(nnodes2,1), c_val*ones(nnodes2,1), G_val*ones(nnodes2,1));
[F2]    = fdm2d_neumann(F2,xnode2,neighb2,NEU2);
[K2,F2] = fdm2d_robin(K2,F2,xnode2,neighb2,ROB2);

PHI2 = K2 \ F2;
Q2   = fdm2d_flux(PHI2, neighb2, xnode2, k_val*ones(nnodes2,1));

iB = round((xB-x1)/dx2) + 1;   % 9
jB = round((yB-y1)/dy2) + 1;   % 2
PB = node2(iB,jB);

TB_fine = full(PHI2(PB));
qB_fine = full(Q2(PB,:));

iA2 = round((xA-x1)/dx2) + 1;
jA2 = round((yA-y1)/dy2) + 1;
A2  = node2(iA2,jA2);

printf('\n8) Malla refinada %dx%d (dx=dy=%.2f m) -> B es ahora el nodo %d\n', nx2, ny2, dx2, PB);
printf('   T(A) malla fina = %.4f C   (malla gruesa: %.4f C)\n', full(PHI2(A2)), full(PHI(A)));
printf('   T(B) malla fina = %.4f C   (interpolado en malla gruesa: %.4f C)\n', TB_fine, TB_coarse);
printf('   q(B) malla fina = [%.4f, %.4f] W/m^2, |q(B)| = %.4f W/m^2\n', qB_fine(1), qB_fine(2), norm(qB_fine));
printf('   q(B) interpolado malla gruesa = [%.4f, %.4f] W/m^2, |q(B)| = %.4f W/m^2\n', qB_coarse(1), qB_coarse(2), norm(qB_coarse));
printf('   Diferencia relativa T(B): %.2f %%\n', 100*abs(TB_fine-TB_coarse)/abs(TB_fine));
printf('   Diferencia relativa |q(B)|: %.2f %%\n', 100*abs(norm(qB_fine)-norm(qB_coarse))/norm(qB_fine));

%% Graficos (malla refinada)
fdm2d_graph_mesh(full(PHI2), Q2, xnode2, icone2, 0, 0);   % Temperatura
fdm2d_graph_mesh(full(PHI2), Q2, xnode2, icone2, 0, 1);   % Flujo de calor
