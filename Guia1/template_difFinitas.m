% =====================================================================
% TEMPLATE - Diferencias finitas 1D con difFinitas.m
%   rho*cp * dT/dt = k * d2T/dx2 - c*T + G(x)
% Correr desde la carpeta que tiene difFinitas.m y dirchlet.m
% (o descomentar el addpath).
% =====================================================================
clear; clc; close all;
% addpath('ruta/a/Tema 1/1D');

%% ---------------- 1. DATOS DEL MODELO --------------------------------
% Pasar la ecuacion a la forma  rho*cp*T_t = k*T'' - c*T + G
%   ej: T_t - T'' = 100          -> rho*cp=1, k=1, c=0, G=100
%   ej: k*T'' + G = 0            -> estacionario, solo k y G
%   ej: T'' - T = 0              -> k=1, c=1, G=0
m.rho = 1;          % densidad        (rho*cp solo importa en transitorio)
m.cp  = 1;          % calor especifico
m.k   = 1;          % conductividad
m.c   = 0;          % reaccion (+c*T resta). Ojo: si la ec. tiene +T, c es negativo
Gfun  = @(x) 100*ones(size(x));    % fuente; ej: @(x) 10*x.^3, @(x) 0*x

%% ---------------- 2. DOMINIO Y MALLA ---------------------------------
L1 = 0;  L2 = 1;    % extremos
Nnod = 41;          % cantidad TOTAL de nodos (ojo: "n elementos" = n+1 nodos)
x  = linspace(L1, L2, Nnod);
dx = x(2) - x(1);
m.G = Gfun(x);      % difFinitas pide G evaluada en cada nodo (vector fila)

%% ---------------- 3. CONDICIONES DE BORDE ----------------------------
% Una fila por borde: [tipo  valor  extra]
%   Dirichlet: [1  T0   -1]      T = T0
%   Neumann  : [2  q    -1]      -k*dT/dn = q   (q > 0 = calor que SALE)
%   Robin    : [3  h  Tinf]      -k*dT/dn = h*(T - Tinf)
% Traducir el enunciado a esas formas (n = normal exterior):
%   derecha:   dT/dn =  dT/dx      izquierda: dT/dn = -dT/dx
%   ej: T'(1) = -T(1)    -> Robin derecha  h=1, Tinf=0   (con k=1)
%   ej: T'(L) = 0        -> Neumann q=0 (aislado)
%   ej: T'(0) = 5 (k=1)  -> izquierda: -k*dT/dn = k*T'(0) = 5 -> q=5
cb = [1  0  -1;     % borde izquierdo (x = L1)
      3  1   0];    % borde derecho   (x = L2)

%% ---------------- 4. CONDICION INICIAL (solo transitorio) ------------
Ti = 10*x.*(1-x).^2;          % T(x,0), vector fila
% Ti = zeros(1,Nnod);         % si arranca de cero
% Conviene que Ti respete las Dirichlet (ej. Ti(1) = T0)

%% ---------------- 5. ESQUEMA TEMPORAL --------------------------------
% et = [esquema  maxPasos  tol  dt]
%   esquema: 0 estacionario | 1 explicito (Forward Euler) | 2 implicito (Backward Euler)
%   maxPasos: maximo de pasos de tiempo
%   tol: corte cuando norm(T^{n+1}-T^n)/norm(T^{n+1}) < tol
%   dt: paso de tiempo (el explicito lo ignora y usa 0.9*dt critico)
% Si te dan lambda = k*dt/dx^2  ->  dt = lambda*dx^2/k   (con rho*cp = 1;
%   en general dt = lambda*rho*cp*dx^2/k)
lambda = 40;
dt  = lambda*(m.rho*m.cp)*dx^2/m.k;
esquema = 2;                  % 0 / 1 / 2
et  = [esquema  1e5  1e-4  dt];

%% ---------------- 6. RESOLVER ----------------------------------------
if esquema == 0
  T = difFinitas(x, m, cb, [0 0 0 0], zeros(1,Nnod));   % vector columna
  Tf = T;
else
  T = difFinitas(x, m, cb, et, Ti);    % matriz: filas = nodos, columnas = tiempos
  Tf = T(:,end);                       % ultimo instante
  if esquema == 1
    dtc = 0.5*(m.rho*m.cp)*dx^2/m.k;   % dt critico explicito (sin reaccion)
    printf('Explicito: dt usado = 0.9*dt_crit = %g\n', 0.9*dtc);
    dt = 0.9*dtc;                      % el codigo usa este dt internamente
  end
  printf('Pasos: %d, t final = %g\n', columns(T)-1, (columns(T)-1)*dt);
end

%% ---------------- 7. POSPROCESO (descomentar lo que pida) ------------
% --- Temperatura en un punto (debe ser nodo de la malla) ---
xp = 0.5;
ip = round((xp - L1)/dx) + 1;          % indice del nodo
printf('T(%g) final = %.4f\n', xp, Tf(ip));

% --- Varios puntos ---
xs = [0.25 0.5 0.75 1];
is = round((xs - L1)/dx) + 1;
disp([xs' Tf(is)]);

% --- Primer instante en que T(xp) supera un valor (transitorio) ---
if esquema ~= 0
  umbral = 23;
  j = find(T(ip,:) > umbral, 1);       % columna j <-> t = (j-1)*dt
  if isempty(j)
    disp('Nunca supera el umbral');
  else
    printf('Primer t con T(%g) > %g: t = %g (paso %d)\n', xp, umbral, (j-1)*dt, j-1);
  end
end

% --- Error contra la solucion exacta (si la dan) ---
Tex = @(x) 25*x.*(3-2*x);
Te  = Tex(x(:));
err_rel = norm(Tf - Te)/norm(Te);      % relativo en norma 2
err_max = max(abs(Tf - Te));           % maximo absoluto
printf('Error relativo = %e, error maximo = %e\n', err_rel, err_max);

% --- Grafico final vs exacta ---
figure; plot(x, Tf, 'b-o', x, Te, 'r--'); grid on;
legend('MDF', 'Exacta'); xlabel('x'); ylabel('T');
