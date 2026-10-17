%% Parcial 1 - 2014 - Ejercicio 2
%% Postprocesamiento de un problema de tension plana ya resuelto por FEM
%% (funciones de forma lineales/bilineales), siguiendo el estilo de
%% Tema 2/Codigos-2D (shape functions fisicas, B = grad(N)).
clear; close all; clc;

%% Datos
E  = 2.1e9;   % [Pa]
nu = 0.3;
t  = 0.005;   % [m] (no se usa para tensiones/deformaciones puntuales)

% D de tension plana
D = E/(1-nu^2) * [1 nu 0; nu 1 0; 0 0 (1-nu)/2];

% Malla (Figura 3)
nodes = [0.0 0.0;   % 1
         1.0 0.0;   % 2
         2.0 0.0;   % 3
         2.0 0.5;   % 4
         1.25 0.5;  % 5
         0.0 0.5;   % 6
         0.0 1.0;   % 7
         1.0 1.0;   % 8
         2.0 1.0];  % 9

% Desplazamientos nodales dados [mm]
u = [0.0; 0.7189; -0.05839; 0.55891; 0.40925; 0.0; 0.0; -0.13365; 0.17009];
v = [0.0; 0.66431; 0.0; 1.7722; 0.89222; 0.0; 0.0; 0.45779; 3.9252];

fprintf('D (plano tension):\n'); disp(D);

%% ---------- Punto (0.75, 0.75) -> elemento triangular (6,5,8) ----------
tri = [6 5 8];
xe = nodes(tri,:);
[Nfun,b,c,A] = triShapeCST(xe);
fprintf('\nElemento triangular %d-%d-%d, Area = %.4f\n',tri(1),tri(2),tri(3),A);

px = 0.75; py = 0.75;
Nv = Nfun(px,py);
fprintf('N(%.2f,%.2f) = [%.4f %.4f %.4f]\n',px,py,Nv(1),Nv(2),Nv(3));

ue = u(tri); ve = v(tri);
uP = Nv'*ue;
vP = Nv'*ve;

% CST: deformaciones constantes
epsX  = (b'*ue)/(2*A);
epsY  = (c'*ve)/(2*A);
epsXY = (c'*ue + b'*ve)/(2*A);
eps_v = [epsX; epsY; epsXY];
sigma = D*eps_v;

fprintf('u(%.2f,%.2f) = %.4f mm   v(%.2f,%.2f) = %.4f mm\n',px,py,uP,px,py,vP);
fprintf('eps = [%.4f %.4f %.4f]\n',eps_v(1),eps_v(2),eps_v(3));
fprintf('sigma [Pa] = [%.4e %.4e %.4e]\n',sigma(1),sigma(2),sigma(3));
fprintf('--> Al ser CST, la deformacion y la tension son CONSTANTES en todo\n');
fprintf('    el elemento (6,5,8); en cualquier punto cercano a (0.75,0.75)\n');
fprintf('    que siga perteneciendo a ese elemento, eps y sigma son los mismos.\n');

%% ---------- Punto (1.5, 0.75) -> elemento cuadrangular (5,4,9,8) ----------
quad = [5 4 9 8];
xe = nodes(quad,:);
C = quadShapeBilineal(xe); % columnas: a,b,c,d de cada Ni

px = 1.5; py = 0.75;
basis   = [1 px py px*py];
Nv      = (basis*C)';                 % [N1;N2;N3;N4]
bcoef   = C(2,:)';                    % b_i (coef. de x)
ccoef   = C(3,:)';                    % c_i (coef. de y)
dcoef   = C(4,:)';                    % d_i (coef. de x*y)

fprintf('\nElemento cuadrangular %d-%d-%d-%d\n',quad(1),quad(2),quad(3),quad(4));
fprintf('N(%.2f,%.2f) = [%.4f %.4f %.4f %.4f]\n',px,py,Nv(1),Nv(2),Nv(3),Nv(4));

ue = u(quad); ve = v(quad);
uP = Nv'*ue;
vP = Nv'*ve;

% Deformaciones (dependen de x o y, bilineales)
epsX  = sum((bcoef + dcoef*py).*ue);
epsY  = sum((ccoef + dcoef*px).*ve);
epsXY = sum((ccoef + dcoef*px).*ue + (bcoef + dcoef*py).*ve);
eps_v = [epsX; epsY; epsXY];
sigma = D*eps_v;

fprintf('u(%.2f,%.2f) = %.4f mm   v(%.2f,%.2f) = %.4f mm\n',px,py,uP,px,py,vP);
fprintf('eps = [%.4f %.4f %.4f]\n',eps_v(1),eps_v(2),eps_v(3));
fprintf('sigma [Pa] = [%.4e %.4e %.4e]\n',sigma(1),sigma(2),sigma(3));
