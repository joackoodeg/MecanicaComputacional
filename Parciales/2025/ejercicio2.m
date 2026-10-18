% Parcial 15/09/2025 - Problema 2: placa 1x1, MDF 2D, esquema explicito
%   k = 10, G = 500, c = 0, rho*cp = 1
%   Superior (y=1): T = 150        (Dirichlet)
%   Inferior (y=0): Robin h = 50, Tinf = 20
%   Izquierdo (x=0): T = 100       (Dirichlet)
%   Derecho (x=1): aislado         (Neumann q = 0)
%   dx = dy = 1/6, T inicial = 125, tolerancia 1e-5
close all; clear all; more off;

%% Malla estructurada 7 x 7 (numeracion: x por fuera, y por dentro, como en los ejercicios)
n = 7;  h = 1/(n-1);
[X, Y] = ndgrid(0:h:1, 0:h:1);              % X(ix,iy), Y(ix,iy)
id = @(ix,iy) (ix-1)*n + iy;                % numero de nodo
xnode = zeros(n*n, 2);
for ix = 1:n
  for iy = 1:n
    xnode(id(ix,iy),:) = [(ix-1)*h, (iy-1)*h];
  end
end
icone = [];
for ix = 1:n-1
  for iy = 1:n-1
    icone(end+1,:) = [id(ix,iy), id(ix+1,iy), id(ix+1,iy+1), id(ix,iy+1)];  % SW SE NE NW
  end
end

%% Condiciones de borde
% Esquinas: (0,1) y (1,1) -> 150 (superior); (0,0) -> 100 (izquierdo, Dirichlet manda sobre Robin);
%           (1,0) -> Robin abajo + Neumann a la derecha
DIR = []; NEU = []; ROB = [];
for ix = 1:n
  DIR(end+1,:) = [id(ix,n), 150];                 % superior
end
for iy = 1:n-1
  DIR(end+1,:) = [id(1,iy), 100];                 % izquierdo (sin la esquina superior)
end
for ix = 2:n
  ROB(end+1,:) = [id(ix,1), 50, 20, 1];           % inferior: h, Tinf, direccion 1 = Sur
end
for iy = 1:n-1
  NEU(end+1,:) = [id(n,iy), 0, 2];                % derecho: q = 0, direccion 2 = Este
end

%% Modelo
model.nnodes = n*n;
model.k   = 10*ones(n*n,1);
model.c   = zeros(n*n,1);
model.G   = 500*ones(n*n,1);
model.rho = 1;  model.cp = 1;
model.ts    = 0;          % [0] explicito, [1] implicito, [2] estacionario
model.maxit = 1e6;
model.tol   = 1e-5;
% Condicion inicial: 125 uniforme, con los nodos Dirichlet ya en su valor
% (si no, el codigo los lleva a su valor de a poco y falsea el tiempo).
model.PHI_n = 125*ones(n*n,1);
model.PHI_n(DIR(:,1)) = DIR(:,2);

%% Resolver
[PHI, Q] = fdm2d(xnode, icone, DIR, NEU, ROB, model);
dt = fdm2d_explicit_delta_t(xnode, model);
nsteps = columns(PHI) - 1;
Tf = PHI(:,end);
Qf = Q(:, end-1:end);                              % flujo del ultimo paso [Qx Qy]

%% 1) Tiempo para alcanzar el estacionario
printf('\n1) dt = %g, pasos = %d, t estacionario = %.4f s\n', dt, nsteps, nsteps*dt);

%% 2) Temperaturas
% (0.5,0.5) y (1,0) son nodos; (0.25,0.25) NO (cae entre nodos) -> interpolacion bilineal
TP = @(x,y) Tf(id(round(x/h)+1, round(y/h)+1));
printf('2) T(0.5,0.5) = %.4f\n', TP(0.5,0.5));
printf('   T(1,0)     = %.4f\n', TP(1,0));
Tgrid = reshape(Tf, n, n);                         % Tgrid(ix,iy)
T025 = interp2(0:h:1, 0:h:1, Tgrid', 0.25, 0.25, 'linear');
printf('   T(0.25,0.25) = %.4f  (interpolacion bilineal)\n', T025);

%% 3) Flujo en el centro
P = id(round(0.5/h)+1, round(0.5/h)+1);
printf('3) q(0.5,0.5) = [%.4f, %.4f], |q| = %.4f\n', Qf(P,1), Qf(P,2), norm(Qf(P,:)));

%% Grafico
fdm2d_graph_mesh(full(Tf), Qf, xnode, icone, 0, 0);
