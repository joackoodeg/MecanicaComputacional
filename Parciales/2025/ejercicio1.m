m.rho=1; m.cp=1; m.k=1; m.c=0;

cb = [1 0 -1;      % izquierda: tipo 1 = Dirichlet, T = 0
      3 1  0];     % derecha:   tipo 3 = Robin, h = 1, Tinf = 0

% 1) estacionario, 4 nodos (para verificar lo hecho a mano)
x4  = linspace(0,1,4);         % nodos 0, 1/3, 2/3, 1
m.G = 100*ones(1,4);           % fuente G = 100 en cada nodo
T4  = difFinitas(x4, m, cb, [0 0 0 0], zeros(1,4))

% 2) transitorio implicito, 41 nodos, lambda = 40
% et = [esquema, max pasos, tolerancia, dt]
%   2    : esquema implicito
%   1e5  : maximo de pasos de tiempo
%   1e-4 : tolerancia para cortar
%   dt   : paso de tiempo
x   = linspace(0,1,41);        % 41 nodos, dx = 1/40
m.G = 100*ones(1,41);
dt  = 40*(1/40)^2;             % lambda = k*dt/dx^2 = 40  ->  dt = 0.025
Ti  = 10*x.*(1-x).^2;          % condicion inicial T(x,0)
T   = difFinitas(x, m, cb, [2 1e5 1e-4 dt], Ti);

j = find(T(21,:) > 23, 1);     % fila 21 = x = 0.5; primera columna con T > 23
t = (j-1)*dt                   % -> 0.65

% 3) temperaturas finales
Tf = T(:,end);                 % ultima columna = estado estacionario alcanzado
Tf([11 21 31 41])'             % nodos en x = 0.25, 0.5, 0.75, 1

% 4) error relativo
Te  = 25*x'.*(3-2*x');         % solucion exacta en los nodos
err = norm(Tf-Te)/norm(Te)     % ||T - Tex|| / ||Tex||
