# Guia 1 - Ejercicio 1 - todos los items (a..g)
# rho*cp*dT/dt = k*T'' - c*T + G   en [L1,L2]
#
# Para cada item:
#   - estacionario con CB de orden 2 (malla uniforme), comparado con la analitica
#   - si hay Neumann/Robin: tambien con CB de orden 1
#   - malla no uniforme (CB orden 2)
#   - si rho*cp ~= 0: Forward Euler, Backward Euler y Crank-Nicolson desde T(x,0)=0
# Una figura por item (+ una con la evolucion temporal en los transitorios).

clear all; close all; clc;

N   = 41;      # nodos por item
tol = 1e-7;    # tolerancia de los esquemas temporales (error relativo entre pasos)
maxIt = 2e5;

# cb: [tipo valor extra] -> tipo 1 Dirichlet(T), 2 Neumann(q), 3 Robin(h, Tinf)
casos = {};
casos{end+1} = struct('nombre','a', 'L',[0 1],  'cb',[1 10 -1; 1 50 -1],   'rhocp',0, 'k',2, 'c',0,  ...
  'G',@(x) 100*ones(size(x)), 'Tex',@(x) -25*x.^2 + 65*x + 10);
casos{end+1} = struct('nombre','b', 'L',[0 2],  'cb',[1 100 -1; 2 0 -1],   'rhocp',0, 'k',1, 'c',1,  ...
  'G',@(x) zeros(size(x)),    'Tex',@(x) 100*exp(-x).*(exp(2*x)+exp(4))/(1+exp(4)));
casos{end+1} = struct('nombre','c', 'L',[1 5],  'cb',[2 2 -1; 1 0 -1],     'rhocp',0, 'k',1, 'c',0,  ...
  'G',@(x) 100*(x-3).^2,      'Tex',@(x) (-25*x.^4 + 300*x.^3 - 1350*x.^2 + 1906*x + 2345)/3);
casos{end+1} = struct('nombre','d', 'L',[0 1],  'cb',[1 10 -1; 3 0.2 50],  'rhocp',0, 'k',1, 'c',1,  ...
  'G',@(x) 50*ones(size(x)),  'Tex',@(x) -36.6897*exp(-x) - 3.3103*exp(x) + 50);
casos{end+1} = struct('nombre','e', 'L',[5 10], 'cb',[3 2 100; 1 50 -1],   'rhocp',1, 'k',2, 'c',0,  ...
  'G',@(x) x.^3,              'Tex',@(x) -x.^5/40 + 1225*x/3 - 4600/3);
casos{end+1} = struct('nombre','f', 'L',[0 1],  'cb',[1 0 -1; 3 2 10],     'rhocp',2, 'k',2, 'c',2,  ...
  'G',@(x) 75*ones(size(x)),  'Tex',@(x) -5/4*exp(-(x+1)).*(exp(x)-1).*(11*exp(x) + (11-30*exp(1))));
casos{end+1} = struct('nombre','g', 'L',[0 1],  'cb',[1 50 -1; 2 5 -1],    'rhocp',1, 'k',2, 'c',-2, ...
  'G',@(x) zeros(size(x)),    'Tex',@(x) 73.2433*sin(x) + 50*cos(x));

# resuelve un caso; manda el plot interno de tp_dif_finitas a la figura 100
function [T, err, nit, Tfull] = resolver(cs, xnode, ordenCB, et)
  figure(100);
  model.rho = cs.rhocp; model.cp = 1;
  model.k = cs.k; model.c = cs.c;
  model.G = cs.G(xnode);
  model.ordenCB = ordenCB;
  Tfull = tp_dif_finitas(xnode, model, cs.cb, et);
  nit = columns(Tfull) - 1;
  T = Tfull(:,end);
  err = max(abs(T - cs.Tex(xnode(:))));
endfunction

printf('\n%4s | %10s %10s %10s | %20s %20s %20s\n', 'item', 'err_ord2', 'err_ord1', 'err_nounif', ...
       'FE (nit, err)', 'BE (nit, err)', 'CN (nit, err)');
printf('%s\n', repmat('-', 1, 100));

for ic = 1:numel(casos)
  cs = casos{ic};
  L1 = cs.L(1); L2 = cs.L(2);
  xu = linspace(L1, L2, N);
  xn = L1 + (L2-L1)*linspace(0, 1, N).^2;     # no uniforme, refinada a la izquierda
  xfino = linspace(L1, L2, 300);
  et0 = [0 0 0 0];
  tieneNR = any(cs.cb(:,1) > 1);              # tiene Neumann o Robin?

  # --- estacionario ---
  [T2, e2] = resolver(cs, xu, 2, et0);
  [Tn, en] = resolver(cs, xn, 2, et0);
  if tieneNR
    [T1, e1] = resolver(cs, xu, 1, et0);
  else
    e1 = NaN;
  endif

  figure(ic); clf; hold on;
  plot(xfino, cs.Tex(xfino), 'r-', 'linewidth', 1.5);
  plot(xu, T2, 'bx', 'markersize', 7);
  leg = {'Analitica', 'Estacionario CB ord.2'};
  if tieneNR
    plot(xu, T1, 'gs', 'markersize', 5);
    leg{end+1} = 'Estacionario CB ord.1';
  endif

  # --- transitorio (solo si rho*cp ~= 0) ---
  sFE = '-'; sBE = '-'; sCN = '-';
  if cs.rhocp ~= 0
    dx = max(diff(xu));
    dt_crit = cs.rhocp*dx^2/(2*cs.k);
    dtFE = 0.9*dt_crit;                        # explicito: por debajo del limite de estabilidad
    dtI  = 20*dt_crit;                         # implicitos: paso mucho mayor
    [TFE, eFE, nFE] = resolver(cs, xu, 2, [1 maxIt tol dtFE]);
    [TBE, eBE, nBE] = resolver(cs, xu, 2, [2 maxIt tol dtI]);
    [TCN, eCN, nCN, Tt] = resolver(cs, xu, 2, [3 maxIt tol dtI]);
    figure(ic);   # resolver deja activa la figura 100
    plot(xu, TFE, 'm^', 'markersize', 5);
    plot(xu, TBE, 'cv', 'markersize', 5);
    plot(xu, TCN, 'ko', 'markersize', 5);
    leg(end+1:end+3) = {sprintf('FE (dt=%.2e, %d it)', dtFE, nFE), ...
                        sprintf('BE (dt=%.2e, %d it)', dtI, nBE), ...
                        sprintf('CN (dt=%.2e, %d it)', dtI, nCN)};
    sFE = sprintf('(%6d, %.2e)', nFE, eFE);
    sBE = sprintf('(%6d, %.2e)', nBE, eBE);
    sCN = sprintf('(%6d, %.2e)', nCN, eCN);

    # evolucion temporal con CN: algunos instantes hasta el estacionario
    figure(10+ic); clf; hold on;
    idx = unique(round(linspace(1, columns(Tt), 8)));
    for j = idx
      plot(xu, Tt(:,j), '-', 'linewidth', 1);
    endfor
    plot(xfino, cs.Tex(xfino), 'r--', 'linewidth', 1.5);
    xlabel('x'); ylabel('T'); grid on;
    title(sprintf('Ej 1%s - evolucion CN desde T=0 (t = %.2f ... %.2f)', cs.nombre, 0, (columns(Tt)-1)*dtI));
  endif

  figure(ic);
  legend(leg);
  xlabel('x'); ylabel('T'); grid on;
  title(sprintf('Ej 1%s - [%g,%g], rho*cp=%g k=%g c=%g, N=%d', cs.nombre, L1, L2, cs.rhocp, cs.k, cs.c, N));

  printf('%4s | %10.3e %10.3e %10.3e | %20s %20s %20s\n', cs.nombre, e2, e1, en, sFE, sBE, sCN);
endfor
close(100);
printf('\n(err = max|T_num - T_analitica| en los nodos; ord1/ord2 = orden de la CB Neumann/Robin)\n');
