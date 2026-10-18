function [T, res] = fem1d_mixto(xbordes, grado, k, c, G, bc, x_eval, verbose)
% FEM1D_MIXTO  Conduccion estacionaria 1D por elementos finitos,
% combinando elementos lineales y cuadraticos.
%
%   PDE:  -d/dx( k dT/dx ) + c T = G(x)      (equivale a k T'' - c T + G = 0)
%
% ENTRADAS
%   xbordes : extremos de los elementos, [x0 x1 ... xNe]  (Ne+1 valores)
%   grado   : grado de cada elemento, 1 = lineal, 2 = cuadratico  (Ne valores)
%   k, c    : conductividad y reaccion. Escalar o un valor por elemento.
%   G       : fuente. Escalar constante o handle @(x) (vectorizado: x.^3).
%   bc      : condiciones de borde, campos bc.izq y bc.der, cada uno uno de:
%               struct('tipo','dirichlet','valor',T0)
%               struct('tipo','neumann','q',q)        q > 0 = calor que SALE
%               struct('tipo','robin','h',h,'Tinf',Tinf)
%   x_eval  : (opcional) puntos donde interpolar T.
%   verbose : (opcional) true para imprimir matrices elementales y globales.
%
% SALIDAS
%   T   : temperaturas nodales
%   res : struct con xnode, icone, K y F (sin condiciones de borde),
%         Kcb y Fcb (con Robin/Neumann), Q (calor que ENTRA por cada
%         borde Dirichlet; negativo = sale) y Teval (T en x_eval).
%
% EJEMPLO (Ejercicio 2 de la guia):
%   bc.izq = struct('tipo','dirichlet','valor',0);
%   bc.der = struct('tipo','robin','h',5,'Tinf',50);
%   [T, r] = fem1d_mixto([-1 0 1], [1 2], 2, 0, @(x) 10*x.^3, bc, [-0.25 0.25], true);

  if nargin < 7, x_eval = []; end
  if nargin < 8, verbose = false; end

  Ne = length(grado);
  if isscalar(k), k = k*ones(1, Ne); end
  if isscalar(c), c = c*ones(1, Ne); end
  if isnumeric(G), Gval = G; G = @(x) Gval*ones(size(x)); end

  % ---------- Nodos y conectividad ----------
  xnode = xbordes(1);
  icone = cell(Ne, 1);
  nod = 1;
  for e = 1:Ne
    xa = xbordes(e);  xb = xbordes(e+1);
    if grado(e) == 1
      xnode = [xnode, xb];
      icone{e} = [nod, nod+1];
      nod = nod + 1;
    else
      xnode = [xnode, (xa+xb)/2, xb];
      icone{e} = [nod, nod+1, nod+2];
      nod = nod + 2;
    end
  end
  Nnod = length(xnode);

  % ---------- Elemento master [-1, 1] y Gauss (4 puntos) ----------
  pg = [-0.861136311594053, -0.339981043584856, 0.339981043584856, 0.861136311594053];
  wg = [ 0.347854845137454,  0.652145154862546, 0.652145154862546, 0.347854845137454];

  % ---------- Matrices elementales y ensamble ----------
  K = zeros(Nnod);  F = zeros(Nnod, 1);
  for e = 1:Ne
    nn = icone{e};
    xa = xbordes(e);  xb = xbordes(e+1);
    J  = (xb - xa)/2;                 % dx = J dpsi
    xc = (xa + xb)/2;
    nl = length(nn);
    Ke = zeros(nl);  Fe = zeros(nl, 1);
    for g = 1:length(pg)
      [N, dN] = forma(pg(g), grado(e));
      dNdx = dN / J;
      Ke = Ke + wg(g) * (k(e)*(dNdx'*dNdx) + c(e)*(N'*N)) * J;
      Fe = Fe + wg(g) * G(xc + pg(g)*J) * N' * J;
    end
    K(nn, nn) = K(nn, nn) + Ke;
    F(nn)     = F(nn) + Fe;
    if verbose
      fprintf('\nElemento %d (grado %d, nodos %s):\n', e, grado(e), mat2str(nn));
      disp('Ke ='); disp(Ke);
      disp('Fe ='); disp(Fe);
    end
  end
  res.K = K;  res.F = F;

  % ---------- Condiciones naturales (Robin / Neumann) ----------
  bordes = {bc.izq, bc.der};
  nodos_borde = [1, Nnod];
  fijos = [];  vfijos = [];
  for b = 1:2
    cb = bordes{b};  n = nodos_borde(b);
    switch lower(cb.tipo)
      case 'robin'
        K(n, n) = K(n, n) + cb.h;
        F(n)    = F(n) + cb.h * cb.Tinf;
      case 'neumann'
        F(n) = F(n) - cb.q;
      case 'dirichlet'
        fijos(end+1)  = n;
        vfijos(end+1) = cb.valor;
      otherwise
        error('Tipo de condicion desconocido: %s', cb.tipo);
    end
  end
  res.Kcb = K;  res.Fcb = F;
  if verbose
    disp('K global (sin condiciones de borde):'); disp(res.K);
    disp('F global (sin condiciones de borde):'); disp(res.F);
    disp('K global (con Robin/Neumann):'); disp(K);
    disp('F global (con Robin/Neumann):'); disp(F);
  end

  % ---------- Dirichlet: se resuelven solo los nodos libres ----------
  hay_robin = any(strcmpi({bc.izq.tipo, bc.der.tipo}, 'robin'));
  if isempty(fijos) && ~hay_robin && all(c == 0)
    error('Sin Dirichlet, sin Robin y sin reaccion: K es singular (T indefinida a menos de una constante).');
  end
  if isempty(fijos)
    T = K \ F;
    res.Q = [];
  else
    libres = setdiff(1:Nnod, fijos);
    T = zeros(Nnod, 1);
    T(fijos) = vfijos;
    T(libres) = K(libres, libres) \ (F(libres) - K(libres, fijos) * T(fijos));
    % ---------- Flujo reactivo en los nodos Dirichlet ----------
    res.Q = K(fijos, :) * T - F(fijos);    % calor que ENTRA por ese borde
  end

  % ---------- Interpolacion ----------
  res.Teval = zeros(size(x_eval));
  for p = 1:length(x_eval)
    xp = x_eval(p);
    e = find(xp >= xbordes(1:end-1) & xp <= xbordes(2:end), 1);
    psi = (2*xp - xbordes(e) - xbordes(e+1)) / (xbordes(e+1) - xbordes(e));
    N = forma(psi, grado(e));
    res.Teval(p) = N * T(icone{e});
  end

  res.xnode = xnode;  res.icone = icone;

  if verbose
    disp('Temperaturas nodales [x, T]:'); disp([xnode', T]);
    for i = 1:length(fijos)
      fprintf('Calor que entra por el nodo Dirichlet %d: %g (negativo = sale)\n', fijos(i), res.Q(i));
    end
    for p = 1:length(x_eval)
      fprintf('T(%g) = %.4f\n', x_eval(p), res.Teval(p));
    end
  end
end

function [N, dN] = forma(psi, grado)
% Funciones de forma y derivadas (respecto de psi) en el master [-1, 1]
  if grado == 1
    N  = [(1 - psi)/2, (1 + psi)/2];
    dN = [-1/2, 1/2];
  else
    N  = [psi*(psi - 1)/2, (1 - psi)*(1 + psi), psi*(psi + 1)/2];
    dN = [psi - 1/2, -2*psi, psi + 1/2];
  end
end
