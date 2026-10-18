function [T] = elemFinitos(xnode, model, cb, et)

  % [T] = elemFinitos(xnode, model, cb, et)
  %
  % Resuelve, por el Metodo de Elementos Finitos (FEM) en 1D, el modelo
  % completo de transferencia de calor:
  %   rho*cp * dT/dt = k * d2T/dx2 - c*T + G
  %
  % xnode: vector con TODAS las coordenadas nodales (incluye nodos medios
  %        si se usan funciones de forma cuadraticas).
  %
  % model: struct con las constantes del modelo
  %   model.k       conductividad termica (escalar)
  %   model.c       coeficiente de reaccion (escalar)
  %   model.rho     densidad (escalar)
  %   model.cp      calor especifico (escalar)
  %   model.G       fuente: escalar constante o function handle @(x)
  %   model.tipoFF  1 = funciones de forma lineales (default)
  %                 2 = funciones de forma cuadraticas
  %                 (si es cuadratico, cada elemento consume 3 nodos
  %                  consecutivos de xnode: [2e-1, 2e, 2e+1])
  %                 Tambien puede ser un VECTOR con un valor por elemento
  %                 para mezclar lineales y cuadraticos, p. ej. [1 2]:
  %                 un lineal consume 2 nodos y un cuadratico 3, y los
  %                 elementos vecinos comparten el nodo extremo.
  %
  % cb: matriz 2x3 de condiciones de borde
  %   cb(1,:) = borde izquierdo (x=L1), cb(2,:) = borde derecho (x=L2)
  %   cb(:,1) = tipo: 1-Dirichlet, 2-Neumann, 3-Robin
  %   cb(:,2) = valor: T (Dirichlet), q (Neumann), h (Robin)
  %   cb(:,3) = -1 (Dirichlet/Neumann) o T_inf (Robin)
  %
  % et: esquema temporal
  %   et(1) = 0 estado estacionario
  %           1 Forward Euler  (explicito, theta=0)
  %           2 Backward Euler (implicito, theta=1)
  %           3 Crank-Nicolson (theta=0.5)
  %   et(2) = cantidad de pasos de tiempo (ignorado si et(1)==0)
  %   et(3) = tolerancia de corte por convergencia al estacionario
  %   et(4) = paso de tiempo dt (ignorado si et(1)==0)
  %
  % T: si et(1)==0, vector columna con la temperatura en cada nodo.
  %    si no, matriz (Nnod x Npasos+1) con la historia temporal,
  %    con T(:,1) la condicion inicial nula.

  Nnod = length(xnode);

  if isfield(model,'tipoFF')
    tipoFF = model.tipoFF;
  else
    tipoFF = 1;
  endif

  if isa(model.G,'function_handle')
    Gfun = model.G;
  else
    Gval = model.G;
    Gfun = @(x) Gval;
  endif

  % Conectividad de elementos segun el tipo de funcion de forma
  if ~isscalar(tipoFF)
    % Malla mixta: un tipo por elemento, el nodo inicial se lleva con un contador
    Nelem  = length(tipoFF);
    conect = cell(Nelem,1);
    nod = 1;
    for e = 1:Nelem
      conect{e} = nod:(nod + tipoFF(e));
      nod += tipoFF(e);
    endfor
    if nod ~= Nnod
      error('xnode tiene %d nodos pero tipoFF necesita %d', Nnod, nod);
    endif
  elseif tipoFF == 1
    Nelem = Nnod - 1;
    conect = zeros(Nelem,2);
    for e = 1:Nelem
      conect(e,:) = [e, e+1];
    endfor
  else
    Nelem = (Nnod - 1) / 2;
    conect = zeros(Nelem,3);
    for e = 1:Nelem
      conect(e,:) = [2*e-1, 2*e, 2*e+1];
    endfor
  endif

  % Ensamblaje
  Kdiff = zeros(Nnod,Nnod);
  S     = zeros(Nnod,Nnod); % matriz de forma (N_i*N_j), reusada para masa y reaccion
  F     = zeros(Nnod,1);

  for e = 1:Nelem
    if iscell(conect)
      nodos_e = conect{e};
    else
      nodos_e = conect(e,:);
    endif
    xe = xnode(nodos_e)(:);
    [Ke, Se, Fe] = elemento_local(xe, model.k, Gfun);
    Kdiff(nodos_e,nodos_e) += Ke;
    S(nodos_e,nodos_e)     += Se;
    F(nodos_e)              += Fe;
  endfor

  M    = (model.rho*model.cp) * S;
  Ktot = Kdiff + model.c * S;

  % Condiciones de Neumann y Robin (actuan sobre Ktot y F)
  for lado = 1:2
    if lado == 1
      nodo = 1;
    else
      nodo = Nnod;
    endif
    tipo = cb(lado,1);
    if tipo == 2
      q = cb(lado,2);
      F(nodo) -= q;
    elseif tipo == 3
      h    = cb(lado,2);
      Tinf = cb(lado,3);
      Ktot(nodo,nodo) += h;
      F(nodo)         += h*Tinf;
    endif
  endfor

  if et(1) == 0
    % ----- Estado estacionario -----
    [Ktot, F] = aplicar_dirichlet(Ktot, F, cb, Nnod);
    T = Ktot\F;
    plot(xnode, T, 'b-x');
    title('Estado estacionario');
  else
    % ----- Esquemas temporales (theta-method) -----
    switch et(1)
      case 1
        theta = 0; % Forward Euler
      case 2
        theta = 1; % Backward Euler
      case 3
        theta = 0.5; % Crank-Nicolson
      otherwise
        error('et(1) invalido: debe ser 0, 1, 2 o 3');
    endswitch

    nPasos = et(2);
    tol    = et(3);
    dt     = et(4);

    % Dirichlet: no se toca Ktot (se combina con M y reemplazar su fila
    % no impondria T = valor). Se resuelven solo los nodos libres y la
    % columna de cada nodo fijo, multiplicada por su valor, pasa al RHS.
    fijos = [];
    if cb(1,1) == 1, fijos(end+1) = 1;    endif
    if cb(2,1) == 1, fijos(end+1) = Nnod; endif
    libres = setdiff(1:Nnod, fijos);

    LHS  = M/dt + theta*Ktot;
    Rmat = M/dt - (1-theta)*Ktot;

    Tant = zeros(Nnod,1); % condicion inicial nula (T(x,0)=0)
    Tant = fijar_dirichlet(Tant, cb, Nnod);
    T = Tant;

    for n = 1:nPasos
      RHS  = Rmat*Tant + F;
      Tact = fijar_dirichlet(zeros(Nnod,1), cb, Nnod);
      if isempty(fijos)
        Tact = LHS\RHS;
      else
        Tact(libres) = LHS(libres,libres) \ (RHS(libres) - LHS(libres,fijos)*Tact(fijos));
      endif

      T = [T, Tact];
      err = norm(Tact-Tant,2)/norm(Tact,2);

      plot(xnode, Tact, '*-b');
      title(sprintf('paso: %d - error: %e', n, err));
      pause(0.01);

      if err < tol
        disp('Se alcanza el estado estacionario segun la tolerancia definida');
        return
      endif

      Tant = Tact;
    endfor
    disp('Se alcanza el maximo de pasos sin llegar al estado estacionario');
  endif

endfunction


function [Ke, Se, Fe] = elemento_local(xe, k, Gfun)
  % Matrices/vector de un elemento (2 o 3 nodos) via cuadratura de Gauss
  nnod_e = length(xe);
  h = xe(end) - xe(1);
  J = h/2;

  % 3 puntos de Gauss (exacto hasta polinomios de grado 5)
  pospg = [-sqrt(3/5), 0, sqrt(3/5)];
  pespg = [5/9, 8/9, 5/9];

  Ke = zeros(nnod_e,nnod_e);
  Se = zeros(nnod_e,nnod_e);
  Fe = zeros(nnod_e,1);

  for g = 1:length(pospg)
    eta = pospg(g);
    w   = pespg(g);

    if nnod_e == 2
      N      = [(1-eta)/2; (1+eta)/2];
      dNdeta = [-1/2; 1/2];
    else
      N      = [eta*(eta-1)/2; 1-eta^2; eta*(eta+1)/2];
      dNdeta = [eta-0.5; -2*eta; eta+0.5];
    endif

    dNdx = dNdeta / J;
    xg   = N' * xe;

    Ke += w*J * k * (dNdx * dNdx');
    Se += w*J * (N * N');
    Fe += w*J * Gfun(xg) * N;
  endfor
endfunction


function [K, F] = aplicar_dirichlet(K, F, cb, Nnod)
  if cb(1,1) == 1
    K(1,:) = 0; K(1,1) = 1; F(1) = cb(1,2);
  endif
  if cb(2,1) == 1
    K(Nnod,:) = 0; K(Nnod,Nnod) = 1; F(Nnod) = cb(2,2);
  endif
endfunction


function [T] = fijar_dirichlet(T, cb, Nnod)
  if cb(1,1) == 1
    T(1) = cb(1,2);
  endif
  if cb(2,1) == 1
    T(Nnod) = cb(2,2);
  endif
endfunction
