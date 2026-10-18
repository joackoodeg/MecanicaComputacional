## FEM 1D con elementos lineales y cuadraticos combinados
## PDE: k*T'' + G(x) = 0
## Cambios respecto de ElementosFinitosCombinado.m:
##  1. El nodo inicial de cada elemento se lleva con un contador (antes se usaba
##     el indice del elemento, que falla si un cuadratico no es el ultimo).
##  2. El largo de cada elemento se calcula a partir de sus nodos.
##  3. Se guarda K antes de Dirichlet para calcular el flujo reactivo.
##  4. Se agrega la interpolacion en puntos arbitrarios (inciso d).
##  5. Neumann con signo explicito: q > 0 = flujo saliente -> F -= q.
##  6. Sin paquete symbolic: el vector de fuente se integra con Gauss
##     (3 puntos, exacto hasta grado 5) en el elemento master.

###### PARAMETROS DE ENTRADA ####################
x_ini = -1;                  # extremo izquierdo
L_ele = 1;                   # largo de cada elemento
k = 2;
G = @(x) 10*x.^3;            # fuente (puede depender de x)
tipo_func_forma = [0, 1];    # 0 = lineal, 1 = cuadratica (uno por elemento)

# Condiciones de borde (comentar las que no se usen)
nodo_dirichlet = 1;  T_dirichlet = 0;
nodo_robin = 4;      h_rob = 5;  T_inf = 50;
# nodo_neumann = 4;  q_neumann = 10;   # q > 0 sale del dominio

x_eval = [-0.25, 0.25];      # puntos donde interpolar (inciso d)
#################################################

###### NODOS Y CONECTIVIDAD #####################
Nelem = length(tipo_func_forma);
xnode = x_ini;
icone = cell(Nelem, 1);
nod = 1;                                  # nodo inicial del elemento actual
for e = 1:Nelem
    xa = xnode(end);
    if tipo_func_forma(e) == 0
        xnode = [xnode, xa + L_ele];
        icone{e} = [nod, nod+1];
        nod += 1;
    else
        xnode = [xnode, xa + L_ele/2, xa + L_ele];
        icone{e} = [nod, nod+1, nod+2];
        nod += 2;
    endif
endfor
Nnod = length(xnode);
disp('xnode:'); disp(xnode);

###### FUNCIONES DE FORMA EN EL MASTER [-1, 1] ##
N_lin  = @(psi) [(1 - psi)/2, (1 + psi)/2];
N_cuad = @(psi) [psi.*(psi - 1)/2, (1 - psi).*(1 + psi), psi.*(psi + 1)/2];
pospg = [-sqrt(15)/5, 0, sqrt(15)/5];     # 3 puntos de Gauss
pespg = [5/9, 8/9, 5/9];

###### MATRICES ELEMENTALES Y ENSAMBLE ##########
K = zeros(Nnod);
F = zeros(Nnod, 1);
for e = 1:Nelem
    nn = icone{e};
    xe = xnode(nn);
    he = xe(end) - xe(1);
    J  = he/2;                            # dx = J dpsi
    xc = (xe(1) + xe(end))/2;             # centro del elemento
    if tipo_func_forma(e) == 0
        Ke = (k/he) * [1 -1; -1 1];
        Nfun = N_lin;
    else
        Ke = (k/(3*he)) * [7 -8 1; -8 16 -8; 1 -8 7];
        Nfun = N_cuad;
    endif
    Fe = zeros(length(nn), 1);
    for pg = 1:length(pospg)
        psi = pospg(pg);
        xg  = xc + psi*J;                 # punto de Gauss en coordenadas reales
        Fe += pespg(pg) * G(xg) * Nfun(psi)' * J;
    endfor
    printf('\nElemento %d (nodos %s):\n', e, mat2str(nn));
    disp('Ke ='); disp(Ke);
    disp('Fe ='); disp(Fe);
    K(nn, nn) += Ke;
    F(nn) += Fe;
endfor
disp('K global (sin condiciones de borde):'); disp(K);
disp('F global (sin condiciones de borde):'); disp(F);

###### CONDICIONES DE BORDE #####################
# Primero las naturales (Robin, Neumann), despues Dirichlet
if exist('nodo_robin', 'var')
    K(nodo_robin, nodo_robin) += h_rob;
    F(nodo_robin) += h_rob * T_inf;
end
if exist('nodo_neumann', 'var')
    F(nodo_neumann) -= q_neumann;
end
K_sinD = K;  F_sinD = F;                 # copia para el flujo reactivo
if exist('nodo_dirichlet', 'var')
    K(nodo_dirichlet, :) = 0;
    K(nodo_dirichlet, nodo_dirichlet) = 1;
    F(nodo_dirichlet) = T_dirichlet;
end
disp('K global con condiciones de borde:'); disp(K);
disp('F global con condiciones de borde:'); disp(F);

###### SOLUCION #################################
T = K \ F;
disp('Temperaturas nodales:'); disp([xnode', T]);

###### FLUJO REACTIVO EN EL NODO DIRICHLET ######
if exist('nodo_dirichlet', 'var')
    Q = K_sinD(nodo_dirichlet, :) * T - F_sinD(nodo_dirichlet);
    printf('Flujo que entra por el nodo Dirichlet: %g (negativo = sale)\n', Q);
end

###### INTERPOLACION (INCISO d) #################
for p = 1:length(x_eval)
    xp = x_eval(p);
    for e = 1:Nelem
        nn = icone{e};  xe = xnode(nn);
        if xp >= xe(1) && xp <= xe(end)
            psi_p = (2*xp - xe(1) - xe(end)) / (xe(end) - xe(1));   # x -> psi
            if tipo_func_forma(e) == 0
                Np = N_lin(psi_p);
            else
                Np = N_cuad(psi_p);
            endif
            printf('T(%g) = %.4f  (elemento %d, N = %s)\n', xp, Np * T(nn), e, mat2str(Np, 4));
            break;
        endif
    endfor
endfor
