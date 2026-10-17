pkg load symbolic
warning('off', 'all');

######PARAMETROS DE ENTRADA####################
x1 = -1;   # Extremo izquierdo
x2 = 1;    # Extremo derecho
h = 1;     # Tamaño de paso (largo de cada elemento)
k=2;
# Definimos un vector de funciones de carga para cada barra
G = { @(x) 10 *(x^3); @(x) 10 *(x^3)}; % Ejemplo de funciones
G_index=1; # Depende de x
##G=[30 30 30 30];
##G_index=0 #No depende de x
tipo_func_forma = [0, 1]; # 0 = función de forma lineal, 1 = función de forma cuadrática
######PARAMETROS DE ENTRADA####################


#######ARMAMOS XNODE###############
xnode = []; # Inicializar vector de nodos
# Construcción del vector de nodos
x_actual = x1;

for i = 1:length(tipo_func_forma)
    if tipo_func_forma(i) == 0
        # Elemento con funciones de forma lineales (solo nodos en los extremos)
        xnode = [xnode, x_actual, x_actual + h];
        x_actual += h;
    elseif tipo_func_forma(i) == 1
        # Elemento con funciones de forma cuadráticas (nodo intermedio)
        xnode = [xnode, x_actual, x_actual + h/2, x_actual + h];
        x_actual += h;
    endif
endfor
# Eliminar nodos repetidos
xnode = unique(xnode);
# Mostrar el vector de nodos
disp("Vector de nodos xnode:");
disp(xnode);
#############################

# Inicializar matrices globales
n_nodos = length(xnode);
K_global = zeros(n_nodos, n_nodos);
M_global = zeros(n_nodos, n_nodos);
F_global = zeros(n_nodos, 1);

# Definir variable simbólica
syms x;

# Recorrer los elementos y calcular sus matrices
for i = 1:length(tipo_func_forma)
    if tipo_func_forma(i) == 0
        # Elemento de función de forma lineal
        x1 = xnode(i);
        x2 = xnode(i + 1);

        K_el = ((k)/h) * [1 -1; -1 1];
        M_el = (h / 6) * [2 1; 1 2];

        N = [(x2 - x) / h; (x - x1) / h];

        F_e = zeros(2, 1);
        if G_index == 1
            F_e(1) = double(int(G{i}(x) * N(1), x, x1, x2));
            F_e(2) = double(int(G{i}(x) * N(2), x, x1, x2));
        endif

        # Ensamblaje
        K_global(i:i+1, i:i+1) += K_el;
        M_global(i:i+1, i:i+1) += M_el;
        F_global(i:i+1) += F_e;

    elseif tipo_func_forma(i) == 1
        # Elemento de función de forma cuadrática
        x1 = xnode(i);
        x2 = xnode(i+2);
        xm = xnode(i+1);

        K_el = (k / (3 * h)) * [7 -8 1; -8 16 -8; 1 -8 7];
        M_el = (h / 30) * [4 2 -1; 2 16 2; -1 2 4];

        N = [(x - xm) * (x - x2) / ((x1 - xm) * (x1 - x2));
             (x - x1) * (x - x2) / ((xm - x1) * (xm - x2));
             (x - x1) * (x - xm) / ((x2 - x1) * (x2 - xm))];

        F_e = zeros(3, 1);
        if G_index == 1
            F_e(1) = double(int(G{i}(x) * N(1), x, x1, x2));
            F_e(2) = double(int(G{i}(x) * N(2), x, x1, x2));
            F_e(3) = double(int(G{i}(x) * N(3), x, x1, x2));
        endif

        # Ensamblaje
        K_global(i:i+2, i:i+2) += K_el;
        M_global(i:i+2, i:i+2) += M_el;
        F_global(i:i+2) += F_e;
    endif
endfor

# Mostrar resultados
disp('Matriz de rigidez global K:');
disp(K_global);
disp('Matriz de masa global M:');
disp(M_global);
disp('Vector de fuerzas global F:');
disp(F_global);



#############################################
##### Condiciones de frontera################
#############################################

% Configuración de nodos para condiciones de frontera
nodo_dirichlet = 1; % Nodo correspondiente a x = -1
temp_dirlichet = 0;

nodo_robin = 4; % Nodo correspondiente a x = 1
h_rob = 5; % Coeficiente convectivo
T_inf = 50; % Temperatura del entorno

% nodo_neumann = 2; % Descomentar si necesitas aplicar Neumann en algún nodo
# flujo_neumann= 10;
% Condición de Dirichlet: T = 0 en x = -1
if exist('nodo_dirichlet', 'var')
    K_global(nodo_dirichlet, :) = 0; % Anular la fila correspondiente
    K_global(nodo_dirichlet, nodo_dirichlet) = 1; % Colocar 1 en la diagonal
    F_global(nodo_dirichlet) = temp_dirlichet; % Ajustar el vector de fuerzas
end

% Condición de Robin: h = 5, T_inf = 50 en x = 1
if exist('nodo_robin', 'var')
    K_global(nodo_robin, nodo_robin) += h_rob; % Sumar h al coeficiente diagonal
    F_global(nodo_robin) += h_rob * T_inf; % Sumar h * T_inf al vector de fuerzas
end

% Condición de Neumann (opcional): flujo = valor en un nodo
if exist('nodo_neumann', 'var')
    F_global(nodo_neumann) += flujo_neumann; % Ajustar el vector de fuerzas
end

#############################################
########### IMPRIMIR STENCIL ################
#############################################

disp('K global con condiciones de contorno')
K_global
disp('F global con condiciones de contorno')
F_global
disp('Matriz de masa global M:');
M_global

###### SOLUCIÓN DEL SISTEMA ####################
# Resolver el sistema de ecuaciones K_global * T = F_global
T = K_global \ F_global;

# Mostrar la solución
disp("Temperaturas en los nodos:");
disp(T);
#################################################
