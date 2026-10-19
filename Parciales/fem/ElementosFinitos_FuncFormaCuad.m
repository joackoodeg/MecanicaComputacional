pkg load symbolic
warning('off', 'all');

# Parámetros de entrada
x1 = 0;   # Extremo izquierdo
x2 = 1;   # Extremo derecho
h = 1;  # Tamaño de paso (largo de cada elemento)

# Generar nodos automáticamente con nodos intermedios
xnode = [];
for i = x1:h:x2-h
    xnode = [xnode, i, i + h/2];
end
xnode = [xnode, x2];  # Agregar el último nodo

# Mostrar nodos generados
disp("Nodos generados:");
disp(xnode);

k = 2;  # Valor de rigidez

# Definir función de carga
G = { @(x) 10 * (x^3) };  # Ejemplo de función
G_index = 1;  # Depende de x

# Inicializar matrices y vectores
n_nodos = length(xnode);
K_global = zeros(n_nodos, n_nodos);
M_global = zeros(n_nodos, n_nodos);
F_global = zeros(n_nodos, 1);

# Definir variable simbólica
syms x;

# Matrices de elemento cuadrático
K_el = (k / (3 * h)) * [7 -8 1; -8 16 -8; 1 -8 7];
M_el = (h / 30) * [4 2 -1; 2 16 2; -1 2 4];

# Número de elementos finitos
n_elem = (length(xnode) - 1) / 2;

# Calcular vector de fuerzas por elemento
F_e = zeros(n_elem, 3);
for i = 1:2:(length(xnode)-2)
    x1 = xnode(i);
    x2 = xnode(i+2);
    xm = xnode(i+1);

    # Funciones de forma cuadráticas
    N = [(x - xm) * (x - x2) / ((x1 - xm) * (x1 - x2));
         (x - x1) * (x - x2) / ((xm - x1) * (xm - x2));
         (x - x1) * (x - xm) / ((x2 - x1) * (x2 - xm))]

    elem = (i + 1) / 2;  # Índice del elemento
    F_e(elem,1) = double(int(G{i}(x) * N(1), x, x1, x2))
    F_e(elem,2) = double(int(G{i}(x) * N(2), x, x1, x2))
    F_e(elem,3) = double(int(G{i}(x) * N(3), x, x1, x2))
end

# Ensamblaje de la matriz global
for i = 1:2:(length(xnode)-2)
    K_global(i:i+2, i:i+2) = K_global(i:i+2, i:i+2) + K_el;
    M_global(i:i+2, i:i+2) = M_global(i:i+2, i:i+2) + M_el;

    elem = (i + 1) / 2;
    F_global(i:i+2) = F_global(i:i+2) + F_e(elem, :)';
end


disp('K global sin condiciones de contorno')
disp(K_global)
disp('F global sin condiciones de contorno')
disp(F_global)
disp('Matriz de masa global M:')
disp(M_global)


#############################################
####### Condiciones de frontera################
###############################################
##
##% Configuración de nodos para condiciones de frontera
##nodo_dirichlet = 1; % Nodo correspondiente a x = -1
##nodo_robin = 5; h_rob=5 ; T_inf = 50; % Nodo correspondiente a x = 1
##
##% nodo_neumann = 2; % Descomentar si necesitas aplicar Neumann en algún nodo
##
##% Condición de Dirichlet: T = 0 en x = -1
##if exist('nodo_dirichlet', 'var')
##    K_global(nodo_dirichlet, :) = 0; % Anular la fila correspondiente
##    K_global(nodo_dirichlet, nodo_dirichlet) = 1; % Colocar 1 en la diagonal
##    F_global(nodo_dirichlet) = 0; % Ajustar el vector de fuerzas
##end
##
##% Condición de Robin: h = 5, T_inf = 50 en x = 1
##if exist('nodo_robin', 'var')
##    K_global(nodo_robin, nodo_robin) += h_rob; % Sumar h al coeficiente diagonal
##    F_global(nodo_robin) += h_rob * T_inf; % Sumar h * T_inf al vector de fuerzas
##end
##
##% Condición de Neumann (opcional): flujo = valor en un nodo
##if exist('nodo_neumann', 'var')
##    flujo = 10; % Valor del flujo
##    F_global(nodo_neumann) += flujo; % Ajustar el vector de fuerzas
##end
##
###############################################
############# IMPRIMIR STENCIL ################
###############################################
##
##disp('K global con condiciones de contorno')
##K_global
##disp('F global con condiciones de contorno')
##F_global
##disp('Matriz de masa global M:');
##M_global


