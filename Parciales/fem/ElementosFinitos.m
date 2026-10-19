pkg load symbolic
warning('off', 'all');

# Parámetros de entrada
x1=-1;
x2=1;
h=1;
xnode = x1:h:x2; % Discretización
k = 2; % Valor de rigidez


# Definimos un vector de funciones de carga para cada barra
G = { @(x) 30 *(x^3); @(x) 30 *(x^3)}; % Ejemplo de funciones
G_index=1; # Depende de x


##G=[30 30 30 30];
##G_index=0 #No depende de x


# Inicializar matrices y vectores
n_nodos = length(xnode);
K = zeros(n_nodos, n_nodos);
M = zeros(n_nodos, n_nodos);
F = zeros(n_nodos, 1);

# Definir variable simbólica
syms x;

# Calcular matrices de elementos
K_el = ((k)/h) * [1 -1; -1 1] % Matriz de rigidez del elemento
M_el = (h / 6) * [2 1; 1 2] % Matriz de masa del elemento

% Inicializar el vector de fuerzas
F_e = zeros(length(xnode)-1,2); % Cada columna corresponderá a un elemento

### Definir funciones de forma como vectores y calcular integrales
for i = 1:(length(xnode)-1)
    x1 = xnode(i)
    x2 = xnode(i + 1)

    % Funciones de forma
    N = [(x2 - x) / h; (x - x1) / h] % N para el elemento

  % Calcular las integrales y almacenarlas en el vector F_e

  #G DEPENDE DE x
    if G_index==1
      F_e(i,1) = double(int(G{i}(x) * N(1), x, x1, x2));
      F_e(i,2) = double(int(G{i}(x) * N(2), x, x1, x2));
    endif

  #G NO DEPENDE DE x
    if G_index==0
        F_e(i,1) = double(int(G(i) * N(1), x, x1, x2))
        F_e(i,2) = double(int(G(i) * N(2), x, x1, x2))
    endif
end

###############ENSANBLAMOS#############################################
#######################################################################
K_global = zeros(n_nodos,n_nodos);
% Ensamblar la matriz de rigidez global K
for i = 1:(length(xnode)-1)
    K_global(i:i+1, i:i+1) += K_el; % Sumar la matriz de rigidez del elemento a la matriz global
end

% Inicializar el vector de fuerzas global
F_global = zeros(n_nodos, 1);
% Ensamblar el vector de fuerzas
for i = 1:(length(xnode)-1)
    F_global(i:i+1) += F_e(i, :)'; % Sumar las fuerzas del elemento al vector global
end

M_global = zeros(n_nodos,n_nodos);
% Ensamblar la matriz de masa global M
for i = 1:(length(xnode)-1)
    M_global(i:i+1, i:i+1) += M_el; % Sumar la matriz de masa del elemento a la matriz global
end

disp('K global sin condiciones de contorno')
K_global
disp('F global sin condiciones de contorno')
F_global
disp('Matriz de masa global M:');
M_global

#############################################
##### Condiciones de frontera################
#############################################

% Configuración de nodos para condiciones de frontera
nodo_dirichlet = 1; % Nodo correspondiente a x = -1
temp_dirlichet = 0;

nodo_robin = 3; % Nodo correspondiente a x = 1
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

#############################################
########### INCISO B) PARCIAL 1 2024 ################
#############################################

T = [10; 39.848; 52.149; 54.173; 52.852]; % Condiciones iniciales
rho = 1;
cp= 1;
delta_t = 0.1; % Paso de tiempo

#### Cálculo de T_n_mas_1
# Para esquema implícito: (M + delta_t * K) * Tn+1 = M*T^n + delta_t * F
T_n_mas_1 = inv(M_global) * ( (delta_t)/(rho*cp) * (F_global-(K_global*T))+ M_global*T)
T_n_mas_1 = inv(M_global) * ((delta_t / (rho * cp)) * (F_global - K_global * T) + M_global * T)

### Calculo de T_n_mas_1
# Para esquema explicito
T_n_mas_1 = T + inv(M_global)*( (delta_t / (rho * cp)) * (K_global * T+ F_global))

#############################################
########### INCISO C) PARCIAL 1 2024 ################
#############################################

### TEMPERATURA EN x=3.2 interpolacion directa
###seleccionamos los nodos 3 y 4 ya que son los que estan entre 3.2
T3 = T_n_mas_1(3); #nodo3=3
T4 = T_n_mas_1(4); #nodo4=3.5

# Calcular T en x = 3.2 mediante interpolación lineal
x_eval = 3.2;
T_x32 = T3 + (x_eval - 3) / (3.5 - 3) * (T4 - T3);

# Mostrar el resultado
disp('Temperatura en x = 3.2 mediante interpolación lineal directa:')
disp(T_x32)

# TEMPERATURA EN x=3.2 interpolacion con funciones de forma
syms x;
N1 = (3.5 - x) / (3.5 - 3);  % Función de forma asociada a T_3
N2 = (x - 3) / (3.5 - 3);    % Función de forma asociada a T_4

# Evaluar las funciones de forma en x = 3.2
x_eval = 3.2;
N_eval = double(subs([N1; N2], x, x_eval));

# Seleccionar los valores nodales correspondientes
#seleccionamos los nodos 3 y 4 ya que son los que estan entre 3.2
T3 = T_n_mas_1(3); #nodo3 = 3
T4 = T_n_mas_1(4); #nodo4 = 3.5

# Calcular T en x = 3.2 usando funciones de forma
T_x3_2 = N_eval' * [T3; T4];

# Mostrar el resultado
disp('Temperatura en x = 3.2 usando funciones de forma:')
disp(T_x3_2)

###############################################
############# INCISO D) PARCIAL 1 2024 ################
###############################################
# Parámetros de entrada
k = 2;  # Conductividad térmica (valor dado en el problema)

# Definir las temperaturas en los nodos 3 y 4
T3 = T_n_mas_1(3)
T4 = T_n_mas_1(4)
dT_dx = (T4 - T3) / (3.5 - 3) # Calcular el gradiente de temperatura entre x = 3 y x = 3.5

# Calcular el flujo de calor en x = 3.2 usando la ley de Fourier
q_x3_2 = -k * dT_dx;
q= -k * (T_n_mas_1(4) - T_n_mas_1(3)) / 0.5

# Mostrar el resultado
disp('Flujo de calor en x = 3.2:')
disp(q_x3_2)

########################################################
############# OBTENER N SOLUCIONES EN EL TIEMPO #########
########################################################
##
##%Explicito
##deltaT_max = (rho * cp * h^2) / (2 * k);  % Calculamos el delta t para asegurar estabilidad, debe ser menor igual asi q restamos 0.5
##epsilon = 1e-5;  % Un valor pequeño para asegurar que deltaT sea menor
##deltaT = deltaT_max - epsilon;
##disp('Paso de tiempo ajustado para estabilidad (deltaT):');
##disp(deltaT);
##
### Condiciones iniciales
##T = [10; 39.848; 52.149; 54.173; 52.852];  % Vector de temperaturas iniciales
##n_nodos = length(T);  % Número de nodos
##n_tiempo = 10;  % Número de pasos de tiempo
##
### Matriz para almacenar las temperaturas en cada paso temporal
##T_tiempo = zeros(n_nodos, n_tiempo);
##T_tiempo(:, 1) = T;  % Inicializar con la condición inicial
##
### Evolución temporal utilizando el esquema explícito con deltaT ajustado
##for n = 1:(n_tiempo - 1)
##    % Aplicar la fórmula para obtener T^{n+1}
##    T_tiempo(:, n+1) = T_tiempo(:, n) + (deltaT / (rho * cp)) * (K_global * T_tiempo(:, n) + F_global);
##
##    % Imprimir la temperatura en el nodo x=3 (nodo 3) para cada paso temporal
##    disp(['Temperatura en el nodo x=3 para el tiempo n=', num2str(n), ':']);
##    disp(T_tiempo(3, n+1));  % Nodo 3 es el asociado a x=3
##end
