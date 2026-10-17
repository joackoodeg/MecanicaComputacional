
% Parametros
rhoCp = 1;  % rho*Cp
k = 2;      % conductividad
G_func = @(x) 100*(4 - x);  % G(x)
Ti = 10;    % Temperatura en x=2 (Dirichlet)
T4 = 10;    % Condición inicial en todos los nodos
dx = 0.5;   % Espaciado en x
dt = 0.01;  % Paso de tiempo
xnode = 2:dx:4; % Nodos espaciales

% Cantidad de nodos
Nnod = length(xnode);

% Condicion inicial: T(x, 0) = 10
T = T4 * ones(Nnod, 1);

% Construccion de la matriz de rigidez espacial y vector de fuerza
sp = 0; % Bandera para trabajar sin sparse
[K, f] = FuncLineales(xnode, sp, k, @(x) G_func(x), Ti, T4);

% Inicializar el tiempo
t_final = 0.5;  % Tiempo final de simulación
num_steps = t_final / dt;

% Ciclo de tiempo para resolver
for t = 1:num_steps
    % Actualizacion en el tiempo usando el esquema explicito
    T_new = T + dt * (K * T + f);

    % Condiciones de borde
    T_new(1) = Ti;  % Dirichlet en x = 2
    T_new(end) = T_new(end-1) + 1 * dx / k;  % Neumann en x = 4

    % Actualizar la temperatura para el siguiente paso de tiempo
    T = T_new;

    % Visualizar la solución en cada paso de tiempo
    plot(xnode, T, 'o-', 'LineWidth', 2);
    title(sprintf('Distribución de temperatura en t = %.3f', t*dt));
    xlabel('Posición (x)');
    ylabel('Temperatura (T)');
    pause(0.1);
end

disp('Distribución final de temperatura:');
disp(T);

% --- Función FuncLineales (adaptada para este problema) ---
function [K, f] = FuncLineales(xnode, sp, k, G_func, Ti, Td)
    Nnod = length(xnode);
    Nelem = Nnod - 1;

    % Derivadas de funciones de forma
    dN_dx = @(h)[-1/h,1/h];

    % Inicialización de matrices locales
    Kele = zeros(Nelem, 2, 2);
    fele = zeros(Nelem, 2);

    % Ensamblaje de cada elemento
    for ele = 1:Nelem
        h = xnode(ele+1) - xnode(ele);
        Kele(ele, :, :) = (k / h) * (dN_dx(h)' * dN_dx(h));  % Matriz de rigidez
        fele(ele, :) = G_func((xnode(ele) + xnode(ele+1)) / 2) * (h / 2);  % Vector de fuerzas
    end

    % Ensamblaje global
    if sp
        rows = [];
        cols = [];
        coef = [];
    else
        K = zeros(Nnod, Nnod);
    end
    f = zeros(Nnod, 1);

    for ele = 1:Nelem
        in_gl = [ele, ele+1];  % Nodos globales involucrados
        if sp
            for il = 1:2
                ig = in_gl(il);
                for jl = 1:2
                    jg = in_gl(jl);
                    rows = [rows; ig];
                    cols = [cols; jg];
                    coef = [coef; Kele(ele, il, jl)];
                end
            end
            K = sparse(rows, cols, coef);
        else
            Klocal = reshape(Kele(ele, :, :), 2, 2);
            K(in_gl, in_gl) = K(in_gl, in_gl) + Klocal;
            f(in_gl) = f(in_gl) + fele(ele, :)';
        end
    end

    % Condiciones de borde Dirichlet en x=2
    K(1, :) = 0;
    K(1, 1) = 1;
    f(1) = Ti;

    % Condición Neumann en x=4 se aplicará fuera de esta función
end

