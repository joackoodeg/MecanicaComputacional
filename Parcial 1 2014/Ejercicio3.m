%% Parcial 1 - 2014 - Ejercicio 3
%% d/dx(kappa dphi/dx) + Q = 0,  x in [0,1]
%% phi(0) = 1 ;  dphi/dx(1) = 1
%% Metodo de Residuos Ponderados - Colocacion puntual, 3 funciones base
%% trigonometricas. Mismo espiritu que Tema 2/Ejemplos-RP (Ejemplo1.m,
%% Ejemplo2.m), pero con ponderacion por colocacion puntual en vez de
%% Galerkin, y derivadas analiticas (no requiere el paquete symbolic).
clear; close all; clc;

kappa = 1;
Q = -10;

% psi(x) satisface AMBAS condiciones de borde por si sola:
%   psi(0) = 1 ;  d(psi)/dx (1) = 1   ->   psi(x) = x + 1
psi   = @(x) x+1;

% Funciones base: Nm(x) = sin((2m-1)*pi*x/2)
%   Nm(0) = 0            -> no rompe phi(0)=1
%   dNm/dx(1) = 0         -> no rompe dphi/dx(1)=1
M = 3;
lambda = @(m) (2*m-1)*pi/2;
N      = @(m,x)  sin(lambda(m)*x);
d2N    = @(m,x) -lambda(m)^2*sin(lambda(m)*x);

% Puntos de colocacion: 3 puntos interiores equiespaciados en (0,1)
xl = [0.25 0.5 0.75];

% Residuo fuerte: R(x) = kappa*phi''(x) + Q
%   phi(x) = psi(x) + sum_m a_m*N_m(x),  psi'' = 0 (psi es lineal)
% Colocacion puntual: R(x_l) = 0, l = 1..3
%   kappa * sum_m [d2N_m(x_l)] a_m = -Q
K = zeros(M,M);
F = zeros(M,1);
for l = 1:M
    for m = 1:M
        K(l,m) = kappa*d2N(m,xl(l));
    end
    F(l) = -Q;
end

fprintf('Puntos de colocacion: '); disp(xl);
fprintf('Matriz K:\n'); disp(K);
fprintf('Vector F:\n'); disp(F);

a = K\F;
fprintf('Coeficientes a:\n'); disp(a');

phi_ap = @(x) psi(x) + a(1)*N(1,x) + a(2)*N(2,x) + a(3)*N(3,x);

fprintf('\nphi_aprox(x) = x + 1');
for m=1:M
    fprintf(' %+.4f*sin(%d*pi*x/2)', a(m), 2*m-1);
end
fprintf('\n');

%% Solucion exacta
% kappa*phi'' = -Q = 10  ->  phi'' = 10
% phi(x) = 5x^2 + C1*x + C2 ; phi(0)=1 -> C2=1 ; phi'(1)=1 -> 10+C1=1 -> C1=-9
phi_exact = @(x) 5*x.^2 - 9*x + 1;

%% Grafico
x_num = 0:0.01:1;
h = figure();
try
    graphics_toolkit(h,'gnuplot'); % mas robusto para exportar PNG sin pantalla
catch
end
plot(x_num,phi_ap(x_num),'r-','LineWidth',1.5); hold on;
plot(x_num,phi_exact(x_num),'b--','LineWidth',1.5);
plot(xl, phi_exact(xl), 'ko','MarkerFaceColor','k');
legend('Aproximada (colocacion, M=3)','Exacta','Puntos de colocacion');
xlabel('x'); ylabel('\phi(x)');
title('Ejercicio 3 - MRP por colocacion puntual');
grid on;
print(h,'Ejercicio3_grafico.png','-dpng');
fprintf('\nGrafico guardado en Ejercicio3_grafico.png\n');
