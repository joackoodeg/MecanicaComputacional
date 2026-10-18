# Guia 1 - Ejercicio 1 - Item b
# rho*cp*dT/dt = k*T'' - c*T + G   en [0,2]
# T(0)=100, q(2)=0, rho*cp=0, k=1, c=1, G=0
# Solucion analitica: T(x) = 100*exp(-x)*(exp(2x)+exp(4))/(1+exp(4))
#
# rho*cp = 0 -> solo estado estacionario
# CB derecha Neumann -> se compara orden 1 y orden 2 (model.ordenCB, tp_dif_finitas)

clear all; close all; clc;

L1 = 0; L2 = 2;
T_ex = @(x) 100*exp(-x).*(exp(2*x)+exp(4))/(1+exp(4));

function [T, err] = resolver_b(xnode, T_ex, ordenCB)
  N = length(xnode);
  model.rho = 0;
  model.cp  = 0;
  model.k   = 1;
  model.c   = 1;
  model.G   = zeros(1,N);
  model.ordenCB = ordenCB;
  cb = [1 100 -1; 2 0 -1];   # Dirichlet izq, Neumann der (q=0)
  et = [0 0 0 0];            # estado estacionario
  T = tp_dif_finitas(xnode, model, cb, et);
  err = max(abs(T(:) - T_ex(xnode(:))));
endfunction

# --- 1) Malla uniforme N=21, CB de orden 1 y 2, comparacion con la analitica ---
N = 21;
xnode = linspace(L1, L2, N);
[T1, e1] = resolver_b(xnode, T_ex, 1);
[T2, e2] = resolver_b(xnode, T_ex, 2);

figure(1); clf; hold on;
xfino = linspace(L1, L2, 200);
plot(xfino, T_ex(xfino), 'r-', 'linewidth', 1.5);
plot(xnode, T1, 'gs', 'markersize', 7);
plot(xnode, T2, 'bx', 'markersize', 8);
legend('Analitica', 'DF - CB orden 1', 'DF - CB orden 2', 'location', 'northeast');
xlabel('x'); ylabel('T'); title('Ej 1b - estacionario, N=21'); grid on;

printf('\nItem b - N=%d uniforme: err CB orden 1 = %.3e | err CB orden 2 = %.3e\n\n', N, e1, e2);

# --- 2) Convergencia: malla uniforme (CB orden 1 y 2) y no uniforme (CB orden 2) ---
Ns = [5 10 20 40 80 160];
eu1 = zeros(size(Ns)); eu2 = eu1; en2 = eu1; hu = eu1; hn = eu1;
figure(2);  # plots internos de tp_dif_finitas (se cierra al final)
for i = 1:length(Ns)
  N = Ns(i);
  xu = linspace(L1, L2, N);
  xn = L1 + (L2-L1)*linspace(0, 1, N).^2;   # no uniforme (refinada a la izquierda)
  hu(i) = max(diff(xu)); hn(i) = max(diff(xn));
  [~, eu1(i)] = resolver_b(xu, T_ex, 1);
  [~, eu2(i)] = resolver_b(xu, T_ex, 2);
  [~, en2(i)] = resolver_b(xn, T_ex, 2);
endfor
close(2);

# orden estimado: p = log(e_i/e_{i+1}) / log(h_i/h_{i+1})
ord = @(e, h) [NaN, log(e(1:end-1)./e(2:end)) ./ log(h(1:end-1)./h(2:end))];
p_u1 = ord(eu1, hu); p_u2 = ord(eu2, hu); p_n2 = ord(en2, hn);

printf('%5s | %10s %10s %6s | %10s %6s | %10s %10s %6s\n', ...
       'N', 'h_unif', 'err_ord1', 'p', 'err_ord2', 'p', 'h_nounif', 'err_ord2', 'p');
for i = 1:length(Ns)
  printf('%5d | %10.3e %10.3e %6.2f | %10.3e %6.2f | %10.3e %10.3e %6.2f\n', ...
         Ns(i), hu(i), eu1(i), p_u1(i), eu2(i), p_u2(i), hn(i), en2(i), p_n2(i));
endfor
