% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales II
% Trabajo Práctico 2: 
%   - Filtrado Adaptativo -
% 2° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Simulación del Algoritmo LMS. Ejercicio 2 - Encontrar el valor de mu en que se vuelve inestable
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
N_avg = 2; % Cantidad de realizaciones Monte-Carlo para calcular esperanzas
N_its = 500; % Cantidad de iteraciones del LMS

a = 0.3; % Coeficiente del AR1 que genera 'u'
sigma_x = 1; % Varianza del ruido del AR1
SNR = 20; % [dB] SNR entre el ruido 'v' del modelo y la entrada 'u' filtrada

% mu = [.0005, .001, .005, .01];
mu = linspace(0.01, 0.014, 6);
% mu = logspace(log10(0.005), log10(0.014), 5);

h_file_name = 'datos/ir_short.mat';
% ------------------- %

load(h_file_name);
h = w0; % Filtro óptimo

M = length(h); % Largo del filtro

% Calculo la Matriz de correlación del proceso AR1
k = 0:(M-1);
r_u = corr_AR1 (k, a, sigma_x);
Ru = toeplitz(r_u);
clear k;
% Calculo la varianza del ruido del modelo de regresión lineal
sigma_v = sqrt(h' * Ru * h / 10^(SNR/10));
Jmin = sigma_v^2;

% Simulación Monte-Carlo para cada 'mu'
% Mismatch:
D = zeros(N_its+1, length(mu));
Misadj = zeros(length(mu), 1);
for k = 1:length(mu)
  D(1, k) = N_avg;
  % Para cada realización:
  for j = 1:N_avg
    % Calculo el proceso de entrada completo (largo: N_its + M )
    u = AR1(N_its + M, a, sigma_x);
    w = zeros(M, 1);
    % Itero N_its veces el LMS
    for i = M:N_its
      u_i = u(i : -1 : i-M+1);
      v_i = normrnd(0, sigma_v);
      d_i = h' * u_i + v_i;
      e_i = d_i - w' * u_i;
      w = LMS(w, mu(k), u_i, d_i); % i+1
      D(i+1, k) = D(i+1, k) + mismatch(h, w);
    end
    Misadj(k) = Misadj(k) + e_i^2;
  end
end
D = D / N_avg;
Misadj = Misadj / N_avg;
Misadj = (Misadj - Jmin) / Jmin;

lambda = eig(Ru);
D_teorico = zeros(size(D));
i = (0:N_its)';
for k = 1:length(mu)
  for j = 1:M
    D_teorico(:, k) = D_teorico(:, k) + mu(k)*Jmin / (2 - mu(k)*lambda(j)) ...
                       + (1 - mu(k)*lambda(j)).^(2*i) * (h(j).^2 - mu(k)*Jmin / (2 - mu(k)*lambda(j)));
  end
end
D_teorico = D_teorico / sum(h.^2);

figure
hold all;
i = (1:N_its+2-M)';
colors = {'b', [.1, .6, .1], 'r', 'm', 'c', 'k'};
for k = 1:length(mu)
  plot(i, 10*log10(D(M:end, k)), 'Color', colors{k})
  str = sprintf('La estimación del Misadjustment para mu = %f, es: %f\n', mu(k), Misadj(k));
  disp(str);
  legend_str{k} = sprintf('mu = %f', mu(k));
end
for k = 1:length(mu)
  plot(i, 10*log10(D_teorico(1:end-M+1, k)), '--', 'Color', colors{k})
  str = sprintf('El Misadjustment teórico para mu = %f, es: %f\n', mu(k), mu(k)*trace(Ru)/2);
  disp(str);
  legend_str{length(mu)+k} = sprintf('mu = %f - Teórico', mu(k));
end
legend(legend_str, 'Location', 'SouthWest');
title('Algoritmo LMS - Límite de estabilidad');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);

str = sprintf('El valor límite de estabilidad es : %f', 2/max(lambda));
disp(str);
% print('-dpng', 'Imagenes/LMS_stability.png');
