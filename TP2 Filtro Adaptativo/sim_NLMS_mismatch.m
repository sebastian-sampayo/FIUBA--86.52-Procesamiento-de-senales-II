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
% Simulación del Algoritmo NLMS. Ejercicios.
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
N_avg = 50; % Cantidad de realizaciones Monte-Carlo para calcular esperanzas
N_its = 20000; % Cantidad de iteraciones del LMS

a = 0.9; % Coeficiente del AR1 que genera 'u'
sigma_x = 1; % Varianza del ruido del AR1
SNR = 20; % [dB] SNR entre el ruido 'v' del modelo y la entrada 'u' filtrada

% mu = linspace(0.0001, 0.001, 5);
mu = [0.1 1];

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
D = zeros(N_its+1, length(mu)); % Mismatch
Misadj = zeros(length(mu), 1); % Misadjustment
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
%       w = LMS(w, mu(k), u_i, d_i); % i+1
      w = NLMS(w, mu(k), u_i, d_i); % i+1
      D(i+1, k) = D(i+1, k) + mismatch(h, w);
    end
  end
end
D = D / N_avg;

% Gráfico del Mismatch
figure
hold all;
i = (1:N_its+2-M)';
colors = {'b', [.1, .6, .1], 'r', 'm', 'k'};
for k = 1:length(mu)
  plot(i, 10*log10(D(M:end, k)), 'Color', colors{k})
  legend_str{k} = sprintf('mu = %f', mu(k));
end

legend(legend_str, 'Location', 'NorthEast');
title('Algoritmo NLMS');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
ylim([-30, 1]);

% print('-dpng', 'Imagenes/LMS_stability_0_9_2.png');
% print('-dpng', 'Imagenes/NLMS_mismatch.png');
