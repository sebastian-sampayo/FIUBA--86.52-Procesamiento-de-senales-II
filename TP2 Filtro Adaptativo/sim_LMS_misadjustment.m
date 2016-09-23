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
% Simulación del Algoritmo LMS. Ejercicio de Misadjustment
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
N_avg = 50; % Cantidad de realizaciones Monte-Carlo para calcular esperanzas
N_its = 20000; % Cantidad de iteraciones del LMS

a = 0.3; % Coeficiente del AR1 que genera 'u'
sigma_x = 1; % Varianza del ruido del AR1
SNR = 20; % [dB] SNR entre el ruido 'v' del modelo y la entrada 'u' filtrada

N_mu = 10;
mu = linspace(5e-4, 1e-2, N_mu);
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
% ---------------------- %

% Simulación Monte-Carlo para cada 'mu'
% Mismatch:
D = zeros(N_its+1, length(mu));
Misadj = zeros(length(mu), 1);
E_e_i = zeros(length(mu), 1);
error = zeros(N_its+1, length(mu));
for k = 1:length(mu)
  D(1, k) = N_avg;
  error(1, k) = N_avg;
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
      error(i+1, k) = error(i+1, k) + e_i.^2;
    end
    E_e_i(k) = E_e_i(k) + e_i^2;
  end
end
E_e_i = E_e_i / N_avg;
error = error / N_avg;
Misadj_i = (error - Jmin) / Jmin;
Misadj = (E_e_i - Jmin) / Jmin;

lambda = eig(Ru);

Misadj_teorico = mu*trace(Ru)/2;

% Grafico el Misadjustment en función del tiempo:
figure
hold all;
k = 1:length(mu);
for k = 1:length(mu)
    plot(Misadj_i(:,k));
    legend_str{k} = sprintf('mu = %f', mu(k));
end
legend(legend_str, 'Location', 'SouthWest');
title('Algoritmo LMS - Misadjustment');
ylabel('Misadjustment(i)');
xlabel('Iteración i');
xlim([1, N_its]);
% print('-dpng', 'Imagenes/LMS_misadjustment_it.png');

% Grafico el Misadjustment promediado:
for k = 1:length(mu)
  Misadj_filtrado(:, k) = filter(ones(1,256)/256, [1], Misadj_i(:,k));
end

figure
hold all;
k = 1:length(mu);
for k = 1:length(mu)
    plot(Misadj_filtrado(:,k));
    legend_str{k} = sprintf('mu = %f', mu(k));
end
legend(legend_str, 'Location', 'SouthWest');
title('Algoritmo LMS - Misadjustment promediado');
ylabel('Promedio móvil del Misadjustment(i)');
xlabel('Iteración i');
xlim([1, N_its]);
% print('-dpng', 'Imagenes/LMS_misadjustment_filtrado.png');


% Grafico el Misadjustment en función de mu
Misadj = Misadj_filtrado(end,:);
figure
hold all;
plot(mu, Misadj);
plot(mu, Misadj_teorico);
legend('Misadjustment medido', 'Misadjustment teórico', 'Location', 'NorthWest');
title('Algoritmo LMS - Misadjustment');
ylabel('Misadjustment');
xlabel('mu');

% disp(str);
% print('-dpng', 'Imagenes/LMS_misadjustment.png');
