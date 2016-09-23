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
% Simulación del Algoritmo RLS. Programa para ajustar los parámetros de los
% algoritmos y obtener la máximo velocidad de convergencia en cada caso.
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
N_avg = 2; % Cantidad de realizaciones Monte-Carlo para calcular esperanzas
N_its = 2000; % Cantidad de iteraciones del LMS

a = 0.9; % Coeficiente del AR1 que genera 'u'
sigma_x = 1; % Varianza del ruido del AR1
SNR = 30; % [dB] SNR entre el ruido 'v' del modelo y la entrada 'u' filtrada

N_param = 3;
% mu = .0075;
% mu = 0.001; % a = .9, SNR = 5
% mu = 0.0015; % a = 0, SNR = 5
% mu = linspace(.001, .002, N_param);
mu = linspace(.0007, .001, N_param);
K = [1, 2, 4];
% mu_apa = .1; % a = .9, SNR = 5
mu_apa = 1; % a = 0, SNR = 5
% lambda = .9995; % a = .9, SNR = 5
% lambda = .9999; % a = 0, SNR = 5
% lambda = [.91, .95, .99];
lambda = linspace(.999, .9999, N_param);
% lambda = linspace(.99, 1, N_param);
% alpha = -.5; % SNR = 5
alpha = 1; % SNR = 30

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
D_apa = zeros(N_its+1, N_param); % Mismatch
D_lms = zeros(N_its+1, N_param);
D_rls = zeros(N_its+1, N_param);
  for k = 1:3
    D_apa(1, k) = N_avg;
%     if (k == 1)
        D_lms(1, k) = N_avg;
        D_rls(1, k) = N_avg;
%     end
    % Para cada realización:
    for j = 1:N_avg
      % Calculo el proceso de entrada completo (largo: N_its + M )
      u = AR1(N_its + M, a, sigma_x);
      w_apa = zeros(M, 1);
      % Inicializo el APA (usa variables persistentes dentro para almacenar
      % los valores previos)
      APA(M, K(k));
%       if (k == 1)
          RLS(M, lambda(k), alpha, r_u(1));
          w_lms = zeros(M, 1);
          w_rls = zeros(M, 1);
%       end
      % Itero N_its veces el LMS
      for i = M:N_its
%         u_i = u(i+M-1 : -1 : i);
        u_i = u(i : -1 : i-M+1);
        v_i = normrnd(0, sigma_v);
        d_i = h' * u_i + v_i;
        w_apa = APA(w_apa, mu_apa, u_i, d_i, K(k)); % i+1        
        D_apa(i+1, k) = D_apa(i+1, k) + mismatch(h, w_apa);
%         if (k == 1)
            w_lms = LMS(w_lms, mu(k), u_i, d_i);
            w_rls = RLS(w_rls, u_i, d_i);

            D_lms(i+1, k) = D_lms(i+1, k) + mismatch(h, w_lms);
            D_rls(i+1, k) = D_rls(i+1, k) + mismatch(h, w_rls);
%         end
      end
    end
  end
D_apa = D_apa / N_avg;
D_lms = D_lms / N_avg;
D_rls = D_rls / N_avg;

% ---------------------- %
% Gráfico del Mismatch

% Gráfico de prueba de parámetros:
% APA
figure
hold all;
i = (1:N_its+2-M)';
for k = 1:N_param
  plot(i, 10*log10(D_apa(M:end, k)));
  legend_str{k} = sprintf('K = %i, mu = %.2f', K(k), mu_apa);
end
legend(legend_str, 'Location', 'NorthEast');
title('Algoritmos APA');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
% LMS
figure
hold all;
i = (1:N_its+2-M)';
for k = 1:N_param
  plot(i, 10*log10(D_lms(M:end, k)));
  legend_str{k} = sprintf('mu = %f', mu(k));
end
legend(legend_str, 'Location', 'NorthEast');
title('Algoritmos LMS');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
% RLS
figure
hold all;
i = (1:N_its+2-M)';
for k = 1:N_param
  plot(i, 10*log10(D_rls(M:end, k)));
  legend_str{k} = sprintf('lambda = %f, alpha = %f', lambda(k), alpha);
end
legend(legend_str, 'Location', 'NorthEast');
title('Algoritmos RLS');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
