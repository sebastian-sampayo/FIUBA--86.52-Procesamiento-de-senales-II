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
% Simulación del Algoritmo APA. Ejercicios.
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
N_avg = 50; % Cantidad de realizaciones Monte-Carlo para calcular esperanzas
N_its = 20000; % Cantidad de iteraciones del LMS

a = 0.95; % Coeficiente del AR1 que genera 'u'
sigma_x = 1; % Varianza del ruido del AR1
SNR = 20; % [dB] SNR entre el ruido 'v' del modelo y la entrada 'u' filtrada

mu = linspace(0, 2, 10);;
K = [1, 2, 4, 8];

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
D = zeros(N_its+1, length(mu), length(K)); % Mismatch
for p = 1:length(K)
  for k = 1:length(mu)
    D(1, k, p) = N_avg;
    % Para cada realización:
    for j = 1:N_avg
      % Calculo el proceso de entrada completo (largo: N_its + M )
      u = AR1(N_its + M, a, sigma_x);
      w = zeros(M, 1);
      % Inicializo el APA (usa variables persistentes dentro para almacenar
      % los valores previos)
      APA(M, K(p));
      % Itero N_its veces el LMS
      for i = M:N_its
        u_i = u(i : -1 : i-M+1);
        v_i = normrnd(0, sigma_v);
        d_i = h' * u_i + v_i;
        e_i = d_i - w' * u_i;
        w = APA(w, mu(k), u_i, d_i, K(p)); % i+1        
        D(i+1, k, p) = D(i+1, k, p) + mismatch(h, w);
      end
    end
  end
end
D = D / N_avg;

% ---------------------- %
% Promedio D para ver mejor el valor límite
for p = 1:length(K)
  for k = 1:length(mu)
    D_filtrado(:, k, p) = filter(ones(1,128)/128, [1], D(:, k, p));
  end
end

% Gráfico del Mismatch filtrado
figure
hold all;
i = (1:N_its+2-M)';
q = 0;
colors = {'r', [.1, .6, .1], 'b', 'm', 'y', 'y', 'c', [.8, .3, .3], [.3, .8, .6], 'k'};
for p = 1:length(K)
  for k = 1:length(mu)
    plot(i, 10*log10(D_filtrado(M:end, k, p)), 'Color', colors{k})
%     q = (p-1)*length(K) + k;
    q = q + 1;
    legend_str{q} = sprintf('K = %i, mu = %.1f', K(p), mu(k));
  end
end

legend(legend_str, 'Location', 'NorthEast');
title('Algoritmo APA');
ylabel('Mismatch [dB]');
xlabel('Iteración i');
xlim([1, N_its]);


% Gráfico del Mismatch en función del mu
figure
hold all;
i = (1:N_its+1)';
for p = 1:length(K)
  plot(mu, 10*log10(D_filtrado(end, :, p)))
  legend_str{p} = sprintf('K = %i', K(p));
end

legend(legend_str, 'Location', 'NorthWest');
title('Algoritmo APA');
ylabel('Mismatch [dB]');
xlabel('mu');
% xlim([1, N_its]);

% print('-dpng', 'Imagenes/APA_mismatch_lim.png');
