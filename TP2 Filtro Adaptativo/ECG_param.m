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
% Ejercicio de aplicación. Filtrado de señal de ECG.
% Con este código se barrieron los parámetros para elegir los valores
% -------------------------------------------------------------------------- %

clear all;
close all;

% ---- Parámetros ---- %
A = 1; % Amplitud de la señal de interferencia g(i)
theta = pi/4; % Fase de la señal de interferencia g(i)
f = 60; % [Hz] Frecuencia de la señal de interferencia g(i)
SNR = 20; % [dB] SNR entre el ruido 'v' del modelo y la interferencia g(i)
FS = 8192; % [Hz] Frecuencia de muestreo
FS_D = 500; % [Hz] Frecuencia de Downsampling

% Con estos parámetros converge para i=1000, y el error final es de -11dB
mu = .05;
M = 500;
% --
N_mu = 4;
mu = linspace(0.0001, 0.05, N_mu);

N_M = 3;
M = round(linspace(100, 1000, N_M));

s_file_name = 'datos/senal_ECG1.mat';
% ------------------- %
% Calculo la varianza del ruido del modelo de regresión lineal
sigma_v = sqrt(A^2 / (2*10^(SNR/10)));

load(s_file_name);
s = ECG_signal1'; % Señal del ECG
s = decimate(s, round(FS/FS_D));
N_s = length(s); % Largo de la señal
n = (0:N_s-1)'; 
T = 1/FS_D; % Período de muestreo
t = n*T;
g = A * sin(2*pi*f*t + theta); % Señal de interferencia, g(i)
u = 0.1 * A * sin(2*pi*f*t); % Señal de referencia, u(i)
v = normrnd(0, sigma_v, N_s, 1); % Ruido blanco
d = s + g + v;
    

% Para cada mu
for k = 1:N_mu
  % Para cada largo M de filtro:
  for j = 1:N_M
    w = zeros(M(j), 1);
    % Itero N_its veces el LMS
    N_its = N_s - M(j) + 1; % Cantidad de iteraciones del LMS
    for i = M(j):N_its
      u_i = u(i : -1 : i-M(j)+1);
      s_hat(i,1) = d(i) - w' * u_i;      
%       error(i, k, j) = mismatch(s(i), s_hat(i));
      error(i, k, j) = mismatch(s(1:i), s_hat(1:i));
      w = LMS(w, mu(k), u_i, d(i)); % i+1
    end
  end
end

% Gráfico del Mismatch
figure
hold all;
colors = {'b', [.1, .6, .1], 'r', 'm', 'k'};
q = 0;
for k = 1:N_mu
  for j = 1:N_M
    plot(10*log10(error(:, k, j)));
    q = q+1;
    legend_str{q} = sprintf('mu = %f, M = %i', mu(k), M(j));
  end
end

legend(legend_str, 'Location', 'NorthEast');
title('ECG - Diferencia entre señal real y estimada');
ylabel('Diferencia entre señal real y estimada, normalizada [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
ylim([-30, 1]);

% print('-dpng', 'Imagenes/ECG_param_2.png');
