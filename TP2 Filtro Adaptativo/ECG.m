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

ej4 = [true, false]; % Con esto en true hago el cambio abrupto de frecuencia
ambiente_ruidoso = [true, false]; % Con esto le indico que agregue ruido al modelo

s_file_name = 'datos/senal_ECG1.mat';
% ------------------- %
% Calculo la varianza del ruido del modelo de regresión lineal
sigma_v = sqrt(A^2 / (2*10^(SNR/10)));

load(s_file_name);
s = ECG_signal1'; % Señal del ECG
s = decimate(s, round(FS/FS_D)); % Hago un downsampling hasta 500Hz.
N_s = length(s); % Largo de la señal
n = (0:N_s-1)'; 
T = 1/FS_D; % Período de muestreo
t = n*T;
g = A * sin(2*pi*f*t + theta); % Señal de interferencia, g(i)
u = 0.1 * A * sin(2*pi*f*t); % Señal de referencia, u(i)
% Filtro Notch:
f_notch = [0 57/(FS_D/2) 59/(FS_D/2) 61/(FS_D/2) 63/(FS_D/2) 1];
a_notch = [1 1 0 0 1 1];
h_notch = firgr(M,f_notch,a_notch)';
h_notch = h_notch(1:M);
figure(1)
[H_notch,w_notch] = freqz(h_notch,1,512);
plot(f_notch,a_notch,w_notch/pi,abs(H_notch))
legend('Ideal','Notch FIR implementado')
xlabel('omega/\pi');
% print('-dpng', 'Imagenes/ECG_ej5_notch.png');

% Itero N_its veces el LMS
N_its = N_s - M + 1; % Cantidad de iteraciones del LMS        
% Hago las pruebas para ambiente ruidoso
% y cambio abrupto de frecuencia
for k = 1:2 % Cambio abrupto de frec
    for j = 1:2 % Ruido
        if ambiente_ruidoso(j)
            v = normrnd(0, sigma_v, N_s, 1); % Ruido blanco
        else
            v = zeros(N_s, 1);
        end
        % Recalculo las señales acá adentro, por si quedan los valores del caso de 70Hz
        f = 60;
        g = A * sin(2*pi*f*t + theta); % Señal de interferencia, g(i)
        u = 0.1 * A * sin(2*pi*f*t); % Señal de referencia, u(i)
        d = s + g + v;
        
        w = zeros(M, 1);
        for i = M:N_its
          % Si se configura el cambio abrupto de frecuencia:
          if (i == round(N_its/2)) && ej4(k)
              f = 70;
              g = A * sin(2*pi*f*t + theta); % Señal de interferencia, g(i)
              u = 0.1 * A * sin(2*pi*f*t); % Señal de referencia, u(i)
              d = s + g + v;
          end
          u_i = u(i : -1 : i-M+1);
          s_hat(i, j, k) = d(i) - w' * u_i;      
          % Calculo el error de toda la señal hasta el instante i
          error(i, j, k) = mismatch(s(1:i), s_hat(1:i, j, k));
          w = LMS(w, mu, u_i, d(i)); % i+1
          % Notch:
          d_i = d(i : -1 : i-M+1);
          s_hat_notch(i, j, k) = h_notch' * d_i;
          error_notch(i, j, k) = mismatch(s(1:i-M/2+1), s_hat_notch(M/2:i, j, k)); % desplazo en M/2 por el retardo de fase
        end
    end
end

% Gráficos
figure(3)
hold all;
figure(2)
hold all;
plot(s)
legend_str2{1} = 'Señal de ECG original';

i = (1:N_its)';
colors = {'b', [.1, .6, .1], 'r', 'k'};
q = 1;
str_con = {'Con', 'Sin'};
for k = 1:2 % Cambio abrupto de frec
  for j = 1:2 % Ruido
    figure(3)
    plot(10*log10(error(:, j, k)), 'Color', colors{k});
    plot(10*log10(error_notch(:, j, k)), 'Color', colors{k+2});
    figure(2)
    plot(s_hat(:, j, k), 'Color', colors{k});
    plot(s_hat_notch(M/2:end, j, k), 'Color', colors{k+2});% desplazo en M/2 por el retardo de fase
    legend_str{q} = sprintf('LMS - %s ruido. %s cambio abrupto de frecuencia.', str_con{j}, str_con{k});
    legend_str{q+1} = sprintf('Notch - %s ruido. %s cambio abrupto de frecuencia.', str_con{j}, str_con{k});
    legend_str2{q+1} = legend_str{q};
    legend_str2{q+2} = legend_str{q+1};
    q = q + 2;
  end
end
figure(3)
legend(legend_str, 'Location', 'NorthEast');
title('ECG - Ejercicio 5');
ylabel('Diferencia entre señal real y estimada, normalizada [dB]');
xlabel('Iteración i');
xlim([1, N_its]);
% ylim([-12, -8]);
% print('-dpng', 'Imagenes/ECG_ej5_error.png');

figure(2)
legend(legend_str2, 'Location', 'NorthEast');
title('ECG - Ejercicio 5 - Señal');
xlabel('Iteración i');
xlim([1, N_its]);
% print('-dpng', 'Imagenes/ECG_ej5_señal.png');

% print('-dpng', 'Imagenes/ECG_ej3.png');
% print('-dpng', 'Imagenes/ECG_ej4.png');
