% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales II
% Trabajo Práctico 3: 
%   - Estimación Espectral -
% 2° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Técnicas no paramétricas
% -------------------------------------------------------------------------- %

clear;
close all;

% Colores para los gráficos
blue = [0, 0, .8];
lblue = [.3, .3, .8];
dgreen = [0, .5, 0];
lgreen = dgreen + [.2, 0, .2];
red = [.8, 0, 0];
lred = red + [0, .2, .2];

print_figs = false;
% ---------------------------------------------------- %
% Ejercicio 1
N = 10000;
x = [ randn(N, 1) + i*randn(N, 1) ] / sqrt(2);

x_cov_b = xcorr(x, 'biased');
x_cov_unb = xcorr(x, 'unbiased');
k = -(N-1):(N-1);
k = k';
x_cov_teo = delta(k);


% ---------------------------------------------------- %
% Gráfico 1
figure(1)
hold all;
plot(k, abs(x_cov_unb));
plot(k, abs(x_cov_b));
plot(k, x_cov_teo);

legend('Estimador insesgado', 'Estimador sesgado', 'Teórica');
% CSCG = Circularly Symmetric Complex Gaussian 
title('Autocovarianza de ruido blanco CSCG'); 
xlabel('k');
ylabel('Autocovarianza');

if (print_figs)
  print('-dpng', 'Imagenes/no_param_1.png');
end
% ---------------------------------------------------- %
% Ejercicio 2
% a)
%
% phi_y(w) = |H(w)|^2 * phi_x(w)
% como phi_x(w) = 1
% phi_y(w) = |H(w)|^2

num_H = [ 1, -1.3817, 1.5632, -0.8843, 0.4096 ];
den_H = [ 1, 0.3544, 0.3508, 0.1736, 0.2401 ];
[phi_y, w] = freqz(num_H, den_H, 512, 'whole');
phi_y = abs(phi_y).^2;

% ---------------------------------------------------- %
% Gráfico 2
figure(2)
hold all;
plot(w/pi, dB(phi_y));

title('Densidad espectral de potencia de y, \phi_y'); 
xlabel('\omega * \pi');
ylabel('PSD [dB]');

if (print_figs)
  print('-dpng', 'Imagenes/no_param_2a.png');
end

% ---------------------------------------------------- %
% b)
% Usamos el estimador sesgado de la autocovarianza para estimar la PSD (periodograma).
J = 100;

phi_y64_mean = 0;
phi_y512_mean = 0;
for i = 1:J
  x64 =  ( randn(64, 1) + 1i*randn(64, 1)   ) / sqrt(2);
  x512 = ( randn(512, 1) + 1i*randn(512, 1) ) / sqrt(2);

  y64 = filter(num_H, den_H, x64);
  y512 = filter(num_H, den_H, x512);

  y64_cov = xcorr(y64, 'biased'); % length = 2N-1
  y512_cov = xcorr(y512, 'biased');
  % length = 2M-1, M = (length+1)/2
  Nfft_64 = length(y64_cov);
  Nfft_512 = length(y512_cov);

  phi_y64_hat(:,i) = abs(fft(y64_cov,   Nfft_64));
  phi_y512_hat(:,i) = abs(fft(y512_cov, Nfft_512 ));
  % Prueba con función periodogram(), da igual
  % Nfft_64 = 64;
  % Nfft_512 =512;
  %[phi_y64_hat(:,i), ~] = periodogram(y64, [], Nfft_64, 1, 'twosided');
  %[phi_y512_hat(:,i), ~] = periodogram(y512, [], Nfft_512, 1, 'twosided');

  phi_y64_mean = phi_y64_mean + phi_y64_hat(:,i);
  phi_y512_mean = phi_y512_mean + phi_y512_hat(:,i);

end
phi_y64_mean = phi_y64_mean/J;
phi_y512_mean = phi_y512_mean/J;

s_64 = 0;
s_512 = 0;
for i = 1:J
  s_64 = s_64 + (phi_y64_hat(:,i) - phi_y64_mean).^2 /(J-1);
  s_512 = s_512 + (phi_y512_hat(:,i) - phi_y512_mean).^2 /(J-1);
end

% ---------------------------------------------------- %
% Gráfico 3
figure(3)
hold all;
w64 = linspace(0, 2, Nfft_64);
w512 = linspace(0, 2, Nfft_512);
plot(w64, dB(abs(phi_y64_mean)), 'b', 'LineWidth', 2);
plot(w64, dB(abs(phi_y64_mean) + sqrt(s_64)), '--', 'Color', lblue);
aux = max(1e-2,abs(phi_y64_mean) - sqrt(s_64)); % por si da negativo, pongo -20dB
plot(w64, dB(aux), '--', 'Color', lblue);
plot(w512, dB(abs(phi_y512_mean)), 'Color', dgreen, 'LineWidth', 2);
plot(w512, dB(abs(phi_y512_mean) + sqrt(s_512)), '--', 'Color', lgreen);
aux = max(1e-2,abs(phi_y512_mean) - sqrt(s_512)); % por si da negativo, pongo -20dB
plot(w512, dB(aux), '--', 'Color', lgreen);
plot(w/pi, dB(phi_y), 'k--', 'LineWidth', 2);

legend('N = 64, E[\phi_y]', ...
       'N = 64, E[\phi_y] + \sigma_{\phi_y}', ...
       'N = 64, E[\phi_y] - \sigma_{\phi_y}', ...
       'N = 512, E[\phi_y]', ...
       'N = 512, E[\phi_y] + \sigma_{\phi_y}', ...
       'N = 512, E[\phi_y] - \sigma_{\phi_y}', ...
       '\phi_y Teórico');
title('Periodograma promedio y, \phi_y'); 
xlabel('\omega * \pi');
ylabel('PSD [dB]');

if (print_figs)
  print('-dpng', 'Imagenes/no_param_2b.png');
end

% ---------------------------------------------------- %
% c)
J = 100;
N = 256;
M = [N, N/4, N/16];
Nfft = 2*N;

phi_y_mean = zeros(Nfft, 3);
phi_y_hat = zeros(Nfft, J, 3);
s = zeros(Nfft, 3);
for i = 1:J
  x = ( randn(N, 1) + 1i*randn(N, 1) ) / sqrt(2);

  y = filter(num_H, den_H, x);

  y_cov = xcorr(y, 'biased'); % length = 2N-1
  % length = 2M-1, M = (length+1)/2

  for j = 1:3
    N_w_BT = 2*M(j)-1;
    k = -(M(j)-1):(M(j)-1);
    k = k';
    w_BT = 1 - abs(k)/M(j);
    aux = y_cov((N-N_w_BT/2):(N-1+N_w_BT/2));
    phi_y_hat(:,i, j) = abs(fft(aux.*w_BT,   Nfft));
    phi_y_mean(:, j) = phi_y_mean(:, j) + phi_y_hat(:,i,j);
  end
end
phi_y_mean = phi_y_mean/J;

for j = 1:3
  for i = 1:J
    s(:,j) = s(:,j) + (phi_y_hat(:,i,j) - phi_y_mean(:,j)).^2 /(J-1);
  end
end

% ---------------------------------------------------- %
% Gráfico 4
figure(4)
hold all;
color(:,:,1) = [blue; lblue; lblue];
color(:,:,2) = [dgreen; lgreen; lgreen];
color(:,:,3) = [red; lred; lred];

for j = 1:3
  w = linspace(0, 2, length(phi_y_mean(:,j)) );
  plot(w, dB(abs(phi_y_mean(:,j))), 'Color', color(1, :, j), 'LineWidth', 2);
  plot(w, dB(abs(phi_y_mean(:,j)) + sqrt(s(:,j))), '--', 'Color', color(2, :, j));
  aux = max(1e-2, abs(phi_y_mean(:,j)) - sqrt(s(:,j))); % por si da negativo, pongo -50dB
  plot(w, dB(aux), '--', 'Color', color(3, :, j));
end
w = linspace(0, 2, length(phi_y) );
plot(w, dB(phi_y), 'k--', 'LineWidth', 2);

legend('M = N, E[\phi_y]', ...
       'M = N, E[\phi_y] + \sigma_{\phi_y}', ...
       'M = N, E[\phi_y] - \sigma_{\phi_y}', ...
       'M = N/4, E[\phi_y]', ...
       'M = N/4, E[\phi_y] + \sigma_{\phi_y}', ...
       'M = N/4, E[\phi_y] - \sigma_{\phi_y}', ...
       'M = N/16, E[\phi_y]', ...
       'M = N/16, E[\phi_y] + \sigma_{\phi_y}', ...
       'M = N/16, E[\phi_y] - \sigma_{\phi_y}', ...
       '\phi_y Teórico');
title('Periodograma promedio y, \phi_y, N = 256'); 
xlabel('\omega * \pi');
ylabel('PSD [dB]');
%ylim([-20, 20]);
xlim([0, 2.5]);

if (print_figs)
  print('-dpng', 'Imagenes/no_param_2c.png');
% print('-dtex', 'Imagenes/no_param_2c.tex');
end


