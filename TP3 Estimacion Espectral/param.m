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
% Técnicas paramétricas
% -------------------------------------------------------------------------- %

clear;
close all;

print_figs = false;

% Colores para los gráficos
dblue = [0, 0, .8];
lblue = [.3, .3, .8];
dgreen = [0, .5, 0];
lgreen = dgreen + [.2, 0, .2];
red = [.8, 0, 0];
lred = red + [0, .2, .2];

sigma_v_from_SNR = @(snr) (sqrt(2 * 10.^(-snr/10)));
% ---------------------------------------------------- %
% Ejercicio 2

N = 20;
L = 3:10;
K = 2;
delta = 0.5;
%delta = delta * 10;
f1 = 0.2;
f2 = f1 + delta/N;
w1 = 2*pi*f1;
w2 = 2*pi*f2;
z1 = exp(1i*w1);
z2 = exp(1i*w2);
SNR = [30 10];
sigma_v = sigma_v_from_SNR(SNR);
phi1 = rand()*2*pi - pi;
phi2 = rand()*2*pi - pi;
n = (1:N)';

% Recorro para cada SNR
for i=1:2
  % Genero el ruido 
  v = [ randn(N, 1) + 1i*randn(N, 1) ] / sqrt(2) * sigma_v(i);

  % Genero la señal 'y'
  y = exp(1i*(w1*n + phi1)) + exp(1i*(w2*n + phi2)) + v;

  % Recorro para cada L
  for j = 1:length(L)
    [w_hat(:,j), z_hat, p_roots, lambda] = music(y, L(j), K);

    w_err(:,j) = 100 * abs(w_hat(:,j)-[w1;w2])./[w1;w2];
    var_err(j) = 100 * abs(lambda(1) - sigma_v(i)^2)/abs(lambda(1));


    % Errores
    buffer = '';
    str = sprintf('\n------------------------------------------------\n');
    buffer = strcat(buffer, str);
    str = sprintf('Parámetros:\n');
    buffer = strcat(buffer, str);
    str = sprintf('N = %d, L = %d, K = %d, delta = %.1f, SNR = %d dB\n\n', N, L(j), K, delta, SNR(i));
    buffer = strcat(buffer, str);
    str = sprintf('w1 = %.2f, w2 = %.2f \nhat{w1} = %.2f, hat{w2} = %.2f\n', w1, w2, w_hat(1,j), w_hat(2,j));
    buffer = strcat(buffer, str);
    str = sprintf('w1_{err} = %.2f %%%%, w2_{err} = %.2f %%%% \n', w_err(1,j), w_err(2,j));
    buffer = strcat(buffer, str);
    % str = sprintf('Autovalores de R: %.2f\n', lambda);
    % buffer = strcat(buffer, str);
    str = sprintf('sigma_v-err = %.2f %%%%\n', var_err(j));
    buffer = strcat(buffer, str);

    fprintf(buffer);

    % Solo para el ejercicio 2, j=1
    if (j==1)
      % Gráficos
      str = sprintf('Raíces del MUSIC. SNR = %i dB', SNR(i));
      zeroplot(z_hat(1), i, str, dblue)
      zeroplot(z_hat(2), i, str, lblue)
      plot(real(z1), imag(z1), '^', 'Color', dgreen);
      plot(real(z2), imag(z2), '^', 'Color', lgreen);

      % Prueba para exportar en formato latex
      texbuf = '';
      texbuf = strcat(texbuf,   '\\begin{table}[H]\n', ...
                                '\\begin{center}\n', ...
                                '%% | w | verdadero | estimado | error |\n', ...
                                '\\begin{tabular}{|c|c|c|c|}\n', ...
                                '\\hline\n', ...
                                'Parámetro & Valor verdadero & Valor estimado & Error \\%% \\\\ \n', ...
                                '\\hline\n');
      str = sprintf('$w_1$ & %.2f & %.2f & %.2f \\\\\\\\ \n \\\\hline \n', w1, w_hat(1,1), w_err(1,1));
      texbuf = strcat(texbuf, str);
      str = sprintf('$w_2$ & %.2f & %.2f & %.2f \\\\\\\\ \n \\\\hline \n', w2, w_hat(2,1), w_err(2,1));
      texbuf = strcat(texbuf, str);
      str = sprintf('$\\\\sigma_v^2$ & %.2g & %.2g & %.2f \\\\\\\\ \n \\\\hline \n', ...
                                     sigma_v(i)^2, lambda(1), var_err);
      texbuf = strcat(texbuf, str);
      str = sprintf('Autovalores de R & \\\\multicolumn{3}{c|}{ [');
      str2 = sprintf('%.2g \\\\; ', lambda);
      str = strcat(str, str2);
      str = strcat(str, '] } \\\\ \n \\hline \n');
      texbuf = strcat(texbuf, str);
      texbuf = strcat(texbuf, '', ...
                              '\\end{tabular}\n', ...
                              '\\end{center}\n', ...
                              '\\end{table}\n');

      fprintf(texbuf)
      if (print_figs)
        tex_file_name = sprintf('Tablas/param_2snr%d.tex',SNR(i));
        fig_file_name = sprintf('Imagenes/param_2snr%d.png',SNR(i));
        print('-dpng', fig_file_name);

        tex_file = fopen(tex_file_name, 'w');
        fprintf(tex_file, texbuf);
        fclose(tex_file);
      end
    end

    % Para SNR=30dB solo L=3
    if (i==1) 
      break
    end
  end


end

% Ejercicio 3
% Error en función de L
% Gráficos
figure(3)
hold all;
plot(L, w_err(1,:));
plot(L, w_err(2,:));
legend('w1', 'w2');
xlabel('L');
ylabel('%');
grid on;
title('Error porcentual en la estimación de w');
disp('--------------------');
w_err
if (print_figs)
  print('-dpng', 'Imagenes/param_3err.png');
end

% Periodograma
y_corr = xcorr(y, 'biased');
Nfft = 1024;
phi_y = abs(fft(y_corr, Nfft));
figure(4)
hold all;
w = linspace(0, 2, Nfft);
plot(w, dB(abs(phi_y)));
plot([w1,w1]/pi, ylim, 'k--');
plot([w2,w2]/pi, ylim, 'k--');
ylabel('PSD [dB]');
xlabel('w \pi');
title('Periodograma de la señal y');
if (print_figs)
  print('-dpng', 'Imagenes/param_3psd.png');
end
% print('-dpng', 'Imagenes/param_3psd.png');

% Error al aumentar L
figure(5)
hold all;
str = 'Error al aumentar L';
zeroplot(p_roots, 5, str, 'b');
plot(real(z1), imag(z1), '^', 'Color', dgreen);
plot(real(z2), imag(z2), '^', 'Color', dgreen);
if (print_figs)
  print('-dpng', 'Imagenes/param_3Lerr.png');
end
% print('-dpng', 'Imagenes/param_3Lerr.png');

