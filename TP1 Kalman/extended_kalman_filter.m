% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales II
% Trabajo Práctico 1: 
%   - Estimación de trayectorias utilizando el filtro de Kalman -
% 2° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Filtro de Kalman Extendido
% -------------------------------------------------------------------------- %
% Función para calcular el filtro de Kalman Extendido a partir de las matrices del
% modelo en variable de estado, condiciones iniciales, mediicones y 
% matrices de correlación de ruido, y función no lineal f(x). Realiza un solo paso (avance de tiempo),
% es decir, parte de x(0|0), y obtiene x(1|0) y x(1|1), devolviendo este último.
% Lo mismo con P.
%
% Uso:
%   [x_hat, P_hat, e_hat] = kalman_filter(x0, P0, A, B, C, Q, R, y, f);
%
% x0: Vector de estado inicial (columna de N dimensiones).
% P0: Matriz de covarianza del estimador, inicial (de NxN).
% A, B, C: Matrices del modelo en variable de estado.
%     A: NxN, B:Nxp, C:qxN
% Q: Matriz de correlación del ruido de proceso (pxp)
% R: Matriz de correlación del ruido de medición (qxq)
% y: Medición nueva. Columna de q dimensiones. 
%   (en caso de no haber medición, utilizar NaN).
% e_hat: Innovación = yk - Ck * x_hat(k|k-1)
% f: Función no lineal de transición de estados.

function [x_hat, P_hat, e_hat] = kalman_filter(x0, P0, A, B, C, Q, R, y, f)

    N = length(x0);

    % Predicción
    x_hat = f(x0); % (k|k-1)

    P_hat = A * P0 * A' + B * Q * B'; % (k|k-1)
    
    e_hat = NaN*y;

    % Actualización (en caso de que haya medición)
    if (~isnan(y))
        % Calculo la innovación para checkear que el filtro anda bien
        e_hat = y - C*x_hat;
        % -----
        K = P_hat * C' * inv(R + C * P_hat * C'); % (Kk)

        x_hat = x_hat + K * (y - C * x_hat); % (k|k)

        P_hat = (eye(N) - K * C) * P_hat; % (k|k)
    end
end
