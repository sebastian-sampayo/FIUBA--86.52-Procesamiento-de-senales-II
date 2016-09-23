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
% Filtro de Kalman - Algoritmo QR-Cholesky
% -------------------------------------------------------------------------- %
% Función para calcular el filtro de Kalman a partir de las matrices del
% modelo en variable de estado, condiciones iniciales, mediicones y 
% matrices de correlación de ruido. Realiza un solo paso (avance de tiempo),
% es decir, parte de x(0|0), y obtiene x(1|0) y x(1|1), devolviendo este último.
% Lo mismo con P.
%
% Uso:
%   [x_hat, P_hat, e_hat] = kalman_filter(x0, P0, A, B, C, Q, R, y);
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

function [x_hat, P_hat, e_hat] = kalman_filter_qr(x0, P0, A, B, C, Q, R, y)

    N = length(x0);
    p = length(y);
    q = size(Q,1);
    if (Q == zeros)
      Qc = Q;
    else
      Qc = chol(Q)'; % esto es para q no pinche el chol()
    end

    % Actualización (en caso de que haya medición)
    if (~isnan(y))
      M = [ chol(R)'          , C*chol(P0)'        , zeros(p,q) ;
            zeros(N,p)       , A*chol(P0)'        , B*Qc       ;
            -y'*inv(chol(R)) , x0'*inv(chol(P0)) , zeros(1,q) ];

      [QR_Q, QR_R] = qr(M');
      QR_R = QR_R';
      Z = QR_R(p+1:p+N, p+1:p+N);
      W2 = QR_R(end, p+1:p+N);
      X = QR_R(1:p, 1:p);
      W1 = QR_R(end, 1:p);

      P_hat = Z*Z';
      x_hat = Z*W2';
      e_hat = - X*W1';
    else
      % Predicción
      x_hat = A * x0; % (k|k-1)
      P_hat = A * P0 * A' + B * Q * B'; % (k|k-1)
      % x_hat = x0;
      % P_hat = P0;
      e_hat = NaN*y;
    end
end
