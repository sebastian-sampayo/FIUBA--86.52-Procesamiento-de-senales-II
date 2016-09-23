% Función para calcular el filtro de Kalman a partir de las matrices del
% modelo en variable de estado, condiciones iniciales, etc
%
% Utiliza función "kalman_filter" e itera en 'k'.
% Los parámetros que se le pasan son iguales salvo:
% y: Matriz de qxM, donde cada columna y(:,i) es una medición en el instante 'i'.
% k: Vector de instantes de tiempo.
% x_hat: es una Matriz de NxM, donde cada columna x_hat(:,i) es el vector de 
%         estados estimado en el instante 'k'.

function [x_hat, P_hat] = kalman_filter_process(x0, P0, A, B, C, Q, R, y, k)

    N = length(x0);
    x_hat(:,1) = x0;
    P_hat(:,:,1) = P0;
    
    % Poner función de 'k(i)' y meter dentro del ciclo.
    Ak = A;%(:,:,i-1);
    Bk = B;
    Ck = C;
    Qk = Q;
    Rk = R;    
    
    for i=2:length(k)
        [x_hat(:,i), P_hat(:,:,i)] = kalman_filter(x_hat(:,i-1), P_hat(:,:,i-1), Ak, Bk, Ck, Qk, Rk, y(:,i));
    end
end