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
% Función para estimar un filtro a través del algoritmo RLS
% -------------------------------------------------------------------------- %
% 
% w(i+1) = w(i) + k(n) + eps'(n)
% 
% Uso:
%   w = APA (w0, lambda, u, d)
% donde:
%   w0: Estimación en el instante anterior
%   lambda: Parámetro del algoritmo RLS
%   u:  Entrada en el instante anterior (largo M)
%   d:  Señal deseada en el instante anterior
%   w:  Estimación en el instante posterior
%
% Utiliza variables persistentes en memoria, con lo cual es necesario 
% inicializarlo antes de empezar a iterar:
%  RLS(M, lambda, alpha, sigma_u)

function w = RLS (w0, u, d, sigma_u)
%  lambda_ = 0.0001;
  persistent P_;
  persistent k_; 
  persistent lambda_;
  
  % Inicialización
  if (nargin == 4)
    M = w0;
    lambda_ = u;
    alpha = d;
    P_ = inv(sigma_u^2 * (1 - lambda_)^alpha * eye(M));
    k_ = zeros(M, 1);
    return;
  end
  
  k_ = 1/lambda_ * P_ * u / ( 1 + 1/lambda_ * u' * P_ * u); % k(n), P(n-1)
  eps = d - w0' * u;
    
  w = w0 + k_ * eps';
  
  % Actualizo el valor de P
  P_ = 1/lambda_ * P_ - 1/lambda_ * k_ * u' * P_; % P(n), P(n-1), k(n)
end
