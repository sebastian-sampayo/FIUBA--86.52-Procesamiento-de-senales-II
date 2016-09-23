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
% Función para estimar un filtro a través del algoritmo APA
% -------------------------------------------------------------------------- %
% 
% w(i+1) = w(i) + mu * A(i) * (A(i)'A(i))^-1 * conj(e(i))
% 
% Uso:
%   w = APA (w0, mu, u, d, K)
% donde:
%   w0: Estimación en el instante anterior
%   mu: Parámetro del algoritmo LMS
%   u:  Entrada en el instante anterior (largo M)
%   d:  Señal deseada en el instante anterior
%   w:  Estimación en el instante posterior
%   K:  Orden del APA
%
% Utiliza variables persistentes en memoria, con lo cual es necesario 
% inicializarlo antes de empezar a iterar:
%  APA(M, K)

function w = APA (w0, mu, u, d, K)
  delta = 0.0001;
  persistent A_;
  persistent d_;
  if (nargin == 2)
    M = w0; K = mu;
    A_ = zeros(M, K);
    d_ = zeros(K, 1);
    return;
  end
  A_ = [u, A_(:, 1:end-1)]; % Matriz de M x K
  d_ = [d; d_(1:end-1)]; % Vector columna
  
  e = d_' - w0' * A_;  % Vector fila
  w = w0 + mu * A_ * inv(A_'*A_ + delta*eye(K)) * e'; % Vector columna
end
