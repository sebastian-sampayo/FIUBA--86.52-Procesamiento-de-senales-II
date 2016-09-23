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
% Función para estimar las frecuencias de sinusoides. Algoritmo MUSIC, raíces.
% -------------------------------------------------------------------------- %
% Algoritmo:
%   - Estimar R de LxL a partir de las N muestras de y.
%   - Diagonalizar el estimado de R y obtener un estimado de G
%   - Armar el vector p(z) = [1 z z^2 ... z^L-1]
%   - Se obtiene el polinomio de orden 2(L-1):
%         p(z^-1)T G G^H p(z)
%       que tiene 2(L-1) raíces, en pares recíprocos conjugados.
%   - Buscar las L-1 raíces que están en el círculo |z| <= 1
%   - De las L-1 raíces obtenidas seleccionar las K raíces más cercanas al círculo 
%     unidad {c_i}^K y estimar:
%         w_hat_k = fase{c_k}, k=1,...,K


function [w, z, p_roots, eigval] = music(y, L, K)
  N = length(y);

  % Estimo R utilizando forward-backwards (técnica de la covarianza, predictor lineal)
  R = fbcorr(y, L);

  % Obtengo G
  [G, D] = eigs(R, L-K, 'sm');
  eigval = eigs(R, L, 'sm');

  % Matriz del producto interno
  Q = G * G';

  p = zeros(2*(L-1) + 1, 1);
  for i = -(L-1) : (L-1)
    p(i+L) = sum(diag(Q, i));
  end

  p_roots = roots(p);
  z = p_roots( abs(p_roots) <= 1);
  % close all
  % zeroplot(p_roots, 5,'music', 'b');
  z = sort(z, 'descend');
  z = z(1:K);

  w = angle(z);
  % A la salida ordeno w y z según w creciente.
  [w, ix] = sort(w);
  z = z(ix);

  
end
