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
% Función para estimar un filtro a través del algoritmo LMS
% -------------------------------------------------------------------------- %
% 
% w(i+1) = w(i) + mu * u(i) * (conj(d(i)) - u(i)' * w(i))
% 
% Uso:
%   w = LMS (w0, mu, u, d)
% donde:
%   w0: Estimación en el instante anterior
%   mu: Parámetro del algoritmo LMS
%   u:  Entrada en el instante anterior
%   d:  Señal deseada en el instante anterior
%   w:  Estimación en el instante posterior

function w = LMS (w0, mu, u, d)
  e = d - w0' * u;
  w = w0 + mu * u * e';
end
