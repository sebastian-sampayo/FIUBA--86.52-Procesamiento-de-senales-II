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
% Función para calcular la autocorrelación de un proceso AR1
% -------------------------------------------------------------------------- %
%
% r(k) = sigma^2 / (1 - a^2) * a^|k|
%
% Uso:
%   r = corr_AR1 (k, a, sigma)
%
% donde:
%  k: Donde se evalua la función
%  a: El coeficiente del AR1.
%  sigma: La varianza del ruido blanco 'x'. 

function r = corr_AR1 (k, a, sigma)
  r = sigma^2 / (1 - a^2) * a.^abs(k);
end
