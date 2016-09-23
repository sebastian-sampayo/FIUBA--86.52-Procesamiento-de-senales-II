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
% Función para pasar a db. 20*log10(x)
% -------------------------------------------------------------------------- %

function d = dB(x)
  d = 10*log10(x);
end
