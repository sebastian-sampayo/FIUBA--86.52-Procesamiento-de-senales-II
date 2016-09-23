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
% Función para calcular el Mismatch entre el filtro óptimo y su estimado
% -------------------------------------------------------------------------- %
%
% D = || h - w ||^2 
%     -------------
%       || h ||^2
%
% D = mismatch(h, w)
%
% donde:
%  h: Filtro óptimo
%  w: Filtro estimado

function D = mismatch(h, w)
  D = sum((h - w).^2) / sum(h.^2);
end
