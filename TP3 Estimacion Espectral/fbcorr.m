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
% Approach forward-backwards para estimar la matriz de autocorrelación.
% -------------------------------------------------------------------------- %
% Utiliza la función lscorr, con el método de la covarianza (LS).
% R: LxL

function R = fbcorr(y, L)
  N = length(y);

  R = lscorr(y, L, 'cov');

  J = rot90(eye(L));

  R = 1/2 * (R + J*R.'*J);
end
