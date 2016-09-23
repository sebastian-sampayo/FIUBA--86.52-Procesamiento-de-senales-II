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
% Algoritmo para estimar la autocorrelación y matriz de correlación mediante LS.
% -------------------------------------------------------------------------- %
%
% Armo un vector con n ceros al principio y n ceros al final.
% Genero la matriz Y de las diapositivas 16 y 17 de "Estimación Espectral-Ep.2.pdf"
% de forma genérica para N1, N2 y n.
% Luego calculo R, r y theta, según diapositivas 15 y 18.
% Parámetros:
% R: Estimado de la matriz de autocorrelación, nxn,
% r: Estimado del vector de autocorrelación, nx1 (a partir de r(1): [r(1) r(2) ...] ).
% theta: Estimado de los coeficientes del AR, nx1.
% n: Cantidad de coeficientes del AR.
% N1 y N2: Límites inferior y superior de la suma de términos de error a minimizar.

% n es la cantidad de coeficientes autoregresivos (1+el grado del polinomio denominador)
% es la dimensión de la matriz R, nxn.

function [R, r, theta] = lscorr(y, n, N1, N2)
  N = length(y);
  
  y_aux = [zeros(n,1); y; zeros(n,1)];

  % N1 y N2 predefinidos: Técnica de la Autocorrelación y de la Covarianza.
  if (nargin == 3)
    if (strcmp(N1, 'autocorr'))
      N1 = 1;
      N2 = N+n;
    elseif (strcmp(N1, 'cov'))
      N1 = n+1;
      N2 = N+1;
    end
  end

  Y = [];
  for i = (N1-1+n):(N2-1+n)
    Y = [Y; transpose(y_aux( i : -1 : (i-n+1) ))];
  end
  
  y_tilde = y_aux((N1+n) : (N2+n));

  R = 1/N * conj(Y)' * conj(Y);

  if (nargout > 1)
    r = 1/N * Y' * y_tilde;
  end
  if (nargout > 2)
    theta = -inv(N*R) * N*r;
  end
end

