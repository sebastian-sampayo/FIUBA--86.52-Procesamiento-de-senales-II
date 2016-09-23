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
% Función para generar un proceso AR1 de condiciones iniciales nulas
% -------------------------------------------------------------------------- %
%
% Genera un proceso AR1 a partir del coeficiente y la varianza del ruido:
%  u(i) = a * u(i-1) + x(i)
%  con u(0) = 0
%
% Uso:
%  u = AR1(N, a, sigma)
%
% donde:
%  N: La cantidad de muestras generadas (largo del vector de salida 'u').
%  a: El coeficiente del AR1.
%  sigma: La varianza del ruido blanco 'x'. 

function u = AR1 (N, a, sigma)
  u = zeros(N, 1);
  % Genero el ruido blanco
  x = normrnd(0, sigma, [N 1]);
  % Calculo la primera muestra: (asumo u(0) = 0)
  % u(1) = x(1);
  % % Genero el proceso AR1
  % for i = 2:N
  %   u(i) = a * u(i-1) + x(i);
  % end
  % Por convolución:
%   h = a.^(0:N-1);
  % u = conv(h, x);
  u = filter([1], [1, -a], x);
%   u = u(1:N);
end
