% -------------------------------------------------------------------------- %
% Facultad de Ingeniería de la Universidad de Buenos Aires
% Procesamiento de Señales II
% Trabajo Práctico 1: 
%   - Estimación de trayectorias utilizando el filtro de Kalman -
% 2° Cuatrimestre de 2015
%
% Sampayo, Sebastián Lucas
% Padrón: 93793
% e-mail: sebisampayo@gmail.com
%
% Ejercicio 7 
% -------------------------------------------------------------------------- %

clear all;
close all;

acel_file_name = 'archivos_tp/Acel.mat';
gyro_file_name = 'archivos_tp/Gyro.mat';
radar_file_name = 'archivos_tp/Radar.mat';
trayectoria_file_name = 'archivos_tp/trayectoria.mat';

load(acel_file_name);
load(gyro_file_name);
load(radar_file_name);
load(trayectoria_file_name);

tiempo_inercial = Acel(:,1);
t_actualizacion = Pradar(:,1);
T = tiempo_inercial(2) - tiempo_inercial(1); % Paso temporal
I = eye(2); % Identidad de R2x2
Z = zeros(2,2); % Matriz nula de R2x2
Ax_b = Acel(:,2) + .1;
Ay_b = Acel(:,3) - .3;
w_b = Gyro(:,2);
Pradar_x = Pradar(:,2);
Pradar_y = Pradar(:,3);
Vradar_x = Vradar(:,2);
Vradar_y = Vradar(:,3);

sigma_pos = 10; % [m]
sigma_vel = .1; % [m/s]
pos_0 = 100; % [m]
vel_0 = .2; % [m/2]
theta_0 = 40 * pi/180; % [rad]
% x0 = [pos_0;
%       pos_0;
%       vel_0;
%       vel_0;
%       1; %cos(theta_0/4);
%       0; %-sin(theta_0/4);
%       0.1;
%       -0.3]; 
% x0 = zeros(8,1);
x0 = .98*[Preal(1,2);
      Preal(1,3);
      Vreal(1,2);
      Vreal(1,3);
      cos(Theta(1,2));
      sin(Theta(1,2));
      .1
      -.3];
P0_z = [ cos(theta_0) , 0               ;
         0              , -sin(theta_0) ]; 
P0 = [ I*pos_0 , Z         , Z    ;
       Z         , I*vel_0 , Z    ;
       Z         , Z       , P0_z ];
P0 = [ P0     , [Z;Z;Z] ;
      [Z,Z,Z] , I/2       ];


% Proceso de muestreo y filtro de kalman
  k = tiempo_inercial;
  N = length(k);
  x_hat(:,1) = x0;
  P_hat(:,:,1) = P0;
  % Matrices constantes a lo largo del proceso
  Bk = eye(6);
  Bk = [Bk     ;
        [Z,Z,Z]];
  Qk = [ I*.2*(T^3)/6 , Z            , Z             ;
         Z            , I*.2*(T^2)/2 , Z             ;
         Z            , Z            , I*.2*(T^2)/2 ];
  %Qk = [ Qk     , [Z;Z;Z] ;
  %      [Z,Z,Z] , Z       ];

  Rk = [ I*sigma_pos^2 , Z             ;
         Z             , I*sigma_vel^2 ]; 

  % Mediciones:
  y = NaN*ones(4,N);
  tr = ismember(k, t_actualizacion);
  y(1,tr') = Pradar_x;
  y(2,tr') = Pradar_y;
  y(3,tr') = Vradar_x;
  y(4,tr') = Vradar_y;
  clear tr;
    
  % En cada iteración se recalculan las matrices del modelo con los nuevos datos
  for i=2:length(k)
    bx = x_hat(7,i-1);
    by = x_hat(8,i-1);
    Ma_k = [ Ax_b(i-1) - bx , Ay_b(i-1) - by  ;
             Ay_b(i-1) - by , -Ax_b(i-1) + bx ];
             
    Mz_k = [ 0         , w_b(i-1) ;
             -w_b(i-1) , 0        ];
%(:,:,i-1);
    z1 = x_hat(5,i-1);
    z2 = x_hat(6,i-1);
    Mb_k = [ -z1 , -z2 ;
              z2 , -z1 ];
    Ak = [ I , I*T , Ma_k*(T^2)/2 , Mb_k*(T^2)/2 ;
           Z , I   , Ma_k*T       , Mb_k*T       ;
           Z , Z   , I + Mz_k*T   , Z            ;
           Z , Z   , Z            , I            ];

    Ck = [ I , Z , Z , Z ;
           Z , I , Z , Z ];

    % r = rank([Ck * Ak; Ck*Ak*Ak]);
    % if (r ~= 8)
    %   disp('No Observable!');
    %   r
    %   break
    % end
    % Aplicamos el filtro de Kalman
    % f2 = @(x0) (...
    %   Ma_k * x0(5:6) );
    %   
    % f = @(x0) (...
    %   [ x0(1:2) + T*x0(3:4) + (T^2)/2*f2(x0);
    %     x0(3:4) + T*f2(x0);
    %     (I + Mz_k)*x0(5:6);
    %     x0(7:8) ]...
    %     );
    Af = [ I , I*T , Ma_k*(T^2)/2 , Z ;
           Z , I   , Ma_k*T       , Z       ;
           Z , Z   , I + Mz_k*T   , Z            ;
           Z , Z   , Z            , I            ];
    f = @(x0) ( Af * x0 );
    [x_hat(:,i), P_hat(:,:,i), e_hat(:,i)] = extended_kalman_filter(x_hat(:,i-1), P_hat(:,:,i-1), Ak, Bk, Ck, Qk, Rk, y(:,i), f);
    % [x_hat(:,i), P_hat(:,:,i), e_hat(:,i)] = kalman_filter(x_hat(:,i-1), P_hat(:,:,i-1), Ak, Bk, Ck, Qk, Rk, y(:,i));

  end

Pos_hat = x_hat(1:2, :)';
Vel_hat = x_hat(3:4, :)';
bias_hat = x_hat(7:8, :)';
% Para calcular el ángulo de orientación, utilizo la inversa del coseno
% y la combino con la estimación del seno para obtener el signo.
Theta_hat = 180/pi*acos(x_hat(5,:))' .* -sign(x_hat(6,:))';

% Realización particular del error = Valor real - Valor estimado
error_real = [Preal(:, 2:3), Vreal(:, 2:3), (Theta(:,2))] ...
             - [Pos_hat, Vel_hat, (Theta_hat)];

P_hat(:,:,end)
disp('Rango de matriz de observabilidad (8x8)');
rank([Ck * Ak; Ck*Ak*Ak])
% var_pos_x = P_hat(1,1,:);
% var_pos_y = P_hat(2,2,:);
% var_vel_x = P_hat(3,3,:);
% var_vel_y = P_hat(4,4,:);
% error_ruido = sum(abs(Pos(:,1)-Pos_gauss(:,1)).^2);
% error_estimacion = sum(abs(Pos(:,1)-Pos_hat(1,:)').^2);
% str = sprintf('Diferencia entre posición real y medición en x: %e', error_ruido);
% disp(str)
% str = sprintf('Diferencia entre posición real y estimación en x: %e', error_estimacion);
% disp(str)

% Gráficos
figure
hold all;
plot(k, Preal(:,2));
plot(k, Pos_hat(:,1));
xlabel('Tiempo [s]');
ylabel('Posición en x [m]');
legend('Posición real', 'Posición estimada');
title('Posición en x');
%print('-dpng', 'Pos_x.png');

figure
hold all;
plot(k, Vreal(:,2));
plot(k, Vel_hat(:,1));
xlabel('Tiempo [s]');
ylabel('Velocidad en x [m/s]');
legend('Velocidad real', 'Velocidad estimada');
title('Velocidad en x');
%print('-dpng', 'Vel_x.png');

figure
hold all;
plot(k, Theta(:,2));
plot(k, Theta_hat);
xlabel('Tiempo [s]');
ylabel('Ángulo de orientación \Theta [°]');
legend('Ángulo real', 'Ángulo estimado');
title('Ángulo de orientación');
%print('-dpng', 'theta.png');

figure
hold all;
plot(Preal(:,2), Preal(:,3))
plot(Pos_hat(:,1), Pos_hat(:,2))
plot(Pradar(:,2), Pradar(:,3))
xlabel('Posición en x [m]');
ylabel('Posición en y [m]');
legend('Posición real', 'Posición estimada', 'Posición medida');
title('Posición');
print('-dpng', 'ej7/Trayectoria.png');

subtitle_str = {{'Posición en x'}, {'Posición en y'}, ...
       {'Velocidad en x'}, {'Velocidad en y'}, ...
       {'cos(\theta)'}, {'-sin(\theta)'}, ...
       {'Sesgo en x'}, {'Sesgo en y'} };
       
figure
k_e = ~isnan(e_hat(1,:));
for i = 1:4
  subplot(2,2,i)
  plot(-300:300, xcorr(e_hat(i,k_e)))
  title(subtitle_str{i});
end
clear k_e;
suptitle('Correlación de Innovaciones')
print('-dpng', 'ej7/Correlacion_de_innovaciones.png');

figure
ylim_m = [0 8; 0 8;
          0 .08; 0 .08;
          -.05 .1; -.05 .1;
          -1 1; -1 1 ];
for i = 1:8
  subplot(4,2,i)
  z = reshape(P_hat(i,i,:), [length(k),1]);
  plot(k, z);
  xlabel('Tiempo [s]');
  ylabel('\sigma ^2');
  ylim(ylim_m(i,:));
  clear z;
  title(subtitle_str{i});
end
suptitle('Errores de estimación');
print('-dpng', 'ej7/Errores_de_estimacion.png');

figure
ylim_m = [-10 10; -10 10;
          -1.5 1.5; -1.5 1.5;
          -40 40];
for i = 1:5
  subplot(3,2,i)
  plot(k, error_real(:,i))
  xlabel('Tiempo [s]');
  ylabel('Error');
  ylim(ylim_m(i,:));
  title(subtitle_str{i});
  if i==5
    title('\theta');
  end
end
suptitle('Realización del error');
print('-dpng', 'ej7/Realizacion_del_error.png');

figure
hold all;
plot(k, bias_hat(:,1));
plot(k, bias_hat(:,2));
xlabel('Tiempo [s]');
legend('en x', 'en y');
title('Sesgo');
print('-dpng', 'ej7/Sesgos.png');

