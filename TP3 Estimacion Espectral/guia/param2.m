clc
close all
clear

N = 20;
L = 8;

% Frecuencias 
K = 2 ;
f1 = 0.2;
delta = 0.5;
f2 = f1 + delta / N;

% Fases
phi1 = rand(1) * 2 * pi - pi;
phi2 = rand(1) * 2 * pi - pi;
% phi2 = -0.2*pi;
% phi1 = 0.5*pi;

% Frecuencias angulares 
w1 = 2 * pi * f1 ;
w2 = 2 * pi * f2 ;

% Ruido 
SNR = 30 ;   
varE = 1 ;  % varianza/energia de UNA exponencial compleja
varV = K*varE * 10^(-SNR/10);   % SRN = 10*log10( K*varE / varV );
v = (randn (N,1) + 1i * randn (N,1))* sqrt(varV/2);

% Señal de entrada. 
k = ( 0:N-1 )';
y = exp(1i * (w1*k  + phi1)) + exp(1i * (w2*k+ phi2)) + v;

R =  cov_FB(y, L);   % Estimo matriz de covarianza
[aves , avas] = eig(R);    % Calculo autovalores
[aves , avas] = sort_eig( aves, avas );  % Los ordeno.

[ w_est, zeros_est] = music( y, L , K );
w1_est = w_est(1);
w2_est = w_est(2);

e_w1 = abs( w1 - w1_est)/w1*100 ;
e_w2 = abs( w2 - w2_est)/w2*100 ;

fprintf('\n\nPisarenko\n');
fprintf('Parametros Iniciales: \n SNR = %ddB\n',SNR );
fprintf(' Delta = %.2f \n w1 = %.3f \n w2 = %.3f \n\n', delta,w1,w2);
fprintf('Frecuencias estimadas:\n w1_est = %.3f  \n w2_est = %.3f\n\n',w1_est,w2_est);
fprintf('Errores relativos de estimación:\n e_w1 = %.1f %% \n e_w2 = %.1f %%\n\n',e_w1,e_w2); 
fprintf('Autovalores Matriz de Convarianza:\n');
disp( avas' );

% Graficos 
figure, zplane(zeros_est,[])
grid on;
title('Diagrama de ceros. Music L = 3. SNR = 30dB');
%print('-dpng','diagramaDeCeros-L3.png');
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Vario la SNR 

fprintf('\n\nVario la SNR para analizar como varia la estimación y los autovalores de la matriz R\n');
% Ruido 
SNR = 10 ;   
varE = 1 ;  % varianza/energia de UNA exponencial compleja
varV = K*varE * 10^(-SNR/10);   % SRN = 10*log10( K*varE / varV );
v = (randn (N,1) + 1i * randn (N,1))* sqrt(varV/2);

% Señal de entrada. 
k = ( 0:N-1 )';
y = exp(1i * (w1*k  + phi1)) + exp(1i * (w2*k+ phi1)) + v;


R =  cov_FB(y, L);   % Estimo matriz de covarianza
[aves , avas] = eig(R);    % Calculo autovalores
[aves , avas] = sort_eig( aves, avas );  % Los ordeno.

[ w_est, zeros_est] = music( y, L , K );
w1_est = w_est(1);
w2_est = w_est(2);

e_w1 = abs( w1 - w1_est)/w1*100 ;
e_w2 = abs( w2 - w2_est)/w2*100 ;

fprintf('Parametros Iniciales: \n SNR = %ddB\n',SNR );
fprintf(' Delta = %.2f \n w1 = %.3f \n w2 = %.3f \n\n', delta,w1,w2);
fprintf('Frecuencias estimadas:\n w1_est = %.3f  \n w2_est = %.3f\n\n',w1_est,w2_est);
fprintf('Errores relativos de estimación:\n e_w1 = %.1f %% \n e_we = %.1f %%\n\n',e_w1,e_w2); 
fprintf('Autovalores Matriz de Convarianza:\n');
disp( avas' );
fprintf('Los autovalores comienzan a alejarse del cero\n');

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Punto 3. Music con L variable de 3 a 10.
clear all

fprintf('\n\nPunto 3. Music con L variable de 3 a 10.\n');
N = 20;
% Frecuencias 
K = 2 ;
f1 = 0.2;
delta = 0.5;
f2 = f1 + delta / N;
% Fases
phi1 = rand(1) * 2 * pi - pi;
phi2 = rand(1) * 2 * pi - pi; 
% Frecuencias angulares 
w1 = 2 * pi * f1 ;
w2 = 2 * pi * f2 ;
% Ruido 
SNR = 10 ;   
varE = 1 ;  % varianza/energia de UNA exponencial compleja
varV = K*varE * 10^(-SNR/10);   % SRN = 10*log10( K*varE / varV );
v = (randn (N,1) + 1i * randn (N,1))* sqrt(varV/2);
% Señal de entrada. 
k = ( 0:N-1 )';
y = exp(1i * (w1*k  + phi1)) + exp(1i * (w2*k+ phi2)) + v;

L_init = 3 ;    % L de MUSIC inicial y final 
L_fin  = 10; 

vec_w_est = zeros(L_fin - L_init, K );  % Guardo todos los w estimados y sus errores 
vec_err_w = zeros(L_fin - L_init, K );

for k=L_init:L_fin 
    
    clear aves avas R w_est zeros_est w1_est w2_est
    L = k ;
    R =  cov_FB(y, L);          % Estimo matriz de covarianza
    [aves , avas] = eig(R);     % Calculo autovalores
    [aves , avas] = sort_eig( aves, avas );  % Los ordeno.

    [ w_est, zeros_est] = music( y, L , K );
    w1_est = w_est(1);
    w2_est = w_est(2);

    e_w1 = abs( w1 - w1_est)/w1*100 ;
    e_w2 = abs( w2 - w2_est)/w2*100 ;
    
    vec_w_est(L-L_init+1,: ) = [ w1_est w2_est ];
    vec_err_w(L-L_init+1,: ) = [ e_w1   e_w2   ];

end

%
fprintf('Parametros Iniciales: \n SNR = %ddB\n',SNR );
fprintf(' Delta = %.2f \n w1 = %.3f \n w2 = %.3f \n\n', delta,w1,w2);
fprintf('Frecuencias estimadas:\n');
fprintf( 'Music L |');  fprintf( '\t%d \t|', L_init:L_fin );        fprintf('\n');
fprintf( 'Est w1  |');  fprintf( '    %.2f \t|', vec_w_est(:,1)');  fprintf('\n');
fprintf( 'Est w2  |');  fprintf( '    %.2f \t|', vec_w_est(:,2)');  fprintf('\n\n');
fprintf('Errores relativos de estimación:\n');
fprintf( 'Music L |');  fprintf( '\t%d \t|', L_init:L_fin );            fprintf('\n');
fprintf( 'Err w1  |');  fprintf( '    %.2f %% \t|', vec_err_w(:,1)');   fprintf('\n');
fprintf( 'Err w2  |');  fprintf( '    %.2f %% \t|', vec_err_w(:,2)'  ); fprintf('\n');
fprintf('\nConclusión: En algunas corridas esporadicas el valor de la estimación se dispara.\n')
fprintf('Pero en general para L>7 se logra estimar ambas frecuencias con un error menor al 10%%\n');



% Graficos de error de estimación en w 
figure, plot( L_init:L_fin , vec_err_w(:,1), 'r');
hold on, plot( L_init:L_fin , vec_err_w(:,2), 'b');
hold on, plot( [ L_init L_fin ], [ 10 10 ], '--k', 'linewidth', 2 );
title('Error de estimación en W vs. orden MUSIC');
xlabel('Orden L');
ylabel('Error porcentual');
legend('Error en w1','Error en w2');
ylim( [ 0 100 ]);
%print('-dpng','errorWvsL.png');



% Periodograma 
y_cor = xcorr(y,'biased');
Nfft = 1024 ; 
PSD = fftshift( abs(fft( y_cor, Nfft )) );
w = linspace( -1 , 1 ,Nfft);

figure
% % plot(w,10*log10( PSD ),'b')
% % ylabel('Amplitud [dB]')
plot(w,20*log10(PSD),'b')
hold on, plot( [w1 w1]/pi, ylim ,'--k');
hold on, plot( [w2 w2]/pi, ylim ,'--k');
ylabel('Amplitud')
title( sprintf( 'Periodograma. w1=%.2f w2=%.2f',w1/pi,w2/pi) );
grid on;
xlabel('Freq. [rad/\pi]')
%print('-dpng','periodogramaMUSIC.png');








