close all
clear

%% Punto 2.1
% Estimadores de la covarianza.

x = (randn (10000,1) + 1i * randn (10000,1)) / sqrt(2); %Ruido blanco circ. simet. con varianza unitaria

sesg    = xcorr(x,'biased'); %estimador de la covarianza con estimador sesgado
insesg  = xcorr(x,'unbiased'); %estimador de la covarianza con estimador insesgado

%Graficos 
figure, plot(abs(insesg))
hold,   plot(abs(sesg),'-r')

title('Estimadores de la covarianza')
legend('Estimador insesgado','Estimador sesgado')
% print('-dpng', 'estimadoresCovarianza.png');


%% Punto 2.2
% Estimador de Blackman-Tukey
%
% usar funcion blackmanTukey2 
%

B = [ 1     -1.3817 1.5632  -0.8843 0.4096  ];
A = [ 1     0.3544  0.3508  0.1736  0.2401  ];   

N_w = 256 ;      % Cantidad de muestras a tomar del periodograma
w = linspace(-pi, pi, N_w );
J = 100 ;        % Iteraciones
N = 64 ;         % Largo de la señal de entrada

% Calculo respuesta real
% [k,v] = tf2latc( B,A );
% Hd = dfilt.latticearma(k,v);
% [ hr, wr ] = freqz( Hd ,N_w );
[hr, wr] = freqz(B, A);

% Periodograma con N = 64
vecPSD_64 = zeros( J , N_w );     % Inicialización. Acá se van a guardar los J vectores de PSD

for j=1:J
    x = (randn (N,1) + 1i * randn (N,1)) / sqrt(2);     % Ruido blanco circ. simet. con varianza unitaria
    y = filter(B , A , x );                             % Filtrdo con H = B/A ; 

     vecPSD_64(j,:) = real( blackmanTukey2( y , w , 'rectangle', 'biased') );
%     y_cov = xcorr(y, 'biased');
%     Nfft = N_w;
%     vecPSD_64(j,:) = abs( fft( y_cov, Nfft) )/Nfft;
%     [aux, ~] = periodogram(y, [], N_w);
%     vecPSD_64(j,:) = aux';
end
PSD_64 = mean( vecPSD_64 , 1 ) ;      % La media de todos los PSD calculados
sigma_PSD_64 = sqrt( sum( ( vecPSD_64 - repmat( PSD_64 , J , 1) ).^2 )/(J-1) );
% sigma_PSD_64 = var( vecPSD_64 ,1 );

% Periodograma con N = 512
N = 512 ;         % Largo de la señal de entrada
vecPSD_512 = zeros( J , N_w );     % Inicialización. Acá se van a guardar los J vectores de PSD

for j=1:J
    x = (randn (N,1) + 1i * randn (N,1)) / sqrt(2);     % Ruido blanco circ. simet. con varianza unitaria
    y = filter(B , A , x );                             % Filtrdo con H = B/A ; 

    vecPSD_512(j,:) = real( blackmanTukey2( y , w , 'rectangle', 'biased') );
end
PSD_512 = mean( vecPSD_512 , 1 ) ;      % La media de todos los PSD calculados 
sigma_PSD_512 = sqrt( sum( (vecPSD_512 - repmat( PSD_512 , J , 1)).^2 )/(J-1) );

%% Grafico

figure, plot( wr, 20*log10( abs( hr ) ),'g', 'linewidth', 2 )
title('PSD real')
xlabel('Freq.');
ylabel('Amplitud [dB]');
xlim([0 pi]);

%print('-dpng','PSDreal.png');

figure, plot( wr, 20*log10( abs( hr ) ),'g', 'linewidth', 2 )
w = linspace(-pi, pi, N_w );
hold on,plot(w, 10*log10( (PSD_64)), 'r', 'linewidth', 2);
hold on, plot(w, 10*log10( (PSD_512)), 'b' ,'linewidth', 2 );
% margenes de sigma 
hold on,plot(w, 10*log10( abs(PSD_64 - sigma_PSD_64 ) ), '--r');
hold on,plot(w, 10*log10( abs(PSD_64 + sigma_PSD_64 ) ), '--r');
hold on,plot(w, 10*log10( abs(PSD_512 - sigma_PSD_512 ) ), '--b');
hold on,plot(w, 10*log10( abs(PSD_512 + sigma_PSD_512 ) ), '--b');
xlim([0 pi]);

title('PSD. Periodograma por Blackman-Tukey')
xlabel('Freq.');
ylabel('Amplitud [dB]');
legend('Real', 'Estimada. N = 64 ', 'Estimada. N = 512 ', ... 
        'Margen N = 64','Margen N = 64','Margen N = 512','Margen N = 512','Location','best' )
%print('-dpng', 'periodogramaBlackmanTukey.png');

%% Punto 2.2 c) 
% Estimador de Barlet 

% Periodograma con N = 256
J = 100 ;
N = 256 ;         % Largo de la señal de entrada
N_w = 256 ;
w = linspace(-pi, pi, N_w );
Nfft = 2^14 ;

vecPSD_BT256 = zeros( J , N_w );     % Inicialización. Acá se van a guardar los J vectores de PSD
vecPSD_B_4   = zeros( J , Nfft);
vecPSD_B_16  = zeros( J , Nfft);

for j=1:J
    x = (randn (N,1) + 1i * randn (N,1)) / sqrt(2);     % Ruido blanco circ. simet. con varianza unitaria
    y = filter(B , A , x );                             % Filtrdo con H = B/A ; 

    vecPSD_BT256(j,:) = real( blackmanTukey2( y , w , 'triangle', 'biased') );
    vecPSD_B_4(j,:) = bartlett( y , N/4 , Nfft );
    vecPSD_B_16(j,:) = bartlett( y , N/16 , Nfft );
    
end;

% Medias y varianzas.

% Normalización, no necesaria.
% % % PSD_BT256 = mean( vecPSD_BT256, 1 );
% % % ma = max( abs( PSD_BT256 ) );
% % % PSD_BT256 = PSD_BT256/ma ;
% % % sigma_PSD_BT256 = sqrt( var( vecPSD_BT256, 1 ) )/ma;
% % % 
% % % PSD_B_4 = mean( vecPSD_B_4, 1 );
% % % ma = max( abs( PSD_B_4) );
% % % PSD_B_4 = PSD_B_4 / ma ;
% % % sigma_PSD_B_4 = sqrt(var( vecPSD_B_4, 1 ) )/ma;
% % % 
% % % PSD_B_16 = mean( vecPSD_B_16, 1 );
% % % ma = max( abs( PSD_B_16) );
% % % PSD_B_16 = PSD_B_16 / ma ;
% % % sigma_PSD_B_16 = sqrt(var( vecPSD_B_16, 1 ))/ma;


PSD_BT256 = mean( vecPSD_BT256, 1 );
sigma_PSD_BT256 = sqrt( var( vecPSD_BT256, 1 ) );

PSD_B_4 = mean( vecPSD_B_4, 1 );
sigma_PSD_B_4 = sqrt(var( vecPSD_B_4, 1 ) );

PSD_B_16 = mean( vecPSD_B_16, 1 );
sigma_PSD_B_16 = sqrt(var( vecPSD_B_16, 1 ));


%% Graficos 

w = linspace(-pi,pi,length(sigma_PSD_BT256));
figure, plot( wr, 20*log10( abs( hr )),'g', 'linewidth', 2  )
hold on, plot( w, 10*log10(PSD_BT256), 'r', 'linewidth', 2 );

w = linspace(-pi,pi,Nfft);
hold on, plot( w, 10*log10(PSD_B_4) ,'b', 'linewidth', 2 );
hold on, plot( w, 10*log10(PSD_B_16),'m', 'linewidth', 2 );

% Margenes 
w = linspace(-pi,pi,length(sigma_PSD_BT256));
hold on, plot( w, 10*log10(PSD_BT256-sigma_PSD_BT256), '--r' );
hold on, plot( w, 10*log10(PSD_BT256+sigma_PSD_BT256), '--r' );

w = linspace(-pi,pi,Nfft);
hold on, plot( w, 10*log10(PSD_B_4-sigma_PSD_B_4) ,'--b' );
hold on, plot( w, 10*log10(PSD_B_4+sigma_PSD_B_4) ,'--b' );
hold on, plot( w, 10*log10(PSD_B_16-sigma_PSD_B_16),'--m');
hold on, plot( w, 10*log10(PSD_B_16+sigma_PSD_B_16),'--m');

% Propiedades 
xlim([0 pi]);
title('PSD. Periodograma')
xlabel('Freq.');
ylabel('Amplitud [dB]');
legend( 'Real', 'Blackman-Tukey', 'Bartlet. M = N/4', 'Bartlet. M = N/16',...
        'Margen Blackman-Tukey', 'Margen Blackman-Tukey',...
        'Margen Bartlet M = N/4', 'Margen Bartlet. M = N/4',...
        'Margen Bartlet. M = N/16', 'Margen Bartlet. M = N/16','Location','NorthWest' ); 

% print('-dpng', 'BTyBartlett.png')







