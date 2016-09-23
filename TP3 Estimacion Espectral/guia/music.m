function [w , z] = music (y, L, K)

    if( L-K < 1 )
        return ;
    end;

    N = length( y );
    R =  cov_FB(y, L);   % Estimo matriz de covarianza
    
    [Vi ,Di] = eig(R);    % Calculo autovalores
    [V , D_avas] = sort_eig( Vi, Di );  % Los ordeno.
   
    G = V(:,1:L-K);         % Me quedo con los aves correspondientes a los menores avas 
    
    A = G * G';             % Armo mi matriz A
    
    poli = zeros( 2*L-1, 1);
    for i= -(L-1):L-1               % Calculo los coeficientes del polinomio.
        poli(i+L) = sum(diag(A,i));
    end

    zeros_p = roots(poli);                      % Busco los ceros del polinomio. 
    zeros_int = zeros_p( abs(zeros_p)<=1 );     % Me quedo con los del interior del circulo unidad.
    
    zeros_int = sort( zeros_int) ; 
    
    zeros_est = zeros_int(end-K+1:end) ; % Me quedo con los K mas proximos al circulo unidad 

    w = sort( angle(zeros_est) );
    z = zeros_est ;