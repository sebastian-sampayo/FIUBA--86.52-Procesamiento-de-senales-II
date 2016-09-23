function [ correlograma ] = blackmanTukey2( xin , wi , win , estimador)

L = length( xin );

if( strcmp( win, 'triangle') )
    W = triang(2*L - 1) ;
else
    W = ones(2*L - 1, 1);
end 

r = xcorr(xin,estimador);   % Estimo la autocorrelación

k = -(L-1):(L-1);
k = k' ;

funCorr = @( w )(  sum( r.*W.*exp(-1i*w*k) ) ) ;

correlograma = zeros(1, length(wi) );   % inicio el vector

for j=1:length(wi)
    correlograma(j) = funCorr( wi(j) ); % calculo el valor de la funcion para los valores que me pide el usuario
end

end
