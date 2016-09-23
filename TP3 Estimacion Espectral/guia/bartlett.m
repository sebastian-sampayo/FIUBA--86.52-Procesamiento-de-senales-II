function correl = bartlett( xin , M , N_puntos)

L_x = length( xin );

L = floor( L_x / M);

vecFFT = zeros(L,N_puntos);

for j = 1:L
    
    x = xin((j-1)*M+1 :j*M );

    vecFFT(j,:) = fftshift( abs( fft( x, N_puntos )).^2 )/M;
end;
    
   correl = mean( vecFFT , 1 );
   
end
    
    
    
    