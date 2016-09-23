function [ V , D ] = sort_eig( Vi , Di )


V = Vi ;
Di = diag( Di ) ;
D = Di ;

N = length( D );

for i = 1:N 
    
    [ mi, pos] = min( Di );
    D(i) = mi ;             % Guardo los minimos en mi nuevo vector y matriz
    V(:,i)  = Vi(:,pos);
    
    Di(pos) = [];           % Borro la posicion del minimo.
    Vi(:,pos) = [] ;
       
end;

end 
