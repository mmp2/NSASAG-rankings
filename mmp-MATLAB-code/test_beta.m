t = 5;
n = 10;
nech = 1;  % equivalent sample size
rmax = n*t-t*(t-1)/2;

xx = 1:0.1:rmax*nech;

bb = betaln( xx, nech*t+1 );

plot( bb );
nfact = gammaln( n )

