function [ sigma, thetaML, costR, cost, R, nsamt ] = ...
		  est_varthet_sdiscr(t, n, sigall, nsam, Q, q, jtied )
%function [ sigma, thetaML, costR, cost, R, nsamt ] = ...
%                 est_varthet_sdiscr(t, n, sigall, nsam, Q, q, jtied )
%
% Takes in a set of sigmas and computes the ML estimates of sigma
% theta restricted to those sigmas.
%
% t = maximum length of permutation
% n = max value in pp (change this...)
% jtied = number of free parameters. 
%	  theta( 1:jtied-1 ) are free, theta( jtied:t ) are tied
%	  all free <=> jtied = 0
% sigall( :, n ) = set of permutations, one per line
% nsam(t) = number samples for each rank j
% q(n,t) = sufficient statistics (normalized)
% Q(n,n,t) = sufficient statistics (normalized)
%
% sigma   = optimal permutation 
% thetaML( 1, t ) = theta estimate
% cost( 1, t ) = lsigma( R(:,:,j) )
% R = sum thetaML(j)*R(:,:,j)
% costR = lsigma( R )
% nsamt = total number samples for the tied parameters, if jtied ~= 0
%	= nsam if jtied == 0 or tmax

nsigmas = size( sigall, 1 );

Rj = zeros( n, n, t );
for jt = 1:t;
    Rj( :,:, jt )= q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt ));
end;

costall = ones( 1, nsigmas )*Inf;
thetall = zeros( nsigmas, t );

for is = 1:nsigmas;
    sigma = sigall( is, : );

   % Estimate theta
   for jt = 1:t;
       cost( jt ) = lsigma( sigma, squeeze( Rj(:,:,jt)) );
   end;
   
   if jtied == 0 || jtied == t
       thetaML = log( 1 + nsam./cost );
       nsamt = nsam;
   else
       thetaML( 1:jtied-1 ) = log( 1 + nsam( 1:jtied-1 )./cost( 1:jtied-1 ) );
       jjj = jtied:t;
       nsamt = sum( nsam( jjj ));
       costt = sum( cost( jjj ));
       thetaML( jjj ) = log( 1 + nsamt/costt );
   end;

   thetall( is, : ) = thetaML;
   costall( is ) = thetaML*cost' - log( 1-exp(-thetaML))*nsam';

end; % for

[ cmin, imin ] = min( costall );

sigma = sigall( imin, : );
thetaML = thetall( imin, : );

% Calculate R, cost 
   
R = zeros( n, n ); 
for jt = 1:t;
    R = R + thetaML( jt )*squeeze( Rj(:,:, jt ));
end;

costR = sum( sum( tril( R( sigma, sigma ), -1)));
for jt = 1:t;
    cost( jt ) = lsigma( sigma, squeeze( Rj(:,:,jt)) );  
end;



