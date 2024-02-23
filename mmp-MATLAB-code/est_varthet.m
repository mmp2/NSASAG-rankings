function [ sigma, thetaML, costR, iterest, R, nsamt ] = est_varthet(t, n, ...
                                                  nsam, Q, q,thetaini, maxiter, jtied )
%function [ sigma, thetaML, costR, iterest, R, nsamt ] = est_varthet(t, n,
%nsam, Q, q, thetaini, maxiter, jtied)
%
% Function version of ../infinite/est_vartheta.m
%
% t = maximum length of permutation
% thetaini = vector of initial values (default 0.1)
% n = max value in pp (change this...)
% nsam(t) = number samples for each rank j
% theta0MAX = (defines prior on theta for BIC. Prior = 1/theta0MAX (later!))
%           = upper bound against thetaML=Inf
% maxiter = max number of iterations (default = 50)
% q(n,t) = sufficient statistics (normalized)
% Q(n,n,t) = sufficient statistics (normalized)
%
% sigma   = optimal permutation <-- found by est_sigma_heur
% thetaML( 1, t ) = theta estimate
% cost( 1, t ) = lsigma( R(:,:,j) )
% iterest = number iterations to convergence
% R = sum thetaML(j)*R(:,:,j)
% costR = lsigma( R ) (not returned)
%
% can work with tied parameters!! (implemented but not used)

thetaML = thetaini;
sigold = zeros( 1, n );  % initial dummy values
sigma = 1:n;

%jtied = 0;  % non-zero if tied parameters
theta0MAX = 10; % upper bound on theta

Rj = zeros( n, n, t );
for jt = 1:t;
    Rj( :,:, jt )= q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt ));
end;

iterest = 0;

while ((sum( abs( sigold - sigma )) > 0)...   % while the estimated permutation
     & (iterest < maxiter))	       % changes
			
   iterest = iterest+1;
   sigold = sigma;

   % Calculate new R
   
   R = zeros( n, n ); 
   for jt = 1:t;
       R = R + thetaML( jt )*squeeze( Rj(:,:, jt ));
   end;

   %  Estimate sigma

   [sigma, costR] = est_sigma_heur( R, q*thetaML', int8( nsam>1 ) );
   % don't search if nsam == 1 
   
   % Estimate theta
   for jt = 1:t;
       cost( jt ) = lsigma( sigma, squeeze( Rj(:,:,jt)) );
   end;
   
   if jtied == 0
       thetaML = log( 1 + nsam./cost );
       nsamt = nsam;
   else
       thetaML( 1:jtied-1 ) = log( 1 + nsam( 1:jtied-1 )./cost( 1:jtied-1 ) );
       jjj = jtied:t;
       nsamt = sum( nsam( jjj ));
       costt = sum( cost( jjj ));
       thetaML( jjj ) = log( 1 + nsamt/costt );
   end;
   thetaML = min( thetaML, theta0MAX );
end; % while
iterest

% Calculate final R, cost

R = zeros( n, n );
for jt = 1:t;
    R = R + thetaML( jt )*squeeze( Rj(:,:, jt ));
end;
costR = sum( sum( tril( R( sigma, sigma ), -1)));
for jt = 1:t;
    cost( jt ) = lsigma( sigma, squeeze( Rj(:,:,jt)) );  % squeeze
end;



