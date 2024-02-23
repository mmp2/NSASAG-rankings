function s = lsigma( sigma, R )

% function s = lsigma( sigma, R )
%
% L_sigma( R ) where R is a matrix and sigma a vector of indices in 
% 	       the rows and colums of R

s = sum(sum(tril(R(sigma, sigma),-1)));

