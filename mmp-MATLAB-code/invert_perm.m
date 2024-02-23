function [ pinv ] = invert_perm( pp, n )

% function [ pinv ] = invert_perm( pp, n )
%
% Computes the inverse of permutation pp or rank n
% The items must be from 1:n
% GLOBAL Q0 (can be larger than n)

global Q0

P = eye( n ); P = P( pp, : );
Q = P'*Q0(1:n,1:n)*P;
pinv = 1+sum( Q, 1 );
if ~isempty( find( sum( P )-1 ))
   pinv = ones( 1, n );
end;
