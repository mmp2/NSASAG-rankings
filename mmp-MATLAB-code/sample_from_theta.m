function [ pp s ] = sample_from_theta( theta, nsamples, t )

% function [ pp s ] = sample_from_theta( theta, nsamples, t )
%
% SAMPLE_FROM_THETA Takes nsamples samples from the infinite model 
% dictated by the given theta, with central permutation = identity
%
% t		     = the length of the top-t permutation
% theta		     = vector of dimension t with >0 values
% nsamples 	     = number of samples	
% pp( nsamples, t )  = sample permutations, one per row
% s( nsamples, t ) = same permutations in the s representation
%		       s( i, j ) = s_j+1 in permutation i

% Can we obtain the Q directly from the V's ?

s = zeros(nsamples, t);  % The permutations table

for j=1:t % For every entry in this vector
    % sample s_j from discrete exponential distribution with parameter theta(j)
    s( :, j ) = floor( random( 'exp', 1/theta( j ), nsamples, 1 ));  
end;

% Construct the permutation based on the s's. The sample was of the form
% representing the number of higher values coming before; this constructs
% the actual permutation in the standard vector representation

pp = zeros( nsamples, t );
iii = repmat( 1:t, [ nsamples, 1 ]);

s = s + 1;
for j = 1:t;    
    sj = s( :, j );
    ssort = sort( pp(:, 1:j-1 ), 2 );
    for j1 = 1:j-1;
	sj = sj + double( ssort( :, j1 ) <= sj );  % can be done faster
    end;
    pp( :, j ) = sj;
end;


