% Samples from infinite Mallows model and generates sufficient statistics


% 	Initializations

%nsamples = 10;
%t = 4;
%theta = 1;

%       Generate permutations and sufficient statistics

[pp, s ] = sample_from_theta( theta*ones( 1, t ), nsamples, t );
items0 = unique( pp' );

n0 = max( max( pp ));
nitems = length( items0 );
sigma1 = 1:2*n0;
%sigma1 = randperm( 2*n0 );  % the central permutation
sigma = sigma1( items0 );

ppsave = pp;
pp = sigma1( pp );
n = max( max( pp ));

[ Q, q ] = make_Qp( pp, nsamples, t );

R = repmat( q, [1, n ])-Q';


