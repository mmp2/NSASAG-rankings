% Generates partial rankings and tests sufficient statistics
%
% Permutation matrix: P(i,j) = 1 if j = pi( i )
% 
global Q0 Q0t

addpath ../B-and-B
%rand( 'state', 1871 );

n = 6;
t = n;  % observed ranks
ndata = 20;


P0 = eye( n );
Q0 = triu( ones( n, n ), 1 );
Q0t = Q0';
nn = 1:n;
nn0 = 0:n-1;

sigma1 = 1:n;   % true central permutation (inverse perm)
sigma1 = randperm( n );

Sig = P0( :, sigma1 );
sigma = sum( Sig*Q0', 2 )+1;

theta0 = ones( 1, n-1 );

% Generate data

[pp vv table ] = sample_from_theta( theta0, ndata );
size(pp);
pp = sigma1( pp );  

% Estimate distribution of V

empiv = zeros( n, n );  % empirical distribution of V

for j=1:n-1;
   tdum = tabulate( vv( :, j ));
   ldum = max( vv( :, j ));
   empiv( j, 1:ldum ) = tdum( :, 2 )'/ndata;
end;
empiv( n, 1 ) = 1;

ppp = pp( :, 1:t );    % partial samples

[ Q, R ] = make_Qp( ppp, ndata, n, t );


% Scoring

Sig2score = Sig;
score = nn0*Sig2score'*(1-R/ndata)+sum(sum(Q0'*Sig2score*Q*Sig2score'))

Sig2score = P0;
scoreP0 = nn0*Sig2score'*(1-R/ndata)+sum(sum(Q0'*Sig2score*Q*Sig2score'))

nrand = 50;
score_rand = zeros( 1, nrand );
for ii = 1:nrand;
    Sig2score = P0(:,randperm( n ));
    score_rand( ii ) = nn0*Sig2score'*(1-R/ndata)+sum(sum(Q0'*Sig2score*Q*Sig2score'));
end;

hp = plot( 1, score, 'ro', 2, scoreP0, 'kx', 1:nrand, score_rand, '.' );
title( 'scores' );

score_rand
sigma1

Q

%figure(1);
%imagesc(Q)
%pause
%plot(sort(R), 'r:o');
%title('sort(R)')

jkfkajads

Qalt = make_Q( pp );

figure(2);
imagesc( Qalt)
title( 'Qalt' );

