%  Estimates the parameter and central permutation from infinite
%  Mallows models data

global Q0

%  Initializations

%nsamples = 1000;
t = 2;
theta = 2*log(2);

nnn = [ 200 500 1000 2000 ];
nn = length( nnn );
niter = 50;


%  Trace

cost_all = zeros( niter, nn );
costb_all = zeros( niter, nn );
costc_all = zeros( niter, nn );
dsigb_all = zeros( niter, nn );
dsigbt_all = zeros( niter, nn );
dsigc_all = zeros( niter, nn );
theta_all = zeros( niter, nn );
thetat_all = zeros( niter, nn );
nitems_all = zeros( niter, nn );

for in = 1:nn;
    
    nsamples = nnn( in )
    for iter = 1:niter;

%   Generate data

[pp, s ] = sample_from_theta( theta*ones( 1, t ), nsamples, t );
items0 = unique( pp' );

n0 = max( max( pp ));
nitems = length( items0 );
%sigma1 = 1:2*n0;
sigma1 = randperm( 2*n0 );  % the central permutation
sigma = sigma1( items0 );

ppsave = pp;
pp = sigma1( pp );
n = max( max( pp ));

[ Q, q ] = make_Qp( pp, nsamples, t );

R = repmat( q, [1, n ])-Q';


%    nitems,n0

    Q0 = triu( ones( n0, n0 ), 1 );


%   Estimate sigma

% FV average rank

  vb = sum( R, 2 );
  [ dum, pib ] = sort( vb', 2, 'descend' );
  pib = pib( 1:nitems );

% CSS (Cohen, Schap, Singer) greedy

  vb = sum( R, 1 );
  pic = zeros( 1, nitems );
  for ii=1:nitems;
     [ dum, jj ] = min( vb );
     imin = find( vb == dum );
     if imin( end )~= jj;           % checks that there are multiple minima
	[ dum, jdum ] = max( q( imin ));  % if yes, takes the one with the
	jj = imin( jdum );		  % largest q
     end;
     pic( ii ) = jj;
     vb = vb - R( jj,: );
     vb( jj ) = inf;
  end;


% Compute costs

  sigma;
  pib;
  pic;
  dsigb = dk( sigma, pib, n, nitems );
  dsigc = dk( sigma, pic, n, nitems );

  dsigbt = dk( sigma(1:t), pib( 1:t ), n, t );

  cost = sum(sum(tril( R( sigma, sigma ), -1)))/nsamples;

  costb = sum( sum( tril( R( pib, pib ), -1)))/nsamples;
  costc = sum( sum( tril( R( pic, pic ), -1)))/nsamples;


% Estimate theta

  thetaML = log( 1 + t/costc );

  trunc = 1;
  % truncated theta estimation
  costbt = sum( sum( tril( R( pib(1:trunc), pib(1:trunc) ), -1)));
  nsamt = sum(q(pib(1:trunc))) - sum( sum( Q( pib( trunc+1:end ), pib( 1:trunc ))));
	% Q is the transpose!!
	% adds all the times the first t items in pib appear, + all the
	% other items that appear in R in their rows

  thetaMLt = log( 1 + nsamt/costbt );

% Trace

  cost_all( iter, in) = cost;
  costb_all( iter, in) = costb;
  costc_all( iter, in) = costc;
  dsigb_all( iter, in) = dsigb;
  dsigbt_all( iter, in) = dsigbt;
  dsigc_all( iter, in) = dsigc;
  theta_all( iter, in) = thetaML;
  thetat_all( iter, in) = thetaMLt;
  nitems_all( iter, in) = nitems;

%imagesc(Q0(pib,pib)-0.5*Q0(sigma,sigma))
%pause
end; % for iter
end; % for in 

%boxplot([cost_all, costb_all, costc_all, dsigb_all, dsigc_all, theta_all])

plot_estimation_nips




