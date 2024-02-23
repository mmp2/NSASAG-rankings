%  Estimates the parameter and central permutation from infinite
%  Mallows models data, variable theta values

global Q0

%  Initializations

israndsigma = 0; % = 1 random central permutation
		 % = 0 central perm is identity
t = 32;
theta = 2*log(2).^(1:0.5:(t+1)/2)

nnn = [ 200 500 1000 2000 ];
nn = length( nnn );
niter = 50;
jtied = t;  % no tied parameters

theta0 = 0.1*ones(1,t);  % ini2=2, ini1=1, ini.1=0.1, ini0=theta

%  Trace

cost_all = zeros( t, niter, nn );
costb_all = zeros( t, niter, nn );
costc_all = zeros( t, niter, nn );
dsigb_all = zeros( niter, nn );
dsigbt_all = zeros( niter, nn );
dsigc_all = zeros( niter, nn );
dsigct_all = zeros( niter, nn );
theta_all = zeros( t, niter, nn );
nitems_all = zeros( niter, nn );
sigmaML_all = zeros( 600, niter, nn );
sigma_all = zeros( 600, niter, nn );

for in = 1:nn;
    
    nsamples = nnn( in )
    for iter = 1:niter;

%   Generate data

   [pp, s ] = sample_from_theta( theta, nsamples, t );
   ppsave = pp;
   items0 = unique( pp' );
   nitems = length( items0 )

   n0 = max( max( pp ));  % needed temporarily, will be reassigned to n

   if israndsigma
      sigma1 = randperm( 2*n0 );  % the central permutation
      sigma = sigma1( items0 );   % a relabeling 
      pp = sigma1( pp );
   else
      sigma = 1:n0;               % central perm is identity
   end;

   n = max( max( pp ));
   if ~exist( 'Q0', 'var' )         % this way, Q0 is the largest so far
      Q0 = triu( ones( n, n ), 1 );   
   else
      if n > size( Q0, 1 );
	 Q0 = triu( ones( n, n ), 1 );   
      end;
   end;

%  Estimate thetaML, sigb, sigc 

   est_vartheta

% Trace

  cost_all( :, iter, in) = cost';
  costb_all( :, iter, in) = costb';
  costc_all( :, iter, in) = costc';
  dsigb_all( iter, in) = dsigb;
  dsigc_all( iter, in) = dsigc;
  dsigbt_all( iter, in) = dsigbt;
  dsigct_all( iter, in) = dsigct;
  theta_all( :, iter, in) = thetaML';
  nitems_all( iter, in) = nitems;
  sigmaML_all( 1:nitems, iter, in ) = sigmaML';
  sigma_all( 1:n, iter, in ) = sigma';

%imagesc(Q0(pib,pib)-0.5*Q0(sigma,sigma))
%pause
end; % for iter
end; % for in 

%boxplot([cost_all, costb_all, costc_all, dsigb_all, dsigc_all, theta_all])

plot_vartheta




