% pp = array of permutations, one/line
% nsamples = number samples (unused)
% t = maximum length of permutation
% theta0 = vector of initial values
% n = max value in pp (change this...)
% theta0MAX = defines prior on theta for BIC. Prior = 1/theta0MAX

theta0MAX = 20;

% Sufficient statistics

[ Q, q ] = make_Qmulti( pp, nsamples, t );  % Q is (n,n,t) array

nsam = sum( q, 1 );   % vector of nsamples in every positions

thetaML = theta0;
piold = zeros( 1, nitems );  % initial dummy values
pic = ones( 1, nitems );      %    "      "      "
pib = ones( 1, nitems );      %    "      "      "
costold = Inf;
iterest = 0;

while ((sum( abs( piold - pib )) > 0)...   % while the estimated permutation
     & (iterest < 30))	       % changes
			
   iterest = iterest+1;
   piold = pib;
   Rb = zeros( n, n );
   Rc = zeros( n, n );

   if iterest > 1
   for jt = 1:t;
       Rb = Rb + thetab( jt )*( q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt )));
       Rc = Rc + thetac( jt )*( q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt )));
   end;

   costac = sum( sum( tril( Rc( pic, pic ), -1)));
   costab = sum( sum( tril( Rb( pib, pib ), -1)));

%  choose which theta

   if costab < costac
      thetaML = thetab;
      R = Rb;
   else
%      disp('**est_vartheta** c chosen')
      thetaML = thetac;
      R = Rc;
   end;
   else      % first iteration
   R = zeros( n, n );
   for jt = 1:t;
       R = R + thetaML( jt )*( q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt )));
   end;
   end;

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

if exist( 'sigma', 'var' )
%  dsigb = dK( sigma, pib, n, nitems );
%  dsigc = dK( sigma, pic, n, nitems );
%  dsigbt = dK( sigma(1:t), pib( 1:t ), n, t );
%  dsigct = dK( sigma(1:t), pic( 1:t ), n, t );

  dsigb = dKdiff( sigma, pib, 0 );
  dsigc = dKdiff( sigma, pic, 0 );

  dsigbt = dKdiff( sigma(1:t), pib( 1:t ), 0 );
  dsigct = dKdiff( sigma(1:t), pic( 1:t ), 0 );
end;

  for jt = 1:t;
    Rj = q(:, jt )*ones( 1, n )-Q(:,:, jt );
    if exist( 'sigma', 'var' )
       cost( jt ) = sum(sum(tril( Rj( sigma, sigma ), -1)));
    end;

    costb( jt ) = sum( sum( tril( Rj( pib, pib ), -1)));
    costc( jt ) = sum( sum( tril( Rj( pic, pic ), -1)));
  end;

% Estimate theta

%  thetaML = log( 1 + nsam./costb )
%  thetaML
if (jtied == t) || (jtied == 0);
  thetab = log( 1 + nsam./costb );
  thetac = log( 1 + nsam./costc );
else
  thetab( 1:jtied-1 ) = log( 1 + nsam( 1:jtied-1 )./costb( 1:jtied-1 ) );
  thetac( 1:jtied-1 ) = log( 1 + nsam( 1:jtied-1 )./costc( 1:jtied-1 ) );
  jjj = jtied:t;
  nsamt = sum( nsam( jjj ));
  costbt = sum( costb( jjj ));
  costct = sum( costc( jjj ));
  thetab( jjj ) = log( 1 + nsamt./costbt );
  thetac( jjj ) = log( 1 + nsamt./costct );
end;

thetab( find( thetab==Inf )) = theta0MAX;
thetac( find( thetac==Inf )) = theta0MAX;

   

%pause
end; % while
iterest
if iterest == 1000
   disp( '**est_vartheta**  Max number iterations reached' );
end;

for jt = 1:t;
    Rb = Rb + thetab( jt )*( q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt )));
    Rc = Rc + thetac( jt )*( q(:, jt )*ones( 1, n )-squeeze( Q(:,:, jt )));
end;

costac = sum( sum( tril( Rc( pic, pic ), -1)));
costab = sum( sum( tril( Rb( pib, pib ), -1)));

%  choose which theta

if costab < costac
   thetaML = thetab;
   cost = costab;
   sigmaML = pib;
   if exist( 'sigma', 'var' )
      dsig = dsigb;
      dsigt = dsigbt;
   end;
else
   thetaML = thetac;
   cost = costac;
   sigmaML = pic;
   if exist( 'sigma', 'var' )
      dsig = dsigc;
      dsigt = dsigct;
   end;
end;

psi1 = 1 - exp( - thetaML);  % 1/psi(theta)
logl = -cost + nsam*log( psi1' );
if jtied == t % all  free params, each param has its onw sample size
   bic = logl - sum( log( nsam ) )/2;
else
   bic = logl - sum( log( [nsam( 1:jtied-1 ) nsamt]))/2;
end;

if jtied == t;
   bicbetter = bic + sum( log( exp( thetaML/2) - exp( -thetaML/2)));
   bicbetter = bicbetter + t*log(sqrt(2*pi)/theta0MAX);
else
   thdum = thetaML(1:jtied);
   bicbetter = bic + sum( log( exp( thdum/2) - exp( -thdum/2)));
   bicbetter = bicbetter + jtied*log(sqrt(2*pi)/theta0MAX);
end;

