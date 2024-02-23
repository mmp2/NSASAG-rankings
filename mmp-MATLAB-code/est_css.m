% pp = array of permutations, one/line
% nsamples = number samples (unused)
% t = maximum length of permutation
% theta0 = vector of initial values
% n = max value in pp (change this...)

% Sufficient statistics

[ Q, q ] = make_Qmulti( pp, nsamples, t );  % Q is (n,n,t) array

nsam = sum( q, 1 );   % vector of nsamples in every positions

% nsam( 1 ) = total number of samples (number of experts that have data)
% Q(i,j)+Q(j,i) = total number of samples for a pair. 
% to add 1/2*(nsam( 1 )-Q(i,j)-Q(j,i)) 

Q = sum( Q, 3 );
Q = (Q - Q'+nsam( 1 ))/2;

thetaML = theta0;
piold = zeros( 1, nitems );  % initial dummy values
pic = ones( 1, nitems );      %    "      "      "
pib = ones( 1, nitems );      %    "      "      "
costold = Inf;
iterest = 0;

% CSS (Cohen, Schap, Singer) greedy

  vb = sum( Q, 1 );
  pic = zeros( 1, nitems );
  for ii=1:nitems;
     [ dum, jj ] = min( vb );
     imin = find( vb == dum );
     if imin( end )~= jj;           % checks that there are multiple minima
	[ dum, jdum ] = max( q( imin ));  % if yes, takes the one with the
	jj = imin( jdum );		  % largest q
     end;
     pic( ii ) = jj;
     vb = vb - Q( jj,: );
     vb( jj ) = inf;
  end;

if exist( 'sigma', 'var' )
  dsigb = dK( sigma, pib, n, nitems);
  dsigc = dK( sigma, pic, n, nitems );

  dsigbt = dK( sigma(1:t), pib( 1:t ), n, t );
  dsigct = dK( sigma(1:t), pic( 1:t ), n, t );
end;


% Estimate theta

  thetaML = log( 1 + nsam./costb )

  thetab = log( 1 + nsam./costb );
  thetac = log( 1 + nsam./costc );



