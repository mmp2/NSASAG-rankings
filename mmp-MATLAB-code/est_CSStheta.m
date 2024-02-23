% estimation of sigma and theta by the CSS method
% used by compare_cohen_univ_css.m


% Construct the CSS Q matrix

[ Q, q ] = make_Qmulti( pp, nsamples, t );

Q = sum( Q, 3 );
q = sum( q, 2 );

R = Q'-Q+nsamples/2*ones( n, n );   % CSS method

%R = q*ones( 1, n )-Q;   % constant theta method

				    
%imagesc(R)
%title ('R')

%pause

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

  pib;
  pic;

  costb = sum( sum( tril( R( pib, pib ), -1)))/nsamples
  costc = sum( sum( tril( R( pic, pic ), -1)))/nsamples

% Estimate theta

  if 1  % always use greedy
%  if costc < costb 
     thetaCSS = log( 1 + t/costc );
     sigmaCSS = pic;
   else
     thetaCSS = log( 1 + t/costb );
     sigmaCSS = pib;
  end;

