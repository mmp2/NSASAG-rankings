function [ Q, q ] = make_Qmulti( ppp, nsamples, t )
% function [ Q, q ] = make_Qmulti( ppp, nsamples, t )
%
% 	 constructs the Q_j, q_j matrices from a set ppp of
%	 partial permutations of size t (or smaller)
%
%        ppp( nsamples, t )
%	 (n : the total number of items n is max(max(ppp), so that
% 	 this function works for top-t rankings of infinite items)
%        nsamples: unused!!
%	 t = number of observed ranks 1:t, (t<n if n finite)
%	     if t variable, then t = dim( ppp, 2 ) and the unobserved
%	     ranks are padded with -1
%	 Q(i, k, j ) = number times i < k (i before k) 
% 		       with rank(i)=j, rank(k)<j 
%		     n x n x t array
%        q( i, j ) = number of times rank(i) = j
%		     n x t matrix		
% 

n = max( max( ppp ));

% 	Compute q 

q = zeros( n, t );
for jj = 1:t;
    tdum = tabulate( ppp(:,jj) );
    if( tdum( 1,1 ) < 1 )  % fill-in-values for unequal lengths
	q( tdum(2:end, 1 ), jj) = tdum( 2:end, 2 );
    else
        q( tdum(:, 1 ), jj) = tdum( :, 2 );
    end;
end;

% 	Compute Q

% destroys ppp inside this loop
% ppsave = ppp;

ppp = ppp(:, t:-1:1);

Q = zeros( n, n, t );
for jj = 1:t-1;
    ppp = sortrows( ppp );
    tdum = tabulate( ppp( :, 1 ));
    idum = tdum( :, 1 );
    tdum = tdum( :, 2 );
    cdum = [0; cumsum( tdum )];
    for jj2 = 2:t-jj+1;
	for ii = 1:length(idum);
            i = idum( ii );
	    if i >= 1   % otherwise it's a -1 value to be eliminated
	       tt = tabulate( ppp( (cdum( ii )+1):cdum(ii+1), jj2 ) );
	       if ~isempty( tt );
                  if tt( 1,1 ) < 1
		     tt = tt( 2:end, : );
                  end;
		  Q( i, tt(:,1), t-jj+1 ) = Q( i, tt(:,1), t-jj+1 ) + tt(:, 2 )';
	       end;
            end;
        end;

    end;
    ppp = ppp( :,2:end );
end;

%   Nice idea but doesn't work here:
%       Q = Q + Q0t( ppi, ppi );  % wrong, should be Q0(pi, pi)
%				  % not the inverse but the direct!!
 
