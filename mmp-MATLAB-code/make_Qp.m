function [ Q, q ] = make_Qp( ppp, nsamples, t )
% function [ Q, q ] = make_Qp( ppp, nsamples, t )
%
% 	 constructs the Q, q matrices from a set ppp of
%	 partial permutations of size t.
%        ppp( nsamples, t )
%	 (n : the total number of items n is max(max(ppp), so that
% 	 this function works for top-t rankings of infinite items)
%	 t = number of observed ranks 1:t, (t<n if n finite)
%	 Q( i, j ) = number times i < j (i before j) in ppp
%		     n x n matrix
%        q( i ) = number of times i present, column vector
% 

n = max( max( ppp ));

% 	Compute q 

q = zeros( n, 1 );
for jj = 1:t;
    tdum = tabulate( ppp(:,jj) );
    q( tdum(:, 1 )) = q( tdum(:, 1 )) + tdum( :, 2 );
end;

% 	Compute Q

% destroys ppp inside this loop
% ppsave = ppp;

Q = zeros( n, n );
for jj = 1:t-1;
    ppp = sortrows( ppp );
    tdum = tabulate( ppp( :, 1 ));
    idum = tdum( :, 1 );
    tdum = tdum( :, 2 );
    cdum = [0; cumsum( tdum )];
    for jj2 = 2:t-jj+1;
	for ii = 1:length(idum);
            i = idum( ii );
	    tt = tabulate( ppp( (cdum( ii )+1):cdum(ii+1), jj2 ) );
	    if ~isempty( tt );
	       Q( i, tt(:,1)) = Q( i, tt(:,1))+tt(:, 2 )';
	    end;
%	    imagesc( Q , [0 nsamples ])
%	    title( ['columns ' num2str( jj ) '  ' num2str( jj2 )]);
%	    pause;
        end;
    end;
    ppp = ppp( :,2:end );
end;

%   Nice idea but doesn't work here:
%       Q = Q + Q0t( ppi, ppi );  % wrong, should be Q0(pi, pi)
%				  % not the inverse but the direct!!
 
