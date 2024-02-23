% NOT FINISHED (2012)
% Computes the marginal probabilities P[ pi(j) = k ] under the FINITE
% Mallows model.
% see ../B-and-B/compute_marginal.m

% Initializations

n1 = n + m;

theta = 0.5*ones( 1, n-1 );

%  Compute P[ Vj ] or P[ Sj ] (from B-and-B/sample_from_theta.m)

vtable = exp(-[theta 1]'*(0:(n-1)));  % Probability table for the V's
vtable = vtable.*hankel(ones(n,1));
vtable = vtable./repmat(sum(vtable, 2), [ 1, n ]);

cdf = cumsum( vtable' )';  % filled with 1's in irrelevant positions

%  Construct F matrices 

FFF = cell( 1, n-1 );  % stores the F matrices
		       % cell j is an (n-j+1) square matrix whose last column=0
		       % therefore ==>  cell j is (n-j+1)x(n-j)
		       % (note: column sums of FF{j} = 1 
for j = 1:n-1;
    FFF{ j } = diag( 1 - cdf( j, 1:n-j+1 ))+diag( cdf( j, 1:n-j), -1 );
    FFF{ j } = FFF{ j }(:, 1:n-j );
end;

%  Recursion to compute the probabilities

Pmarginal = cell( 1, n );    % stores the current marginals, 
			     % computed recursively backwards from n to 1
			     % at stage j cells j:n are not empty
			     % each stores a n-j+1 marginal column vector

Qmarginal = zeros( n, n );   % stores the expectation of Qij
Pmarginal{n} = [ 1 ];

for j = n-1:-1:1;
    Pmarginal{ j } = vtable( j, 1:n-j+1 )';
    for k = j+1:n;
        Qmarginal( j, k ) = cdf( j, 1:n-j )*Pmarginal{ k };
	Pmarginal{ k } = FFF{ j }*Pmarginal{ k }; 
    end;
end;

pmarg = zeros( n, n );
for j = 1:n;
    pmarg( :, j ) = Pmarginal{ j };
end;

Qmarginal = Qmarginal + tril( ones( n, n ), -1 ) - Qmarginal';

% Interesting plots

if 0
clf
plot( 1:n, diag( pmarg ), 'r.', 1:n-1, diag( pmarg, 1 ), 'b.', 1:n-1,diag(pmarg, -1 ), 'go', 1:n-2, diag(pmarg, -2 ), 'yo', 1:n-2, diag( pmarg, 2 ), 'c.', 1:n-3, diag(pmarg, 3 ), 'm.')
legend( 'diag', '+1', '-1', '-2', '+2', '+3', 'North' )
title( 'P marginal' );

pause

plot( 1:n, diag( Qmarginal ), 'r.', 1:n-1, diag( Qmarginal, 1 ), 'b.', 1:n-1,diag(Qmarginal, -1 ), 'go', 1:n-2, diag(Qmarginal, -2 ), 'yo', 1:n-2, diag( Qmarginal, 2 ), 'c.', 1:n-3, diag(Qmarginal, 3 ), 'm.')
legend( 'diag', '+1', '-1', '-2', '+2', '+3', 'North' )
title( 'Q' );

pause
end; % if 0

if 0
% check the asymptotic formula

n2 = round( n/2 );
pasympt = pmarg( :, n2)';

ii = 1:9;
for j = 1:8;
%    ptest = [ sum( pasympt( 1:(n2-j))) pasympt( n2-j+1:end ) ];
%    not so good approximation: first value too large, next ones too small
%     ptest =  [pasympt( (n2-j):-1:1)) pasympt( n2-j+2:end ) ];
    semilogy( ii, pmarg( j, ii ), 'r.', ii, ptest( ii ), 'bo' );
    title( num2str( j ))
    pause
end;
end; % if 0

%  Verification

if 0

nsample = 5000;   % Verified for nsample up to 50000
pp = sample_from_theta( theta, nsample );
Qemp = make_Q( pp )/nsample;

pmarg_emp = zeros( n, n );
for j = 1:n;
    dum = ( pp == j );
    pmarg_emp( :, j ) = sum( dum, 1 )'/nsample;
end;

max( max( abs( pmarg-pmarg_emp)))

subplot( 121 )
imagesc( pmarg );
title( 'pmarg - computed' );

subplot( 122 )
imagesc( pmarg-pmarg_emp );
title( 'pmarg - pmarg empirical' ); colorbar;

pause

max( max( abs( Qmarginal-Qemp)))

subplot( 121 )
imagesc( Qmarginal );
title( 'Qmarginal - computed' );

subplot( 122 )
imagesc( Qmarginal-Qemp );
title( 'Qmarginal - Qemp' ); colorbar;


end; % if 0


% make pretty figure

clf

imagesc( Qmarginal, [0 1 ] );
set( gca, 'XTick', []);
set( gca, 'YTick', []);
%set( gca, 'FontSize', fontsize );

%print( '-dpsc', [ 'fig-Qmarginal-theta' num2str( theta( 1 ), 1 ) '.eps' ] );

pause

imagesc( pmarg ); colorbar;
%imagesc( pmarg, [0 1 ] );
set( gca, 'XTick', []);
set( gca, 'YTick', []);

%print( '-dpsc', [ 'fig-Pmarginal-theta' num2str( theta( 1 ), 1 ) '.eps' ] );

%print( '-dpsc', [ 'fig-Pmarginal-theta' num2str( theta( 1 ), 1 ) '-unnormcolor.eps' ] );

if 1   % doesn't work! neither of two formulas agrees with Q
pause

%   Compute pi_ij by Mallows' formula  --> see my notes 
%
%   pi_ij = depends only on k=j-i > 0
%   phi = exp( -theta ) 
%   only for constant theta

phi = exp( theta( 1 ));   % assumes theta constant
%Qmal = (phi/(1+phi))^2*ones( 1, n-1 );
Qmal = zeros( 1, n-1 );
for k = 2:n;
   %Qmal( k-1) = Qmal( k-1 )*(phi^(2*k-1)-k*phi+(k-1)/phi)/(1-phi^k)/(1-phi^(k-1));
   Qmal( k-1) = k*(phi^(2*k+2)+1)/(phi^(2*k+2)-1)-(k-1)*(phi^(2*k)+1)/(phi^(2*k)-1); % this is from (19) =2p_k-1
end;
Qmal = (Qmal+1)/2; % this is from (19) =2p_k-1

Qmal
Qk = Qmarginal(1,2:end);
Qk

end;
