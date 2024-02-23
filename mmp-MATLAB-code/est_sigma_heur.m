function [ sigma, cost ] = est_sigma_heur( R, q, islocal_search )

% [ sigma, cost ] = est_sigma_heur( R, q, islocal_search )
%
% Estimates the optimal sigma by heuristics FV and greedy with
% local search (if islocal_search==1).
% sigma = argmin( tril( R( sigma, sigma )))

    global ncand ienq picand costcand isfullcand islocalmin LMAX
        
% 1. Find initial solutions
    
    % FV average rank

    n = size( R, 1 );
    vb = sum( R, 2 );
    [ dum, pib ] = sort( vb', 2, 'descend' );
    pib = pib( 1:n );

    % CSS (Cohen, Schap, Singer) greedy

    vb = sum( R, 1 );
    pic = zeros( 1, n );
    for ii=1:n;
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

    costb = sum( sum( tril( R( pib, pib ), -1)));
    costc = sum( sum( tril( R( pic, pic ), -1)));

    if islocal_search    
        %2. Local search
        % We define a queue in which to store solutions and their
        % costs. For each permutation in the queue, we try all its
        % neighbors and enqueue the ones that decrease the cost
        % w.r.t the current permutation (not the global cost). If
        % the current cost cannot be decreased we mark it as a
        % local minimum. We stop when all elements of the queue are
        % local minima.
        %
        % The queue is circular. Enqueuing is sequential, dequeuing
        % is random access.
	% BUG: (small) sometimes ncand grows by 1 for no reason at all!

        
        % Set up the queue

        LMAX = 4*n;  % queue capacity
	NSEARCHMAX = 1e4;

        picand = zeros( n, LMAX );   % queue: permutations
        costcand = zeros( 1, LMAX ); % queue: cost
        isfullcand = zeros( 1, LMAX ); % queue: is this location
                                       % occupied?
        islocalmin = zeros( 1, LMAX ); % queue: is local min
        
        ncand = 0;  % number of elements stored
        ienq = 0;   % pointer for enqueuing
        ideq = 0;   % pointer for dequeuing (not really needed)
        
        enqueue_cand( pib, costb );
        enqueue_cand( pic, costc );
        
        nsearched = 2;
        while 1
            icrt = find( isfullcand & ~islocalmin ); % TODO: more efficient
            if ( isempty( icrt ) | (nsearched > NSEARCHMAX) ) break; end;   %%  LOOP EXIT
            icrt = icrt( 1 ); 
            costold = costcand( icrt );
            piold = picand( :, icrt )';
            for ii = 1:n-1; % all neighbors
                nsearched = nsearched + 1;
                tau = [1:ii-1 ii+1 ii ii+2:n]; % transposition ii
                pi = piold( tau );
                cost = sum(sum(tril(R(pi,pi), -1)));
                if cost < costold
%		   disp( [num2str(nsearched) '  cost=' num2str( ...
%                          cost ) ' ncand=' num2str( ncand ) ' ' num2str( ...
%                          isfullcand( icrt ))]);
                    if isfullcand( icrt ) % remove cand, because
                                          % a better neighbor found
                        dequeue_cand( icrt ); 
                    end;
                    if ncand == LMAX
                        compress_cand;
                    end;
                    if ncand < LMAX
                        enqueue_cand( pi, cost );
                    else % expand queue or exit
		         if LMAX < 3000
			    picand = [ picand zeros( n, LMAX )]; 
			    costcand = [ costcand zeros( 1, LMAX ) ];
			    isfullcand = [ isfullcand zeros( 1, LMAX ) ]; 
			    islocalmin = [ islocalmin zeros( 1, LMAX ) ]; 
			    LMAX = 2*LMAX;
			    ienq = ncand;
			    LMAX
			 else
			    break;  % EXIT LOOP
			 end;
                    end;
                end;
            end;
            if isfullcand( icrt ) % no neigh was better,therefore...
                islocalmin( icrt ) = 1;
            end; % for
        end; % while
        
        ilocalmin = find( isfullcand & islocalmin );
        [ cost, imin ] = min( costcand( ilocalmin ));
        sigma = picand( :, imin )';
%nsearched
%imin
        
    else  % No local search
        if costb < costc
            sigma = pib;
        else
            sigma = pic;
        end;
        cost = min( costb, costc );
    end;
    
function enqueue_cand( pi, cost )
    
    global ncand ienq picand costcand isfullcand islocalmin LMAX
      % ASSUMES: queue is NOT FULL. Crashes if called with full queue.
    
    % search for a free place
 
    ienq = mod( ienq, LMAX )+1;
%    while (isfullcand( ienq ) & (ienq < LMAX))
%        ienq = ienq+1;
%    end;
    ii = find( isfullcand( ienq:end )==0 );
    if isempty( ii )
       ienq = min( find( isfullcand( 1:ienq-1 )==0 ));
    else
       ienq = min( ii );
    end; 
       
    picand( :, ienq ) = pi';
    costcand( ienq ) = cost;
    isfullcand( ienq ) = 1;
    islocalmin( ienq ) = 0;
    ncand = ncand+1;
    
function dequeue_cand( idq )
    
    global ncand ienq picand costcand isfullcand islocalmin LMAX
    
    isfullcand( idq ) = 0;
    ncand = ncand-1;

function compress_cand
    
    global ncand ienq picand costcand isfullcand islocalmin LMAX
    
    %keyboard
    il = find( islocalmin & isfullcand );	
    [cost, imin ] = min( costcand( il )) ;

    % deletes all local minima which are not optimal
    idel = setdiff( il, il( imin ) );
    islocalmin( idel ) = 0;
    isfullcand( idel ) = 0;
    ncand = ncand - length( idel );
    ienq = 0;
    
    