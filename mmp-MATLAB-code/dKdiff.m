function [ d ] = dKdiff( pi1, pi2, isCritch )
 
% function [ d ] = dKdiff( pi1, pi2, isCritch )
% 
% Extension of Kendall distance to partial rankings with different
% sets of elements
%
% pi1, pi2 = row vectors of partial rankings with any integer entries,
%	     any lengths
% isCritch = 1  use Critchlow's formula
%	   = 0  computes min distance between sets

mm = max( max(pi1), max(pi2));

a = intersect( pi1, pi2 );
inda1 = find( ismember( pi1, a ));
inda2 = find( ismember( pi2, a ));
indb1 = find( ~ismember( pi1, a ));
indb2 = find( ~ismember( pi2, a ));

na = length( a );
nb1 = length( indb1 );
nb2 = length( indb2 );

da = dK( pi1( inda1 ), pi2( inda2 ), mm, na );
db = nb1*nb2 ;
d1 = nb1*length( pi1 )-sum( indb1 )- (1-isCritch)*nb1*(nb1-1)/2;
d2 = nb2*length( pi2 )-sum( indb2 )- (1-isCritch)*nb2*(nb2-1)/2;

d = da + db + d1 + d2;

