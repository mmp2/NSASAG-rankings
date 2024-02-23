global Q0

Q0 = triu( ones( 100 ), 1 );

isCritch = 0;

pi1 = [ 1 3 4 6 2 5]
pi2 = [ 2 7 1 3 8 9 ]

sig=randperm(30);
pi1 = sig( pi1 )
pi2 = sig( pi2 )


dKdiff( pi1, pi2, 0 )

