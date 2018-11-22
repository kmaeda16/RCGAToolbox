function f = k_tablet(x)

n = length(x);
k = 0.25 * n;

f = sum( x(1:k) .^ 2 ) + sum( ( 100 * x(k+1:n) ) .^ 2 );
