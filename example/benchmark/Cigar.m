function f = Cigar(x)

f = x(1) ^ 2 + sum( ( 1000 * x(2:end) ) .^ 2 );
