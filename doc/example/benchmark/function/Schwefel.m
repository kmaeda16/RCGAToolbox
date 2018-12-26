function f = Schwefel(x)

n = length(x);
f = 418.9828873 * n;
for i = 1 : n
    f = f + x(i) * sin( sqrt( abs( x(i) ) ) );
end
