function ret = rand()
persistent a;
if isempty(a)
    a = 0;
end
a  = a + 0.12;
if a > 1
    a = 0;
end
ret = a;
