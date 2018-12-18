function ret = randn()
persistent a;
if isempty(a)
    a = -1;
end
a  = a + 0.1;
if a > 1
    a = -1;
end
ret = a;
