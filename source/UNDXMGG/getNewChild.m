function c = getNewChild(p1, p2, p3)
% Handler of UNDX function
% Checking wheather children are in search region 

maxitr = 100;

for i = 1 : maxitr
    c = UNDX(p1,p2,p3);
    c.f = Inf;
    c.g = Inf(size(p1.g));
    c.phi = Inf;
    if min( 0 <= c.gene & c.gene <= 1 )
        return;
    end
end

c.gene(1<c.gene) = 1;
c.gene(c.gene<0) = 0;
