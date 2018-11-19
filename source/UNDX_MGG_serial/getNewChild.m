function c = getNewChild(p1, p2, p3, SearchRegion)
% Handler of UNDX function
% Checking wheather children are in search region 

maxitr = 100;
flg = 1;

for i = 1:maxitr
    c = UNDX(p1,p2,p3);
    if 0 <= min(c.gene) && max(c.gene) <= 1
        flg = 0;
        break;
    end
end

if flg == 1
    for i = 1:length(c.gene)
        if 1 < c.gene(i)
            c.gene(i) = 1;
        elseif c.gene(i) < 0
            c.gene(i) = 0;
        end
    end
end

x = decodeGene2Variable(c,SearchRegion);
c.fitness = getFitness(x);
