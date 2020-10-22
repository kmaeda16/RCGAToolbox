function c = RCGAgetNewChild(p1, p2, p3)
% RCGAgetNewChild generates a child with UNDX and checks wheather it is
% within the search space.
% 
% [SYNTAX]
% c = RCGAgetNewChild(p1, p2, p3)
% 
% [INPUT]
% p1 :  First main parent (chrom structure).
% p2 :  Second main parent (chrom structure).
% p3 :  Sub-parent (chrom structure).
% 
% [OUTPUT]
% c  :  Generated child (chrom structure).


%% Default number of trials of generating children
maxitr = 100; % You can change this line


%% Checking whether c is valid
for i = 1 : maxitr
    c = RCGA_UNDX(p1,p2,p3);
    c.f = Inf;
    c.g = Inf(size(p1.g));
    c.phi = Inf;
    if max(isnan(c.gene))
        warning('Generated child has nan in its genes!');
    end
    if min( 0 <= c.gene & c.gene <= 1 )
        return;
    end
end


%% If genes out of bounds, put them bounds
c.gene(1<c.gene) = 1;
c.gene(c.gene<0) = 0;
