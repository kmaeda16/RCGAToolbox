figure('Position',[100 400 600 250]);

X = readtable('Summary.xlsx');

Y = reshape(X{:,'Mean'},[3,6]);
Y = Y';
% Y = flip(Y);
b = bar(Y*1e+3);
% b = barh(Y*1e+3);
b(1,1).FaceColor = [0 0.4470 0.7410];
b(1,2).FaceColor = [0.8500 0.3250 0.0980];
b(1,3).FaceColor = [0.9290 0.6940 0.1250];

hold on;


Z = reshape(X{:,'SD'},[3,6]);
Z  = Z';

dev = 0.225;
x = [ 1 - dev, 1, 1 + dev;
    2 - dev, 2, 2 + dev;
    3 - dev, 3, 3 + dev;
    4 - dev, 4, 4 + dev;
    5 - dev, 5, 5 + dev;
    6 - dev, 6, 6 + dev;
    ];
er = errorbar(x,Y*1e+3,Z*1e+3);
er(1,1).Color = [0 0 0];
er(1,2).Color = [0 0 0];
er(1,3).Color = [0 0 0];
er(1,1).LineStyle = 'none';
er(1,2).LineStyle = 'none';
er(1,3).LineStyle = 'none';

xlim([0.5 6.5]);
ylim([0 70]);

set(gca,'FontSize',12,'FontName','Arial');
ylabel('Time (ms)','FontSize',12,'FontName','Arial');
set(gca,'TickLength',[0.015 0.015]);
% set(gca,'XTick',10.^(-10:1:10));
% set(gca,'YTick',10.^(-12:2:10));
