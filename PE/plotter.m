function plotter(generation,t_sim,x_sim,t_exp,x_exp,t_label,x_label,statename)

plot(t_exp,x_exp,'o','LineWidth',2);
set(gca,'FontSize',10,'FontName','Arial');
legend(statename);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(t_sim,x_sim,'-','LineWidth',2);
hold off;

xlabel(t_label);
ylabel(x_label);

T = [ t_sim; t_exp ];
xlim( [ min(T) max(T) ] );
ylim( [ min(min(x_exp)) max(max(x_exp)) ] );
title(sprintf('Generation = %d',generation));
drawnow;
