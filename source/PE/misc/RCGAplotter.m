function RCGAplotter(T_sim,Y_sim,T_exp,Y_exp,t_label,y_label,statename)
% RCGAplotter plots simulation and experimental data
% 
% [SYNTAX]
% RCGAplotter(generation,t_sim,x_sim,t_exp,x_exp,t_label,x_label,statename)
% 
% [INPUT]
% T_sim     :  Simulated time vector (column)
% Y_sim     :  Simulated variable vector (column)
% T_exp     :  Experimental time vector (column)
% Y_exp     :  Experimental variable vector (column)
% t_label   :  Label for time (x-axis in the plot)
% y_label   :  Label for variables (y-axis in the plot)
% statename :  Cell array for state variable names


persistent Y_prv;
if isequal(Y_prv,Y_sim)
    return;
else
    Y_prv = Y_sim;
end

plot(T_exp,Y_exp,'o','LineWidth',2);
set(gca,'FontSize',10,'FontName','Arial');
legend(statename);
hold on;
ax = gca;
ax.ColorOrderIndex = 1;
plot(T_sim,Y_sim,'-','LineWidth',2);
hold off;

xlabel(t_label);
ylabel(y_label);

T = [ T_sim; T_exp ];
xlim( [ min(T) max(T) ] );
ylim( [ min(min(Y_exp)) max(max(Y_exp)) ] );


