function RCGAplotter(T_sim, Y_sim, T_exp, Y_exp, t_label, y_label, statename)
% RCGAplotter plots simulation and experimental data.
% 
% [SYNTAX]
% RCGAplotter(T_sim, Y_sim, T_exp, Y_exp, t_label, y_label, statename)
% 
% [INPUT]
% T_sim     :  Simulated time vector (row: time).
% Y_sim     :  Simulated state matrix (row: time; column: state).
% T_exp     :  Experimental time vector (row: time).
% Y_exp     :  Experimental state matrix (row: time; column: state).
% t_label   :  Label for time (x-axis in the plot).
% y_label   :  Label for states (y-axis in the plot).
% statename :  Cell array of state names.


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
