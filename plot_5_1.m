alpha = 10.0;

% Brownian dynamics
dt     = [0.001 0.002 0.003 0.004];
x      = dt / alpha;
epot   = [-0.810393 -0.806041 -0.802471 -0.801589];
err    = [ 0.0011554 0.00107722 0.00100416 0.00123882];

% Monte Carlo, plotted at dt/alpha = 0
x_mc    = 0;
epot_mc = -0.814095;
err_mc  = 0.00058108;

figure;
errorbar(x, epot, err, 'o-', 'LineWidth', 1.2, 'MarkerFaceColor', 'auto');
hold on;
errorbar(x_mc, epot_mc, err_mc, 's', 'LineWidth', 1.2, 'MarkerFaceColor', 'auto');
hold off;

xlim([-0.00003 0.00043]);
xlabel('\Delta t / \alpha');
ylabel('Potential energy per particle');
title('Brownian dynamics vs Monte Carlo (N = 64, \rho = 0.3, T = 1.0)');
legend('Brownian dynamics (\alpha = 10)', 'Monte Carlo', 'Location', 'southeast');
grid on;
