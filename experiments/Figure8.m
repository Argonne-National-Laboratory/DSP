%% Figure 5

d1 = csvread('output/suc8-3-1_duals.csv');
p1 = csvread('output/suc8-3-1_primals.csv');
d2 = csvread('output/suc8-21_duals.csv');
p2 = csvread('output/suc8-21_primals.csv');
d3 = csvread('output/suc8-2-1_duals.csv');
p3 = csvread('output/suc8-2-1_primals.csv');

d2 = d2(1);
p2 = p2(1);

ninf1 = sum(p1 > 1.0e+20);
ninf2 = sum(p3 > 1.0e+20);

LW = 1;
plot((ninf1+1):size(p1,1), p1((ninf1+1):size(p1,1),:), 'o-', 'Color', [0 0.4 0], 'linewidth', LW);
hold on;
plot(1:size(d1,1), d1, '*--', 'Color', [0 0.4 0], 'linewidth', LW);
plot((ninf2+1):size(p3,1), p3((ninf2+1):size(p3,1),:), 'ro-', 'linewidth', LW);
plot(1:size(d3,1), d3, 'r*--', 'linewidth', LW);
plot(1:size(p2,1), p2, 'bo-', 'linewidth', LW);
plot(1:size(d2,1), d2, 'b*--', 'linewidth', LW);
ylabel('Objective Value');
xlabel('Number of Iterations');
legend('DDSub (Upper Bounds)', ...
    'DDSub (Lower Bounds)', ...
    'DSP-P1 (Upper Bounds)', ...
    'DSP-P1 (Lower Bounds)', ...
    'DSP+P1 (Upper Bounds)', ...
    'DSP+P1 (Lower Bounds)', ...
    'Location', 'NorthWest');
hold off;

% save
hgexport(gcf, 'Figure8.eps');


