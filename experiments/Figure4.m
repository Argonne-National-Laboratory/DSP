%% Figure 4

subobj1 = csvread('output/sslp_15_45_15_s0_sub_objvals.csv');
subobj2 = csvread('output/sslp_15_45_15_s2_sub_objvals.csv');
subobj3 = csvread('output/sslp_15_45_15_s3_sub_objvals.csv');
dual1 = csvread('output/sslp_15_45_15_s0_duals.csv');
dual2 = csvread('output/sslp_15_45_15_s2_duals.csv');
dual3 = csvread('output/sslp_15_45_15_s3_duals.csv');

LW = 1.0;

h = figure();
clf();
set(h, 'position', [0, 0, 900, 450]);

h1 = subplot(121);
set(h1, 'position', [0.08, 0.12, 0.45, 0.85]);
plot(1:size(subobj3,1), subobj3, '--', 'Color', [0 0.4 0], 'linewidth', LW);
hold on;
plot(1:size(dual3,1), dual3, '-', 'Color', [0 0.4 0], 'linewidth', LW);
plot(1:size(subobj1,1), subobj1, 'r--', 'linewidth', LW);
plot(1:size(dual1,1), dual1, 'r-', 'linewidth', LW);
plot(1:size(subobj2,1), subobj2, 'b--', 'linewidth', LW);
plot(1:size(dual2,1), dual2, 'b-', 'linewidth', LW);
set(gca, 'FontSize', 15);
xlabel('Number of Iterations');
ylabel('Objective Value');
xlim([1 99]);
ylim([-300 -255]);
hold off;

h2 = subplot(122);
set(h2, 'position', [0.54, 0.12, 0.45, 0.85]);
plot(1:size(subobj3,1), subobj3, '--', 'Color', [0 0.4 0], 'linewidth', LW);
hold on;
plot(1:size(dual3,1), dual3, '-', 'Color', [0 0.4 0], 'linewidth', LW);
plot(1:size(subobj1,1), subobj1, 'r--', 'linewidth', LW);
plot(1:size(dual1,1), dual1, 'r-', 'linewidth', LW);
plot(1:size(subobj2,1), subobj2, 'b--', 'linewidth', LW);
plot(1:size(dual2,1), dual2, 'b-', 'linewidth', LW);
set(gca, 'FontSize', 15);
legend('DDSub (Lower Bound)', ...
    'DDSub (Best Lower Bound)', ...
    'DDCP (Lower Bound)', ...
    'DDCP (Best Lower Bound)', ...
    'DSP (Lower Bound)', ...
    'DSP (Best Lower Bound)', ...
    'Location', 'SouthEast')
xlabel('Number of Iterations');
% ylabel('Objective Value');
% set(h2,'ytick',[]);
xlim([100 550]);
ylim([-300 -255]);
%set(h2,'xticklabel',[200 300 400 500]);
set(h2,'yticklabel',[]);
hold off;

% save
hgexport(gcf, 'Figure4.eps');
