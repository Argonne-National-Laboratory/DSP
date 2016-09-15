%% Figure 6
clc
clear

gra = csvread('output/sslp_15_45_15_3_dualvars.csv');
spx = csvread('output/sslp_15_45_15_0_dualvars.csv');
ipm = csvread('output/sslp_15_45_15_2_dualvars.csv');

h = figure;
% set(h, 'position', [0, 0, 500, 300]);
plot(1:size(gra,1), gra, 'Color', [0 0.4 0]);
hold on;
plot(1:size(spx,1), spx, 'r');
plot(1:size(ipm,1), ipm, 'b');
legend('DDSub', 'DDCP', 'DSP');
xlabel('Number of Iterations');
ylabel('||\lambda^{k+1} - \lambda^k||_2');
xlim([1 100]);
hold off;
