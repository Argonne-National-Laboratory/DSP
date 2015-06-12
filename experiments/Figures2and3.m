%% Figures 2 and 3
clc
probs={...
    'dcap233_200',...
    'dcap233_300',...
    'dcap233_500',...
    'dcap243_200',...
    'dcap243_300',...
    'dcap243_500',...
    'dcap332_200',...
    'dcap332_300',...
    'dcap332_500',...
    'dcap342_200',...
    'dcap342_300',...
    'dcap342_500',...
    'sslp_5_25_50',...
    'sslp_5_25_100',...
    'sslp_15_45_5',...
    'sslp_15_45_10',...
    'sslp_15_45_15',...
    'sslp_10_50_50',...
    'sslp_10_50_100'
    };

nprobs = length(probs);

iter = zeros(nprobs, 3);
gap  = zeros(nprobs, 3);
time = zeros(nprobs, 3);

i = 1;
for p = probs
    fname0 = sprintf('output/%s_s0.csv',p{1});
    fname2 = sprintf('output/%s_s2.csv',p{1});
    fname3 = sprintf('output/%s_s3.csv',p{1});
    j = 1;
    for f = {fname3, fname0, fname2}
        if exist(f{1}, 'file')
            d = importdata(f{1});
            iter(i,j) = d(2);
            ub   = d(3);
            lb   = d(4);
            time(i,j) = d(5);
            gap(i,j)  = abs(ub - lb) / abs(ub) * 100;
        end
        j = j + 1;
    end
    i = i + 1;
end

figure();
axes('position', [0.08 0.05 0.92 0.9]);
bar(gap);
ylim([0 10]);
ylabel('Optimality Gap (%)');
Xlabs = strrep(probs, '_', '\_');
xticklabel_rotate(1:size(gap,1), 90, Xlabs);
legend('DDSub', 'DDCP', 'DSP');

% save
hgexport(gcf, 'Figure2.eps');

rel_time = time ./ repmat(time(:,3), 1, 3);
rel_iter = iter ./ repmat(iter(:,3), 1, 3);

figure();
axes('position', [0.06 0.05 0.92 0.9]);
rel = [rel_time(:,1:2) rel_iter(:,1:2)];
bar(rel(:,[1,3,2,4]));
set(gca,'YScale','log')
Xlabs = strrep(probs, '_', '\_');
xticklabel_rotate(1:size(gap,1), 90, Xlabs);
legend('DDSub (Solution Time)', 'DDSub (# of Iterations)', ...
    'DDCP (Solution Time)', 'DDCP (# of Iterations)', ...
    'location', 'southeast');

% save
hgexport(gcf, 'Figure3.eps');