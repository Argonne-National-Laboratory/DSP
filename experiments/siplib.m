%% Basic results
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

rel_time = time ./ repmat(time(:,1), 1, 3);
figure();
axes('position', [0.08 0.05 0.92 0.9]);
bar(rel_time(:,2:3));
%ylim([0 1]);
%ylabel('Solution Time Relative to Subgradient Method');
Xlabs = strrep(probs, '_', '\_');
xticklabel_rotate(1:size(time,1), 90, Xlabs);
legend('DDCP', 'DSP');

rel_iter = iter ./ repmat(iter(:,1), 1, 3);
figure();
axes('position', [0.08 0.05 0.92 0.9]);
bar(rel_iter(:,2:3));
ylim([0 1]);
%ylabel('Number of Iterations Relative to Subgradient Method');
Xlabs = strrep(probs, '_', '\_');
xticklabel_rotate(1:size(gap,1), 90, Xlabs);
legend('DDCP', 'DSP');
