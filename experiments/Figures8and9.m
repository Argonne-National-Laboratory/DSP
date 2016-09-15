%% Figures 8 and 9
clc
clear
ncores=[5,10,25,50,100,250];
probs={...
    'sslp_10_50_50',...
    'sslp_10_50_100',...
    'sslp_10_50_500',...
    'sslp_10_50_1000',...
    'sslp_10_50_2000'
    };
nprobs = length(probs);
nscales=[4,5,6,6,6];
times=zeros(5,6,nprobs);

% Read results
strong_scales=zeros(6,nprobs);
for i = 1:nprobs
    for j = 1:nscales(i)
        fname = sprintf('output/%s_scale%d.csv',probs{i},ncores(j));
        if exist(fname, 'file')
            d = importdata(fname);
            times(:,j,i) = d(5:9);
        end
    end
    strong_scales(1:nscales(i),i) = times(1,1,i) ./ times(1,1:nscales(i),i);
end

% Figure 8
figure('position',[100 100 700 400]);
loglog(ncores, ncores ./ ncores(1), '-*k')
hold on
loglog(ncores, strong_scales, '-*')
hold off
xlim([3 500])
ylim([0.1,100])
set(gca,'fontsize',14)
set(gca,'ygrid','on')
set(gca,'xtick',ncores)
xlabel('Cores')
ylabel('Speedup')
strprobs = strrep(probs, '_', '\_');
legend({'Linear Speedup', strprobs{1:5}}, 'location', 'eastoutside')

% Figure 9
profile = times(:,:,5);
comm_time = profile(1,:) - sum(profile(2:5,:),1);
profile = [profile(1,:); comm_time; profile(2:5,:)]';

norm_profile = zeros(size(profile));
for i = 1:size(profile,1)
    norm_profile(i,:) = profile(i,:) / profile(i,1);
end

figure('position',[100 100 700 400]);
set(gca,'fontsize',14)
bar(norm_profile(:,2:6)*100,'stacked')
legend('Communication','Master Solution','Lower Bounding',...
    'Upper Bounding','Cut Generation','location','northeastoutside')
xlabel('Cores')
ylabel('Time Spent (%)')
set(gca,'XTickLabel',{'5','10','25','50','100','250'})