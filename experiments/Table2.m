%% Table 2
% This MATLAB script will generate a latex file with Table 2.

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

fileID = fopen('Table2.tex', 'w');

fprintf(fileID, '\\documentclass[11pt, oneside]{article}\n');
fprintf(fileID, '\\usepackage{geometry,graphicx,amssymb}\n');
fprintf(fileID, '\\geometry{letterpaper}\n');
fprintf(fileID, '\\begin{document}\n');
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Computational results for SIPLIB instances using different dual decomposition methods.}\n');
fprintf(fileID, '\\label{tab:siplib:recovery}\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\begin{tabular}{lrlrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '           &           &        &            & Upper & Lower  & Gap  & Wall time \\\\\n');
fprintf(fileID, '  Instance & Scenarios & Method & Iter & Bound  & Bound & (\\%%) & (sec.)    \\\\ \n');
fprintf(fileID, '  \\hline\n');
for p = probs
    fname0 = sprintf('output/%s_s0.csv',p{1});
    fname1 = sprintf('output/%s_s2.csv',p{1});
    fname3 = sprintf('output/%s_s3.csv',p{1});
    i = 0;
    for f = {fname3, fname0, fname1}
        i = i + 1;
        if exist(f{1}, 'file')
            d = importdata(f{1});
            
            % Parse file name
            split_idx = find(p{1} == '_', 1, 'last');
            instance_name = p{1}(1:(split_idx-1));
            num_scenarios = p{1}((split_idx+1):length(p{1}));
            fprintf(fileID, '  {\\tt %s} & %s', strrep(instance_name, '_', '\_'), num_scenarios);
            
            if i == 1
                fprintf(fileID, ' & Subgradient');
            elseif i == 2
                fprintf(fileID, ' & CPM');
            else
                fprintf(fileID, ' & IPM');
            end
            
            gap = abs(d(3) - d(4)) / abs(d(3)) * 100;
            fprintf(fileID, ' & %d & %.2f & %.2f & %.2f & ', ...
                d(2), d(3), d(4), gap);
            if d(5) < 21600
                fprintf(fileID, '%.0f', d(5));
            else
                fprintf(fileID, '$>$ 21600');
            end
            fprintf(fileID, ' \\\\\n');
        end
    end
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');
fprintf(fileID, '\\end{document}\n');
fclose(fileID);