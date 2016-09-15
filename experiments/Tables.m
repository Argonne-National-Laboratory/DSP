%% Tables 2, 4, 5, and 6

clc
clear
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
    'sslp_10_50_50',...
    'sslp_10_50_100',...
    'sslp_10_50_500',...
    'sslp_10_50_1000',...
    'sslp_10_50_2000'...
    'sslp_15_45_5',...
    'sslp_15_45_10',...
    'sslp_15_45_15'
    };

fileID = fopen('Tables.tex', 'w');

fprintf(fileID, '\\documentclass[11pt, oneside]{article}\n');
fprintf(fileID, '\\usepackage[top=1in, bottom=1in, left=1in, right=1in]{geometry}\n');
fprintf(fileID, '\\usepackage{graphicx,amssymb}\n');
fprintf(fileID, '\\geometry{letterpaper}\n');
fprintf(fileID, '\\begin{document}\n');

% Table 2

fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Computational results for SIPLIB instances using different dual decomposition methods.}\n');
fprintf(fileID, '\\label{tab:siplib:recovery}\n');
fprintf(fileID, '\\centering\\scriptsize\n');
fprintf(fileID, '\\begin{tabular}{lrlrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '  Instance & Scen ($r$) & Method & Iter & UB  & LB & Gap (\\%%) & Time (s)    \\\\ \n');
fprintf(fileID, '  \\hline\n');
for p = probs
    fname0 = sprintf('output/%s_0.csv',p{1});
    fname1 = sprintf('output/%s_2.csv',p{1});
    fname3 = sprintf('output/%s_3.csv',p{1});
    i = 0;
    for f = {fname3, fname0, fname1}
        i = i + 1;
        if exist(f{1}, 'file')
            d = importdata(f{1});
            
            % Parse file name
            split_idx = find(p{1} == '_', 1, 'last');
            instance_name = p{1}(1:(split_idx-1));
            num_scenarios = p{1}((split_idx+1):length(p{1}));
            if i == 1
                fprintf(fileID, '  {\\tt %s} & %s', strrep(instance_name, '_', '\_'), num_scenarios);
            else
                fprintf(fileID, '            & ');
            end
            
            if i == 1
                fprintf(fileID, ' & DDSub');
            elseif i == 2
                fprintf(fileID, ' & DDCP');
            else
                fprintf(fileID, ' & DSP');
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
    fprintf(fileID, ' \\hline');
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');


scen = [4,8,16,32,64];

% Table 4
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Numerical results for stochastic unit commitment problem under DSP.}\n');
fprintf(fileID, '\\label{tab:suc:scale}\n');
fprintf(fileID, '\\centering\\scriptsize\n');
fprintf(fileID, '\\begin{tabular}{rrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '  Scen ($r$) & Iter & LB & UB & Gap (\\%%) & Time (s) \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scen
    f = sprintf('output/suc%d_21.csv',s);
    if exist(f,'file')
        d = importdata(f);
        iter = d(1) + 1;
        ub = d(2);
        lb = d(3);
        time = d(4);
        gap = (ub-lb)/ub*100;
    else
        fprintf('File %s does not exist.\n', f0);
    end
    fprintf(fileID, '  %d & %d & %.1f & %.1f & ', s, iter, lb, ub);
    if gap < 0.01
        fprintf(fileID, '$<$');
    end
    fprintf(fileID, '%.2f & ', gap);
    if time > 21600
        fprintf(fileID, '$>$');
    end
    fprintf(fileID, '%.0f \\\\\n', time);
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');


% Table 5
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Upper bounds and lower bounds of the stochastic unit commitment problem resulting from DSP with and without Procedure 1 at the first iteration.}\n');
fprintf(fileID, '\\label{tab:procedure1}\n');
fprintf(fileID, '\\centering\\scriptsize\n');
fprintf(fileID, '\\begin{tabular}{rrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '       & \\multicolumn{2}{r}{DSP-P1} & \\multicolumn{2}{r}{DSP+P1} & \\\\\n');
fprintf(fileID, '  \\cline{2-5}\n');
fprintf(fileID, '  Scen ($r$) & LB & UB & LB & UB & LB Improved (\\%%) \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scen
    f0 = sprintf('output/suc%d_2-1_oneiter.csv',s);
    f1 = sprintf('output/suc%d_21_oneiter.csv',s);
    if exist(f0,'file')
        d = importdata(f0);
        ub0 = d(2);
        lb0 = d(3);
    else
        fprintf('File %s does not exist.\n', f0);
    end
    if exist(f1,'file')
        d = importdata(f1);
        ub1 = d(2);
        lb1 = d(3);
    else
        fprintf('File %s does not exist.\n', f0);
    end
    improve = (lb1-lb0)/lb0*100;
    fprintf(fileID, '  %d & %.1f', s, lb0);
    if ub0 > 1e+10
        fprintf(fileID, ' & $\\infty$');
    else
        fprintf(fileID, ' & %.1f', ub0);
    end
    fprintf(fileID, ' & %.1f', lb1);
    if ub1 > 1e+10
        fprintf(fileID, ' & $\\infty$');
    else
        fprintf(fileID, ' & %.1f', ub1);
    end
    fprintf(fileID, ' & %.2f \\\\\n', improve);
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');


% Table 6
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Computational results from DDCP (Algorithm 2) and DSP (Algorithm 3) with and without Procedure 1.}\n');
fprintf(fileID, '\\label{tab:procedure1}\n');
fprintf(fileID, '\\centering\\scriptsize\n');
fprintf(fileID, '\\begin{tabular}{rrrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '  Scen ($r$) & Iter & Method & LB & UB & Gap (\\%%) & Time (s) \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scen
    f0 = sprintf('output/suc%d_0-1.csv',s);
    f1 = sprintf('output/suc%d_01.csv',s);
    f2 = sprintf('output/suc%d_2-1.csv',s);
    f3 = sprintf('output/suc%d_21.csv',s);
    if exist(f0,'file')
        d = importdata(f0);
        iter = d(1) + 1;
        ub = d(2);
        lb = d(3);
        time = d(4);
        gap = (ub-lb)/ub*100;
    else
        fprintf('File %s does not exist.\n', f0);
    end
    fprintf(fileID, '  %d & %d & DDCP-P1 & %.1f & ', s, iter, lb);
    if ub > 1e+10
        fprintf(fileID, '$\\infty$ & ');
    else
        fprintf(fileID, '%.1f & ', ub);
    end
    if gap > 1.0
        fprintf(fileID, '$\\infty$ & ');
    else
        if gap < 0.01
            fprintf(fileID, '$<$');
        end
        fprintf(fileID, '%.2f & ', gap);
    end
    if time > 21600
        fprintf(fileID, '$>21600$ \\\\\n');
    else
        fprintf(fileID, '%.0f \\\\\n', time);
    end
    
    if exist(f1,'file')
        d = importdata(f1);
        iter = d(1) + 1;
        ub = d(2);
        lb = d(3);
        time = d(4);
        gap = (ub-lb)/ub*100;
    else
        fprintf('File %s does not exist.\n', f1);
    end
    fprintf(fileID, '   & %d & DDCP+P1 & %.1f & ', iter, lb);
    if ub > 1e+10
        fprintf(fileID, '$\\infty$ & ');
    else
        fprintf(fileID, '%.1f & ', ub);
    end
    if gap > 1.0
        fprintf(fileID, '$\\infty$ & ');
    else
        if gap < 0.01
            fprintf(fileID, '$<$');
        end
        fprintf(fileID, '%.2f & ', gap);
    end
    if time > 21600
        fprintf(fileID, '$>21600$ \\\\\n');
    else
        fprintf(fileID, '%.0f \\\\\n', time);
    end
    
    if exist(f2,'file')
        d = importdata(f2);
        iter = d(1) + 1;
        ub = d(2);
        lb = d(3);
        time = d(4);
        gap = (ub-lb)/ub*100;
    else
        fprintf('File %s does not exist.\n', f2);
    end
    fprintf(fileID, '   & %d & DSP-P1 & %.1f & ', iter, lb);
    if ub > 1e+10
        fprintf(fileID, '$\\infty$ & ');
    else
        fprintf(fileID, '%.1f & ', ub);
    end
    if gap > 1.0
        fprintf(fileID, '$\\infty$ & ');
    else
        if gap < 0.01
            fprintf(fileID, '$<$');
        end
        fprintf(fileID, '%.2f & ', gap);
    end
    if time > 21600
        fprintf(fileID, '$>21600$ \\\\\n');
    else
        fprintf(fileID, '%.0f \\\\\n', time);
    end
    
    if exist(f3,'file')
        d = importdata(f3);
        iter = d(1) + 1;
        ub = d(2);
        lb = d(3);
        time = d(4);
        gap = (ub-lb)/ub*100;
    else
        fprintf('File %s does not exist.\n', f3);
    end
    fprintf(fileID, '   & %d & DSP+P1 & %.1f & ', iter, lb);
    if ub > 1e+10
        fprintf(fileID, '$\\infty$ & ');
    else
        fprintf(fileID, '%.1f & ', ub);
    end
    if gap > 1.0
        fprintf(fileID, '$\\infty$ & ');
    else
        if gap < 0.01
            fprintf(fileID, '$<$');
        end
        fprintf(fileID, '%.2f & ', gap);
    end
    if time > 21600
        fprintf(fileID, '$>21600$ \\\\\n');
    else
        fprintf(fileID, '%.0f \\\\\n', time);
    end
    fprintf(fileID, '  \\hline\n');
end
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');


% Table 7
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Numerical results for the extensive form of the stochastic unit commitment problems}\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\begin{tabular}{rrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '            & Number of   &             &             &           & Wall time \\\\\n');
fprintf(fileID, '  Scenarios & B\\&C nodes & Upper Bound & Lower Bound & Gap (\\%%) & (sec.)    \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scen
    fname = sprintf('output/suc%d-de.csv',s);
    if exist(fname, 'file')
        d = importdata(fname);
        fprintf(fileID, '  %d & %d', s, d(5));
        if d(1) > 1.0e+20
            fprintf(fileID, ' & $\\infty$');
        else
            fprintf(fileID, ' & %.1f', d(1));
        end
        fprintf(fileID, ' & %.1f', d(2));
        if d(1) > 1.0e+20
            fprintf(fileID, ' & $\\infty$');
        else
            gap = abs(d(1) - d(2)) / abs(d(1)) * 100;
            if gap < 0.01
                fprintf(fileID, ' & $<0.01$');
            else
                fprintf(fileID, ' & %.2f', gap);
            end
        end
        if d(3) > 21600
            fprintf(fileID, ' & $>21600$');
        else
            fprintf(fileID, ' & %.0f', d(3));
        end
        fprintf(fileID, ' \\\\\n');
    end
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');


fprintf(fileID, '\\end{document}\n');
fclose(fileID);

%%
clc
clear
