%% Table 5

clc
scenarios=[4,8,16,32,64];

fileID = fopen('Table5.tex', 'w');

fprintf(fileID, '\\documentclass[11pt, oneside]{article}\n');
fprintf(fileID, '\\usepackage{geometry,graphicx,amssymb}\n');
fprintf(fileID, '\\geometry{letterpaper}\n');
fprintf(fileID, '\\begin{document}\n');
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Numerical results for stochastic unit commitment problem under dual decomposition method.}\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\begin{tabular}{rrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '            & Number of   &             &             &           & Wall time \\\\\n');
fprintf(fileID, '  Scenarios & B\\&C nodes & Upper Bound & Lower Bound & Gap (\\%%) & (sec.)    \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scenarios
    fname = sprintf('output/suc%d-21.csv',s);
    if exist(fname, 'file')
        d = importdata(fname);
        fprintf(fileID, '  %d & %d & %.1f & %.1f', s, d(1)+1, d(2), d(3));
        gap = abs(d(2) - d(3)) / abs(d(2)) * 100;
        if gap < 0.01
            fprintf(fileID, ' & $<0.01$');
        else
            fprintf(fileID, ' & %.2f', gap);
        end
        if d(4) > 21600
            fprintf(fileID, ' & $>21600$');
        else
            fprintf(fileID, ' & %.0f', d(4));
        end
        fprintf(fileID, ' \\\\\n');
    end
end
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '\\end{tabular}\n');
fprintf(fileID, '\\end{table}\n');
fprintf(fileID, '\\end{document}\n');
fclose(fileID);
