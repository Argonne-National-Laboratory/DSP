%% Table 5

clc
scenarios=[4,8,16,32,64];

fileID = fopen('Table5.tex', 'w');

fprintf(fileID, '\\documentclass[11pt, oneside]{article}\n');
fprintf(fileID, '\\usepackage{geometry,graphicx,amssymb}\n');
fprintf(fileID, '\\geometry{letterpaper}\n');
fprintf(fileID, '\\begin{document}\n');
fprintf(fileID, '\\begin{table}[htpb]\n');
fprintf(fileID, '\\caption{Numerical results for the extensive form of the stochastic unit commitment problems}\n');
fprintf(fileID, '\\centering\n');
fprintf(fileID, '\\begin{tabular}{rrrrrr}\n');
fprintf(fileID, '  \\hline\n');
fprintf(fileID, '            & Number of   &             &             &           & Wall time \\\\\n');
fprintf(fileID, '  Scenarios & B\\&C nodes & Upper Bound & Lower Bound & Gap (\\%%) & (sec.)    \\\\\n');
fprintf(fileID, '  \\hline\n');
for s = scenarios
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
