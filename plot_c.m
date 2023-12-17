opts = delimitedTextImportOptions("NumVariables", 3);
opts.DataLines = [1, Inf];
opts.Delimiter = "\t";
opts.VariableNames = ["VarName1", "VarName2", "VarName3"];
opts.VariableTypes = ["double", "double", "double"];
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

tbl = readtable("./test_evap.txt", opts);
t = tbl.VarName1;
c_average = tbl.VarName2;
c_total = tbl.VarName3;

clear opts tbl

%% 
[fitresult, gof] = createFit(t, c_total);

%% 
function [fitresult, gof] = createFit(t, c_open)

    [xData, yData] = prepareCurveData( t, c_open );

    ft = fittype( 'exp1' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Lower = [-Inf -Inf];
    opts.StartPoint = [0 -0.2565304424498];
    opts.Upper = [Inf Inf];
    
    [fitresult, gof] = fit( xData, yData, ft, opts );
    
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'c_open vs. t', 'exp decay fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
    xlabel( 't', 'Interpreter', 'none' );
    ylabel( 'c_open', 'Interpreter', 'none' );
    grid on
end
