Y = xlsread('inputdata.xlsx','range of the data');  
L = % window length;   
N = % number of SSA components;   
M = % Number of days forecasted;     
outfile = 'output_file.xlsx';

F = ssa_forecast_analysis(Y, L, N, M, outfile);
