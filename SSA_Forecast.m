function [F_out] = ssa_forecast_analysis(Y, L, N, M, outfile)
% Performs SSA decomposition, exploratory analysis, and forecasting.
% Inputs:
%   Y       - Time series (column vector)
%   L       - Embedding dimension (window length)
%   N       - Number of components for forecasting
%   M       - Number of steps to forecast
%   outfile - Excel file (.xlsx) to write output

%% Checking input format
if isrow(Y)
    Y = Y'; 
end

T = length(Y);
K = T - L + 1;

%% 1. Trajectory Matrix Construction
X = zeros(L, K);
for i = 1:K
    X(:,i) = Y(i:i+L-1);
end

%% 2. Covariance Matrix (Toeplitz-like estimate)
covX = xcorr(Y, L-1, 'unbiased');
C = toeplitz(covX(L:end));  % L x L covariance matrix

%% 3. Eigen-Decomposition
[RHO, LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);
[LAMBDA, idx] = sort(LAMBDA, 'descend');
RHO = RHO(:, idx);  % sorting eigenvectors accordingly

% Cumulative Variance
total_lambda = sum(LAMBDA);
cum_var = LAMBDA / total_lambda;
cum_var_pct = cum_var * 100;

%% 4. Principal Components
PC = X.' * RHO;

%% 5. Reconstructed Components via anti-diagonal averaging
RC = zeros(T, L);
for m = 1:L
    outer = PC(:,m) * RHO(:,m)';
    outer = flipud(outer);
    for n = 1:T
        RC(n,m) = mean(diag(outer, -(K-1)+n));
    end
end

%% 6. Forecasting via Linear Recurrent Formula
% Linear recurrent coefficients
A = zeros(L-1, 1);
for i = 1:N
    A = A + RHO(L, i) * RHO(1:L-1, i);
end
v = norm(RHO(L, 1:N));
A = A / (1 - v^2);
A = flipud(A); % for correct indexing in recurrence

% Initial reconstructed signal (using first N components)
G = sum(RC(:,1:N), 2);

% Forecasting
F = zeros(T + M, 1);
F(1:T) = G;
for i = T+1:T+M
    F(i) = A' * F(i-L+1:i-1);
end
F_out = F; % original + forecast

%% 7. Exporting to Excel 
writematrix(C, outfile, 'Sheet', 'Covariance', 'Range', 'A1');
writematrix(LAMBDA, outfile, 'Sheet', 'Eigenvalues', 'Range', 'A1');
writematrix(cum_var_pct, outfile, 'Sheet', 'Cumulative Variance (%)', 'Range', 'A1');
writematrix(PC, outfile, 'Sheet', 'Principal Components', 'Range', 'A1');
writematrix(RC, outfile, 'Sheet', 'Reconstructed Components', 'Range', 'A1');
writematrix(F_out, outfile, 'Sheet', 'Forecast Output', 'Range', 'A1');

% Plot 1: SSA decomposition 
figure;
subplot(3,1,1); plot(Y); title('Original Temperature Time Series'); ylabel('Temp');
subplot(3,1,2); plot(LAMBDA, '-o'); title('Eigenvalues'); ylabel('\lambda');
subplot(3,1,3); plot(cumsum(cum_var_pct), '-s'); title('Cumulative Variance (%)');
xlabel('Component'); ylabel('%');

% Plot 2: Forecasted Values vs Original Values
figure;
plot(1:T, Y, 'b', 'LineWidth', 1.5); hold on;
plot(1:T, F_out(1:T), 'g--', 'LineWidth', 1.5); % reconstructed
plot(T+1:T+M, F_out(T+1:end), 'r', 'LineWidth', 2); % forecasted
legend('Original', 'Reconstructed', 'Forecast');
xlabel('Time'); ylabel('Temperature');
title('SSA-Temperature Forecast');
grid on;

end