# Temperature Time Series Forecasting using Singular Spectrum Analysis (SSA)

This study uses **Singular Spectrum Analysis (SSA)** for exploratory analysis and forecasting of time series univariate atmospheric temperature data. 



## Key Objectives:

- Analyze complex temperature time series data.
- Identify components like trend and oscillations.
- Forecast future temperature values using SSA-based methods.



## Methodology:

This study applies **Singular Spectrum Analysis (SSA)**, a non-parametric spectral decomposition method, to the temperature time series. The SSA method helps identify:

- Trend
- Seasonal/cyclic components
- Noise

The signal is reconstructed using selected components and extended through **linear recurrent relations** to generate forecasts.



## Model Pipeline:

1. **Input**: Read the temperature data from the input file.
2. **SSA Decomposition**:
   - Construct the trajectory matrix.
   - Calculate the covariance matrix.
   - Perform eigen-decomposition.
3. **Component Analysis**:
   - Principal component extraction.
   - Reconstruction through anti-diagonal averaging.
4. **Forecasting**:
   - Estimate recurrent coefficients from leading components.
   - Forecast future values using these coefficients.
5. **Output**:
   - Save covariance matrix, eigenvalues, reconstructed signal, and forecasted values to the output file.
   - Plot the original series, eigenvalues, cumulative variance, and the forecasted values.



## Challenges Addressed:

- **Non-stationarity**: SSA naturally handles trends and seasonal shifts.
- **Noise filtering**: Reconstructed components help in removing noise-dominated eigenvectors.
- **Short-term prediction**: Linear recurrent coefficients help in short-term forecasting.



## Results:

The SSA model successfully:
- Extracted interpretable components from the raw temperature series.
- Reconstructed the denoised signal.
- Produced accurate short-term forecasts.
- Provided eigenvalue plots and variance explanation for component selection.



## Impact:

The SSA's is non-parametric in nature. So, it's great for situations where we don't know much about the system or when the data changes quickly and often.



## Technology and Tools:

- **Language**: MATLAB
- **Functions Used**: `eig`, `toeplitz`, `xcorr`, `xlsread`, `writematrix`, matrix algebra
- **Visualization**: MATLAB's built-in `plot`, `subplot`
- **File I/O**: Reads/writes Excel `.xlsx` files



**Note:** The dataset used in this study is confidential. Hence, the dataset, output and the plots are not included.

