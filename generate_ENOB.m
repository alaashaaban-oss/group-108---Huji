%% ENOB extraction for 4-bit ADC
% Assumptions:
% - exported clock is the effective clock for the final outputs
% - sample outputs 5 ps after the rising edge
% - out<3> = MSB, out<0> = LSB

clear; clc; close all;

%% User settings
fs = 8e9;                  % sampling frequency [Hz]
logic_threshold = 0.5;     % threshold for deciding 0/1 [V]
sample_delay = 5e-12;      % sample 5 ps after the rising clock edge
M = 1024;                  % number of ADC samples used for ENOB
skip_edges = 10;           % ignore startup cycles

%% Read CSV files
clkData = readmatrix('clk.csv','NumHeaderLines',1);
d0Data  = readmatrix('out0.csv','NumHeaderLines',1);
d1Data  = readmatrix('out1.csv','NumHeaderLines',1);
d2Data  = readmatrix('out2.csv','NumHeaderLines',1);
d3Data  = readmatrix('out3.csv','NumHeaderLines',1);

t_clk = clkData(:,1); clk = clkData(:,2);
t_d0  = d0Data(:,1);  d0  = d0Data(:,2);
t_d1  = d1Data(:,1);  d1  = d1Data(:,2);
t_d2  = d2Data(:,1);  d2  = d2Data(:,2);
t_d3  = d3Data(:,1);  d3  = d3Data(:,2);

%% Find rising edges with threshold-crossing interpolation
clk_bin = clk > logic_threshold;
edge_idx = find(clk_bin(1:end-1)==0 & clk_bin(2:end)==1);

if isempty(edge_idx)
    error('No rising edges found in clk waveform.');
end

t_edge = zeros(size(edge_idx));

for k = 1:length(edge_idx)
    i = edge_idx(k);

    t1 = t_clk(i);
    t2 = t_clk(i+1);
    v1 = clk(i);
    v2 = clk(i+1);

    if v2 == v1
        t_edge(k) = t2;
    else
        t_edge(k) = t1 + (logic_threshold - v1) * (t2 - t1) / (v2 - v1);
    end
end

%% Build sampling times
t_sample_all = t_edge + sample_delay;

% Keep only sample times inside all waveform ranges
tmin = max([t_clk(1), t_d0(1), t_d1(1), t_d2(1), t_d3(1)]);
tmax = min([t_clk(end), t_d0(end), t_d1(end), t_d2(end), t_d3(end)]);
t_sample_all = t_sample_all(t_sample_all >= tmin & t_sample_all <= tmax);

if numel(t_sample_all) < skip_edges + M
    error('Not enough valid sample points. Run longer simulation or reduce M / skip_edges.');
end

t_sample = t_sample_all(skip_edges + (1:M));

%% Interpolate each output bit at the sampling times
v0 = interp1(t_d0, d0, t_sample, 'linear');
v1 = interp1(t_d1, d1, t_sample, 'linear');
v2 = interp1(t_d2, d2, t_sample, 'linear');
v3 = interp1(t_d3, d3, t_sample, 'linear');

%% Convert analog outputs to logic bits
b0 = v0 > logic_threshold;
b1 = v1 > logic_threshold;
b2 = v2 > logic_threshold;
b3 = v3 > logic_threshold;

%% Build decimal ADC code
codes = double(8*b3 + 4*b2 + 2*b1 + b0);

%% Remove DC for analysis
x = codes - mean(codes);

%% Use MATLAB sinad()
% Returns SINAD in dBc
[SINAD_dB, totalNoiseDist_dB] = sinad(x, fs);

%% Compute ENOB
ENOB = (SINAD_dB - 1.76) / 6.02;

%% Display results
fprintf('Sampling delay used       = %.1f ps\n', sample_delay*1e12);
fprintf('Skipped initial edges     = %d\n', skip_edges);
fprintf('SINAD (MATLAB sinad)      = %.3f dB\n', SINAD_dB);
fprintf('Noise+Dist power          = %.3f dB\n', totalNoiseDist_dB);
fprintf('ENOB                      = %.3f bits\n', ENOB);

%% Plot sampled ADC codes
figure;
plot(t_sample*1e9, codes, 'o-');
xlabel('Time (ns)');
ylabel('ADC code');
title('Sampled 4-bit ADC Output Codes');
grid on;

%% Plot sampled bit voltages
figure;
plot(t_sample*1e9, v3, 'o-', ...
     t_sample*1e9, v2, 'o-', ...
     t_sample*1e9, v1, 'o-', ...
     t_sample*1e9, v0, 'o-');
xlabel('Time (ns)');
ylabel('Sampled bit voltage (V)');
title('Sampled Output Bit Voltages');
legend('out<3>','out<2>','out<1>','out<0>');
grid on;

%% Plot thresholded bits
figure;
stairs(t_sample*1e9, b3, 'LineWidth', 1.2); hold on;
stairs(t_sample*1e9, b2, 'LineWidth', 1.2);
stairs(t_sample*1e9, b1, 'LineWidth', 1.2);
stairs(t_sample*1e9, b0, 'LineWidth', 1.2);
xlabel('Time (ns)');
ylabel('Bit value');
title('Thresholded Output Bits');
legend('b3','b2','b1','b0');
grid on;

%% Optional: let MATLAB plot the SINAD spectrum
figure;
sinad(x, fs);
title('FFT Spectrum of the 4-Bit Flash ADC Output');
set(gcf,'Position',[100 100 700 400]);