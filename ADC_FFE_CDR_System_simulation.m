%% =========================================================
% FOUR REAL 8-GS/s ADCs -> 32-GS/s INTERLEAVED -> 41-TAP RX-FFE
%
% This script uses the ACTUAL Virtuoso outputs of your 4 ADCs.
%
% Architecture:
%   ADC0 @   0.00 ps
%   ADC1 @  31.25 ps
%   ADC2 @  62.50 ps
%   ADC3 @  93.75 ps
%
% Each ADC runs at 8 GHz.
% Interleaved output rate = 4 x 8 GS/s = 32 GS/s.
%
% VERIFIED timing from your ADC ENOB / transfer check:
%   - read final ADC outputs 5 ps after each lane's RISING clock edge
%   - out<3> is MSB, out<0> is LSB
%   - conversion latency = 3 full 8-GHz cycles = 375 ps
%
% The script:
%   1) reads all 16 real ADC output waveforms
%   2) samples each lane at its own rising edge + 5 ps
%   3) removes the 3-cycle ADC latency
%   4) merges all four lanes chronologically -> 32 GS/s
%   5) aligns the real ADC stream with the same PRBS31 reference
%   6) trains the teammate's 41-tap RX-FFE
%   7) reports BER before/after FFE
%   8) measures FFE boost near 16 GHz
%   9) plots actual ADC codes, eye before/after, and FFE taps
%
% IMPORTANT:
% Change ONLY the FILE-NAME section if your exported CSV names differ.
% =========================================================

clear;
clc;
close all;

%% =========================================================
% 0) SYSTEM PARAMETERS
% =========================================================

num_bits = 30000;

data_rate = 32e9;                 % 32 GT/s PAM2
UI = 1/data_rate;                 % 31.25 ps

adc_lane_rate = 8e9;              % each physical ADC lane
adc_lane_period = 1/adc_lane_rate; % 125 ps

num_adc_lanes = 4;
aggregate_sample_rate = 32e9;

adc_resolution_bits = 4;
adc_total_levels = 2^adc_resolution_bits;
adc_mid_code = adc_total_levels/2; % 8

logic_threshold = 0.5;

% VERIFIED from your ADC testing
output_sample_delay = 5e-12;      % read output 5 ps after rising edge
adc_latency_cycles = 3;           % 3 x 125 ps = 375 ps
skip_startup_edges = 10;

% Four interleaved lane phases
lane_phase_ps = [0 31.25 62.50 93.75];
lane_phase = lane_phase_ps*1e-12;

% RX-FFE
num_ffe_taps = 41;
training_data_fraction = 0.80;

% Alignment search
max_symbol_shift = 2000;

% FFE regularization candidates, matching teammate's approach
alpha_list = [ ...
    0, ...
    1e-6, ...
    1e-5, ...
    1e-4, ...
    1e-3, ...
    1e-2, ...
    3e-2, ...
    1e-1, ...
    3e-1, ...
    1, ...
    3, ...
    10, ...
    30, ...
    100];

min_desired_boost_db = 6;
max_desired_boost_db = 12;

max_eye_traces = 600;

rng(1);

%% =========================================================
% 0.1) FILE NAMES
% =========================================================
%
% ADC0 signals:
%   out0, out1, out2, out3
%
% ADC1 signals:
%   out_1<0>, out_1<1>, out_1<2>, out_1<3>
%
% ADC2 signals:
%   out_2<0>, out_2<1>, out_2<2>, out_2<3>
%
% ADC3 signals:
%   out_3<0>, out_3<1>, out_3<2>, out_3<3>
%
% For easy MATLAB use, export/rename them to the names below.
% =========================================================

base_clock_file = 'clk.csv';

adc_files = {
    {'out0.csv',     'out1.csv',     'out2.csv',     'out3.csv'};
    {'out_1_0_.csv',  'out_1_1_.csv',  'out_1_2_.csv',  'out_1_3_.csv'};
    {'out_2_0_.csv',  'out_2_1_.csv',  'out_2_2_.csv',  'out_2_3_.csv'};
    {'out_3_0_.csv',  'out_3_1_.csv',  'out_3_2_.csv',  'out_3_3_.csv'};
};

% Optional files for the analog eye-before-FFE plot.
vinp_file = 'Vin_p.txt';
vinn_file = 'Vin_n.txt';

% Optional original S4P file for channel + FFE frequency-domain plot.
channel_file = 'Cable_BKP_24dB_0p575m_thru1.s4p';

%% =========================================================
% 1) CHECK REQUIRED FILES
% =========================================================

if ~isfile(base_clock_file)
    error('Missing base clock file: %s',base_clock_file);
end

for lane = 1:num_adc_lanes
    for bit = 1:4
        if ~isfile(adc_files{lane}{bit})
            error(['Missing file: %s\n' ...
                   'Edit the FILE NAMES section if your export name is different.'], ...
                   adc_files{lane}{bit});
        end
    end
end

%% =========================================================
% 2) GENERATE SAME PRBS31 REFERENCE AS TEAMMATE'S CODE
% =========================================================

binary_data = generate_prbs31(num_bits);

% PAM2 / NRZ:
% 0 -> -1
% 1 -> +1
symbols_ref = 2*binary_data - 1;

%% =========================================================
% 3) READ BASE 8-GHz CLOCK AND FIND RISING EDGES
% =========================================================

clkData = clean2(readmatrix(base_clock_file,'NumHeaderLines',1));

t_clk = clkData(:,1);
clk   = clkData(:,2);

clk_bin = clk > logic_threshold;

edge_idx = find( ...
    clk_bin(1:end-1)==0 & ...
    clk_bin(2:end)==1);

if isempty(edge_idx)
    error('No rising clock edges found in %s.',base_clock_file);
end

base_rising_edges = ...
    crossing_times(t_clk,clk,edge_idx,logic_threshold);

measured_lane_period = ...
    median(diff(base_rising_edges));

measured_lane_rate = ...
    1/measured_lane_period;

fprintf('\n========================================================\n');
fprintf(' REAL 4-WAY ADC CLOCKING\n');
fprintf('========================================================\n');
fprintf('Measured lane period = %.3f ps\n', ...
    measured_lane_period*1e12);
fprintf('Measured lane rate   = %.6f GHz\n', ...
    measured_lane_rate/1e9);

fprintf('Expected lane phases = ');
fprintf('%.2f ',lane_phase_ps);
fprintf('ps\n');

%% =========================================================
% 4) READ EACH REAL ADC LANE
% =========================================================
%
% We use the verified timing:
%
%   read ADC output at lane rising edge + 5 ps
%
% and assign that code to the analog sample taken 3 full
% 8-GHz cycles earlier:
%
%   t_input = t_lane_rising - 3*Tclk
%
% Then all four lanes are merged using t_input.
% This automatically produces the correct chronological 32-GS/s sequence.
% =========================================================

all_sample_times = [];
all_adc_codes = [];
all_lane_ids = [];

lane_code_cell = cell(num_adc_lanes,1);
lane_input_time_cell = cell(num_adc_lanes,1);

for lane = 1:num_adc_lanes

    fprintf('\nReading ADC lane %d...\n',lane-1);

    % -----------------------------------------------------
    % Read D0,D1,D2,D3 for this lane
    % -----------------------------------------------------

    D0 = clean2(readmatrix(adc_files{lane}{1},'NumHeaderLines',1));
    D1 = clean2(readmatrix(adc_files{lane}{2},'NumHeaderLines',1));
    D2 = clean2(readmatrix(adc_files{lane}{3},'NumHeaderLines',1));
    D3 = clean2(readmatrix(adc_files{lane}{4},'NumHeaderLines',1));

    % -----------------------------------------------------
    % Generate this lane's phase-shifted rising edges.
    %
    % lane 0: base edge
    % lane 1: base edge + 31.25 ps
    % lane 2: base edge + 62.50 ps
    % lane 3: base edge + 93.75 ps
    % -----------------------------------------------------

    lane_edges = ...
        base_rising_edges + lane_phase(lane);

    % Read final digital outputs 5 ps after this edge.
    t_read = ...
        lane_edges + output_sample_delay;

    % Output code corresponds to input sampled 3 cycles earlier.
    t_input_sample = ...
        lane_edges - ...
        adc_latency_cycles*measured_lane_period;

    % -----------------------------------------------------
    % Keep only valid waveform range
    % -----------------------------------------------------

    tmax_lane = min([ ...
        D0(end,1), ...
        D1(end,1), ...
        D2(end,1), ...
        D3(end,1)]);

    valid = ...
        t_read >= max([D0(1,1),D1(1,1),D2(1,1),D3(1,1)]) & ...
        t_read <= tmax_lane & ...
        t_input_sample >= 0;

    t_read = t_read(valid);
    t_input_sample = t_input_sample(valid);

    % Ignore startup edges after valid-range filtering.
    if length(t_read) <= skip_startup_edges
        error('Not enough valid conversions in ADC lane %d.',lane-1);
    end

    t_read = t_read(skip_startup_edges+1:end);
    t_input_sample = t_input_sample(skip_startup_edges+1:end);

    % -----------------------------------------------------
    % Read bit voltages at rising edge + 5 ps
    % -----------------------------------------------------

    v0 = interp1(D0(:,1),D0(:,2),t_read,'linear');
    v1 = interp1(D1(:,1),D1(:,2),t_read,'linear');
    v2 = interp1(D2(:,1),D2(:,2),t_read,'linear');
    v3 = interp1(D3(:,1),D3(:,2),t_read,'linear');

    good = ...
        isfinite(v0) & ...
        isfinite(v1) & ...
        isfinite(v2) & ...
        isfinite(v3);

    v0 = v0(good);
    v1 = v1(good);
    v2 = v2(good);
    v3 = v3(good);

    t_input_sample = ...
        t_input_sample(good);

    % Logic threshold
    b0 = v0 > logic_threshold;
    b1 = v1 > logic_threshold;
    b2 = v2 > logic_threshold;
    b3 = v3 > logic_threshold;

    % VERIFIED mapping:
    % out<3> = MSB, out<0> = LSB
    lane_codes = double( ...
        8*b3 + ...
        4*b2 + ...
        2*b1 + ...
        b0);

    lane_codes = lane_codes(:);
    t_input_sample = t_input_sample(:);

    lane_code_cell{lane} = lane_codes;
    lane_input_time_cell{lane} = t_input_sample;

    fprintf('  valid conversions = %d\n',length(lane_codes));
    fprintf('  code range        = %d to %d\n', ...
        min(lane_codes),max(lane_codes));

    % Add lane data to aggregate list
    all_sample_times = ...
        [all_sample_times; t_input_sample];

    all_adc_codes = ...
        [all_adc_codes; lane_codes];

    all_lane_ids = ...
        [all_lane_ids; ...
         (lane-1)*ones(length(lane_codes),1)];
end

%% =========================================================
% 5) MERGE / INTERLEAVE THE 4 ADC LANES
% =========================================================

[adc32_sample_times,sort_idx] = ...
    sort(all_sample_times);

adc32_codes = ...
    all_adc_codes(sort_idx);

adc32_lane = ...
    all_lane_ids(sort_idx);

% Remove any accidental duplicate times.
% Duplicates would indicate incorrect phase generation.
dt = diff(adc32_sample_times);

if any(dt <= 0)
    [adc32_sample_times,unique_idx] = ...
        unique(adc32_sample_times,'stable');

    adc32_codes = adc32_codes(unique_idx);
    adc32_lane = adc32_lane(unique_idx);
end

% Measure aggregate sampling interval/rate.
dt = diff(adc32_sample_times);

median_aggregate_period = ...
    median(dt);

measured_aggregate_rate = ...
    1/median_aggregate_period;

fprintf('\n========================================================\n');
fprintf(' INTERLEAVED REAL ADC STREAM\n');
fprintf('========================================================\n');
fprintf('Total interleaved samples = %d\n', ...
    length(adc32_codes));
fprintf('Median sample spacing      = %.3f ps\n', ...
    median_aggregate_period*1e12);
fprintf('Measured aggregate rate    = %.3f GS/s\n', ...
    measured_aggregate_rate/1e9);
fprintf('ADC code range             = %d to %d\n', ...
    min(adc32_codes),max(adc32_codes));

if abs(measured_aggregate_rate-aggregate_sample_rate) / ...
        aggregate_sample_rate > 0.02

    warning(['Aggregate ADC sample rate is not close to 32 GS/s. ' ...
             'Check the four clock phases.']);
end

% Save the real 32-GS/s ADC stream for later use.
writematrix( ...
    [adc32_sample_times adc32_codes adc32_lane], ...
    'real_adc32_interleaved.csv');

%% =========================================================
% 5.1) PLOT INTERLEAVED ADC CODES
% =========================================================

figure('Name','Real 4-Way Interleaved ADC Output','Color','w');

Nplot = min(250,length(adc32_codes));

stairs(1:Nplot, ...
       adc32_codes(1:Nplot), ...
       'LineWidth',1.2);

xlabel('32-GS/s Sample / Symbol Index');
ylabel('ADC Code');
title('Actual 4-bit ADC Output After 4-Way Interleaving');

ylim([-1 16]);
yticks(0:15);
grid on;

%% =========================================================
% 5.2) VERIFY LANE ORDER / TIMING
% =========================================================

figure('Name','ADC Interleaving Check','Color','w');

Ntim = min(80,length(adc32_codes));

stem( ...
    adc32_sample_times(1:Ntim)*1e9, ...
    adc32_lane(1:Ntim), ...
    'filled');

xlabel('Equivalent Analog Sampling Time [ns]');
ylabel('ADC Lane');
yticks(0:3);
title('Four ADC Lanes in Chronological Sampling Order');
grid on;

%% =========================================================
% 6) CENTER REAL ADC CODES FOR DIGITAL FFE
% =========================================================

ffe_input_full = ...
    double(adc32_codes) - adc_mid_code;

ffe_input_full = ...
    ffe_input_full(:);

%% =========================================================
% 7) TIME-BASED ALIGNMENT -- ALLOW NEGATIVE PRBS SHIFT
% =========================================================
%
% The channel-delayed Vin_p/Vin_n files are not aligned to PRBS
% symbol 1 at t = 0.  Searching only for a POSITIVE reference
% shift (as in the old script) gives a false match.
%
% Here the PRBS index is derived from the ACTUAL analog sampling
% time.  Search for the channel/PWL delay (including NEGATIVE
% symbol shifts) using ONLY the early 10-30 ns calibration interval.
% Then use one fixed alignment for the entire record.
%
% No new Virtuoso run is required.
% =========================================================

if ~isfile(vinp_file) || ~isfile(vinn_file)
    error(['This corrected reference-alignment step requires ' ...
           'Vin_p.txt and Vin_n.txt in the MATLAB folder.']);
end

Vp_align = clean2(readmatrix(vinp_file));
Vn_align = clean2(readmatrix(vinn_file));

vp_at_sample = interp1(Vp_align(:,1),Vp_align(:,2), ...
                       adc32_sample_times,'linear',NaN);
vn_at_sample = interp1(Vn_align(:,1),Vn_align(:,2), ...
                       adc32_sample_times,'linear',NaN);
vin_diff_at_sample = vp_at_sample - vn_at_sample;

% A 0-based symbol index from physical source time, not ADC array
% position.  round() is appropriate for clock edges nominally at
% integer multiples of the 31.25-ps UI.
physical_symbol = round(adc32_sample_times/UI);

% The first few ns contain the channel/PRBS-startup transient.
calibration_mask = adc32_sample_times >= 10e-9 & ...
                   adc32_sample_times < 30e-9 & ...
                   isfinite(vin_diff_at_sample);

if nnz(calibration_mask) < 200
    error('Insufficient analog samples in the 10-30 ns calibration interval.');
end

best_corr = -inf;
best_shift = NaN;
best_polarity = 1;

for signed_shift = -350:350
    ref_index = physical_symbol + signed_shift + 1; % MATLAB 1-based

    good = calibration_mask & ...
           ref_index >= 1 & ref_index <= numel(symbols_ref);

    if nnz(good) < 200
        continue;
    end

    analog_segment = vin_diff_at_sample(good);
    reference_segment = symbols_ref(ref_index(good));

    xa = analog_segment(:) - mean(analog_segment(:));
    rr = reference_segment(:) - mean(reference_segment(:));
    rho = sum(xa.*rr)/(sqrt(sum(xa.^2)*sum(rr.^2))+eps);

    if abs(rho) > best_corr
        best_corr = abs(rho);
        best_shift = signed_shift;
        best_polarity = sign(rho);
    end
end

if isnan(best_shift)
    error('Could not find the analog waveform / PRBS alignment.');
end

% Apply the fixed shift to EVERY ADC sample.  Discard pre-channel
% startup and samples outside the reference length.
reference_index = physical_symbol + best_shift + 1;
keep = adc32_sample_times >= 10e-9 & ...
       reference_index >= 1 & ...
       reference_index <= numel(symbols_ref) & ...
       isfinite(vin_diff_at_sample);

ffe_input_signal = best_polarity * ffe_input_full(keep);
ref_bits = symbols_ref(reference_index(keep));
aligned_sample_times = adc32_sample_times(keep);

ffe_input_signal = ffe_input_signal(:);
ref_bits = ref_bits(:);
num_avail = numel(ref_bits);

if num_avail < 300
    error('Too few valid, time-aligned samples for FFE training.');
end

corr_adc = corrcoef(ffe_input_signal, ref_bits);
corr_analog = corrcoef(best_polarity*vin_diff_at_sample(keep), ref_bits);

fprintf('\n========================================================\n');
fprintf(' CORRECTED PHYSICAL-TIME / PRBS ALIGNMENT\n');
fprintf('========================================================\n');
fprintf('Signed PRBS shift       = %+d symbols\n', best_shift);
fprintf('Approx. source delay   = %.3f ns\n', -best_shift*UI*1e9);
fprintf('Calibration correlation = %.4f\n', best_corr);
fprintf('Analog / PRBS corr.     = %.4f\n', corr_analog(1,2));
fprintf('ADC / PRBS corr.        = %.4f\n', corr_adc(1,2));
fprintf('Usable samples          = %d\n', num_avail);
fprintf('========================================================\n');

% Correlation versus ABSOLUTE input-file time.
for ns = 10:10:80
    segment = aligned_sample_times >= ns*1e-9 & ...
              aligned_sample_times < (ns+10)*1e-9;
    if nnz(segment) >= 50
        ca = corrcoef(ffe_input_signal(segment), ref_bits(segment));
        fprintf('%.0f-%.0f ns ADC/PRBS corr = %.4f\n', ...
                ns,ns+10,ca(1,2));
    end
end

%% =========================================================
% 8) BUILD 41-TAP RX-FFE MATRIX
% =========================================================

ffe_center_tap = ...
    ceil(num_ffe_taps/2);

padded_in = ...
    [zeros(ffe_center_tap-1,1); ...
     ffe_input_signal; ...
     zeros(num_ffe_taps-ffe_center_tap,1)];

U_matrix = ...
    zeros(num_avail,num_ffe_taps);

for i = 1:num_avail

    U_matrix(i,:) = ...
        padded_in(i:i+num_ffe_taps-1).';
end

num_train = ...
    floor(training_data_fraction*num_avail);

% Ensure reasonable validation set with a 100-ns simulation.
if num_avail-num_train < 300
    num_train = floor(0.70*num_avail);
end

if num_train <= num_ffe_taps
    error('Not enough samples to train the 41-tap FFE.');
end

U_train = ...
    U_matrix(1:num_train,:);

d_train = ...
    ref_bits(1:num_train);

val_idx = ...
    (num_train+1):num_avail;

%% =========================================================
% 9) REGULARIZED LEAST-SQUARES FFE TRAINING
%     SAME 16-GHz BOOST IDEA AS TEAMMATE'S CODE
% =========================================================

R = ...
    U_train.'*U_train;

regularization_scale = ...
    trace(R)/num_ffe_taps;

best_score = -inf;
best_weights = [];
best_output = [];
best_alpha = NaN;
best_lambda = NaN;
best_eye_opening_pct = -inf;
best_boost_db = NaN;
best_ber = inf;

ffe_sample_rate = aggregate_sample_rate;
data_nyquist_hz = data_rate/2; % 16 GHz

fprintf('\n========================================================\n');
fprintf(' REAL ADC -> RX-FFE REGULARIZATION SEARCH\n');
fprintf('========================================================\n');

for ai = 1:length(alpha_list)

    alpha = alpha_list(ai);

    lambda_reg = ...
        alpha*regularization_scale;

    % Regularized least-squares FFE
    w_try = ...
        (R + lambda_reg*eye(num_ffe_taps)) \ ...
        (U_train.'*d_train);

    y_try = ...
        U_matrix*w_try;

    % Correct possible FFE polarity inversion.
    if mean(y_try(ref_bits>0)) < ...
       mean(y_try(ref_bits<0))

        y_try = -y_try;
        w_try = -w_try;
    end

    % Normalize equalized level separation to ~2.
    y_train = ...
        y_try(1:num_train);

    r_train = ...
        ref_bits(1:num_train);

    pos_train = ...
        y_train(r_train>0);

    neg_train = ...
        y_train(r_train<0);

    level_sep = ...
        mean(pos_train)-mean(neg_train);

    if level_sep > 0
        gain = 2/level_sep;
    else
        gain = 1;
    end

    y_try = y_try*gain;
    w_try = w_try*gain;

    % Validation performance
    y_val = ...
        y_try(val_idx);

    ref_val = ...
        ref_bits(val_idx);

    decisions_val = ...
        y_val >= 0;

    bits_val = ...
        ref_val > 0;

    ber_val = ...
        sum(decisions_val ~= bits_val) / ...
        length(bits_val);

    pos_val = ...
        y_val(ref_val>0);

    neg_val = ...
        y_val(ref_val<0);

    p5_pos = ...
        local_percentile(pos_val,5);

    p95_neg = ...
        local_percentile(neg_val,95);

    eye_opening = ...
        p5_pos-p95_neg;

    eye_opening_pct = ...
        100*eye_opening/2;

    % FFE boost at 16 GHz
    [H_try,f_try] = ...
        freqz(w_try,1,8192,ffe_sample_rate);

    H_try_db = ...
        20*log10(abs(H_try)+eps);

    H_try_relative_db = ...
        H_try_db-H_try_db(1);

    [~,idx16] = ...
        min(abs(f_try-data_nyquist_hz));

    boost_16_db = ...
        H_try_relative_db(idx16);

    fprintf(['alpha = %-8.1e | ' ...
             'BER = %.3e | ' ...
             'eye = %7.2f %% | ' ...
             'FFE boost@16G = %7.2f dB\n'], ...
             alpha, ...
             ber_val, ...
             eye_opening_pct, ...
             boost_16_db);

    boost_ok = ...
        boost_16_db >= min_desired_boost_db && ...
        boost_16_db <= max_desired_boost_db;

    if boost_ok

        score = ...
            eye_opening_pct - ...
            1000*ber_val;

        if score > best_score

            best_score = score;
            best_weights = w_try;
            best_output = y_try;
            best_alpha = alpha;
            best_lambda = lambda_reg;
            best_eye_opening_pct = eye_opening_pct;
            best_boost_db = boost_16_db;
            best_ber = ber_val;
        end
    end
end

%% =========================================================
% 9.1) FALLBACK IF NO 6-12 dB CANDIDATE
% =========================================================

if isempty(best_weights)

    warning(['No candidate produced 6-12 dB boost at 16 GHz. ' ...
             'Selecting the best BER/eye solution instead.']);

    best_score = -inf;

    for ai = 1:length(alpha_list)

        alpha = alpha_list(ai);

        lambda_reg = ...
            alpha*regularization_scale;

        w_try = ...
            (R + lambda_reg*eye(num_ffe_taps)) \ ...
            (U_train.'*d_train);

        y_try = ...
            U_matrix*w_try;

        if mean(y_try(ref_bits>0)) < ...
           mean(y_try(ref_bits<0))

            y_try = -y_try;
            w_try = -w_try;
        end

        y_train = y_try(1:num_train);
        r_train = ref_bits(1:num_train);

        level_sep = ...
            mean(y_train(r_train>0)) - ...
            mean(y_train(r_train<0));

        if level_sep > 0
            gain = 2/level_sep;
        else
            gain = 1;
        end

        y_try = y_try*gain;
        w_try = w_try*gain;

        y_val = y_try(val_idx);
        ref_val = ref_bits(val_idx);

        ber_val = ...
            sum((y_val>=0) ~= (ref_val>0)) / ...
            length(ref_val);

        pos_val = y_val(ref_val>0);
        neg_val = y_val(ref_val<0);

        eye_opening_pct = ...
            100 * ...
            (local_percentile(pos_val,5) - ...
             local_percentile(neg_val,95)) / 2;

        score = ...
            eye_opening_pct - ...
            1000*ber_val;

        [H_try,f_try] = ...
            freqz(w_try,1,8192,ffe_sample_rate);

        H_try_db = ...
            20*log10(abs(H_try)+eps);

        H_try_relative_db = ...
            H_try_db-H_try_db(1);

        [~,idx16] = ...
            min(abs(f_try-data_nyquist_hz));

        boost_16_db = ...
            H_try_relative_db(idx16);

        if score > best_score

            best_score = score;
            best_weights = w_try;
            best_output = y_try;
            best_alpha = alpha;
            best_lambda = lambda_reg;
            best_eye_opening_pct = eye_opening_pct;
            best_boost_db = boost_16_db;
            best_ber = ber_val;
        end
    end
end

weights = ...
    best_weights;

equalized_output = ...
    best_output;

fprintf('\n========================================================\n');
fprintf(' SELECTED REAL-ADC RX-FFE\n');
fprintf('========================================================\n');
fprintf('Selected alpha          = %.3e\n',best_alpha);
fprintf('Selected lambda         = %.3e\n',best_lambda);
fprintf('Validation BER          = %.3e\n',best_ber);
fprintf('Validation eye opening  = %.2f %%\n',best_eye_opening_pct);
fprintf('FFE boost near 16 GHz   = %.2f dB\n',best_boost_db);

%% =========================================================
% 10) FIXED FFE DECISION TIMING
% =========================================================
% Do not optimize a decision delay on the full validation/test record.
% The centered FFE is trained against the correctly aligned reference.
best_decision_delay = 0;
fprintf('FFE decision delay = %d symbols (fixed)\n',best_decision_delay);

%% =========================================================
% 11) BER BEFORE AND AFTER REAL ADC + FFE
% =========================================================

reference_bits = ...
    ref_bits > 0;

decided_before = ...
    ffe_input_signal >= 0;

errors_before = ...
    sum(decided_before ~= reference_bits);

ber_before_ffe = ...
    errors_before/length(reference_bits);

decided_after = ...
    equalized_output >= 0;

errors_after = ...
    sum(decided_after ~= reference_bits);

ber_after_ffe = ...
    errors_after/length(reference_bits);

fprintf('\n========================================================\n');
fprintf(' REAL 4-ADC + FFE BER RESULTS\n');
fprintf('========================================================\n');
fprintf('BER before RX-FFE = %.3e (%d errors)\n', ...
    ber_before_ffe,errors_before);
fprintf('BER after RX-FFE  = %.3e (%d errors)\n', ...
    ber_after_ffe,errors_after);

%% =========================================================
% 12) DIGITAL EYE OPENING AFTER FFE
% =========================================================

pos_samples = ...
    equalized_output(ref_bits>0);

neg_samples = ...
    equalized_output(ref_bits<0);

p5_pos = ...
    local_percentile(pos_samples,5);

p95_neg = ...
    local_percentile(neg_samples,95);

actual_eye_opening = ...
    p5_pos-p95_neg;

actual_eye_opening_pct = ...
    100*actual_eye_opening/2;

fprintf('Digital eye opening after FFE = %.2f %%\n', ...
    actual_eye_opening_pct);

%% =========================================================
% 13) FFE FREQUENCY RESPONSE AT 32 GS/s
% =========================================================

[Hffe,fffe] = ...
    freqz(weights,1,8192,aggregate_sample_rate);

Hffe_db = ...
    20*log10(abs(Hffe)+eps);

Hffe_relative_db = ...
    Hffe_db-Hffe_db(1);

[~,idx_ffe_16] = ...
    min(abs(fffe-data_nyquist_hz));

ffe_boost_16_db = ...
    Hffe_relative_db(idx_ffe_16);

figure('Name','Real ADC RX-FFE Frequency Response','Color','w');

plot( ...
    fffe/1e9, ...
    Hffe_relative_db, ...
    'LineWidth',1.5);

hold on;

plot( ...
    fffe(idx_ffe_16)/1e9, ...
    ffe_boost_16_db, ...
    'o', ...
    'MarkerSize',8, ...
    'LineWidth',1.5);

xline(16,'--','16 GHz');

xlabel('Frequency [GHz]');
ylabel('Relative FFE Gain [dB]');
title('41-Tap RX-FFE Response Using Actual 4-ADC Data');

xlim([0 16]);
grid on;

%% =========================================================
% 14) FFE TAP COEFFICIENTS
% =========================================================

figure('Name','Real ADC RX-FFE Tap Coefficients','Color','w');

stem( ...
    1:num_ffe_taps, ...
    weights, ...
    'filled');

xlabel('Tap Number');
ylabel('Coefficient');
title('41-Tap RX-FFE Coefficients Using Actual ADC Outputs');
grid on;

%% =========================================================
% 15) DIGITAL EYE BEFORE / AFTER FFE
% =========================================================
%
% Because the interleaved ADC now gives one sample per UI,
% a symbol-rate digital eye before/after FFE is valid.
% =========================================================

figure( ...
    'Name','Digital Eye Before and After Real ADC RX-FFE', ...
    'Color','w', ...
    'Position',[100 100 900 700]);

% ---------------------------------------------------------
% BEFORE FFE: actual interleaved ADC samples
% ---------------------------------------------------------

subplot(2,1,1);
hold on;

num_before_traces = ...
    min(max_eye_traces,length(ffe_input_signal)-2);

start_before = ...
    max(1, ...
        floor(0.65*length(ffe_input_signal)) - ...
        num_before_traces);

end_before = ...
    min( ...
        start_before+num_before_traces, ...
        length(ffe_input_signal)-2);

for k = start_before:end_before

    plot( ...
        [-1 0 1], ...
        [ffe_input_signal(k), ...
         ffe_input_signal(k+1), ...
         ffe_input_signal(k+2)], ...
        'LineWidth',0.5);
end

yline(0,'--k');

title('Digital Symbol Eye BEFORE RX-FFE - Actual 4 ADCs');
xlabel('Symbol Time [UI]');
ylabel('ADC Code - 8');
grid on;

% ---------------------------------------------------------
% AFTER FFE
% ---------------------------------------------------------

subplot(2,1,2);
hold on;

num_after_traces = ...
    min(max_eye_traces,length(equalized_output)-2);

start_after = ...
    max(1, ...
        floor(0.65*length(equalized_output)) - ...
        num_after_traces);

end_after = ...
    min( ...
        start_after+num_after_traces, ...
        length(equalized_output)-2);

for k = start_after:end_after

    plot( ...
        [-1 0 1], ...
        [equalized_output(k), ...
         equalized_output(k+1), ...
         equalized_output(k+2)], ...
        'LineWidth',0.5);
end

yline(0,'--k');

title('Digital Symbol Eye AFTER RX-FFE - Actual 4 ADCs');
xlabel('Symbol Time [UI]');
ylabel('Equalized Symbol Level');
grid on;

%% =========================================================
% 16) OPTIONAL: ANALOG EYE BEFORE + DIGITAL EYE AFTER
%     SAME STYLE AS TEAMMATE'S ORIGINAL FIGURE
% =========================================================

if isfile(vinp_file) && isfile(vinn_file)

    VinP = clean2(readmatrix(vinp_file));
    VinN = clean2(readmatrix(vinn_file));

    t_analog = VinP(:,1);

    VinN_interp = ...
        interp1( ...
            VinN(:,1), ...
            VinN(:,2), ...
            t_analog, ...
            'linear', ...
            'extrap');

    vdiff_analog = ...
        VinP(:,2)-VinN_interp;

    peak_v = ...
        max(abs(vdiff_analog));

    if peak_v > 0
        vdiff_analog = ...
            vdiff_analog/peak_v;
    end

    % Original PWL waveform was generated at 512 GS/s:
    source_sample_rate = ...
        data_rate*16;

    samples_per_UI = ...
        round(source_sample_rate/data_rate); % 16

    samples_per_eye = ...
        2*samples_per_UI;

    % Start after initial transient.
    start_index = ...
        max(1,round(10e-9*source_sample_rate));

    num_analog_traces = ...
        min( ...
            max_eye_traces, ...
            floor((length(vdiff_analog)- ...
                   start_index- ...
                   samples_per_eye) / ...
                  samples_per_UI));

    eye_before = ...
        zeros(samples_per_eye+1,num_analog_traces);

    idx = start_index;

    for k = 1:num_analog_traces

        if idx+samples_per_eye <= length(vdiff_analog)

            eye_before(:,k) = ...
                vdiff_analog(idx:idx+samples_per_eye);
        end

        idx = idx+samples_per_UI;
    end

    time_axis_before = ...
        linspace(-1,1,samples_per_eye+1);

    figure( ...
        'Name','Analog Eye Before / Digital Eye After Real ADC FFE', ...
        'Color','w', ...
        'Position',[100 100 900 700]);

    subplot(2,1,1);

    plot( ...
        time_axis_before, ...
        eye_before, ...
        'LineWidth',0.5);

    title('Eye BEFORE Digital RX-FFE');
    xlabel('Time [UI]');
    ylabel('Normalized Amplitude');
    xlim([-1 1]);
    ylim([-1.2 1.2]);
    grid on;

    subplot(2,1,2);
    hold on;

    for k = start_after:end_after

        plot( ...
            [-1 0 1], ...
            [equalized_output(k), ...
             equalized_output(k+1), ...
             equalized_output(k+2)], ...
            'LineWidth',0.5);
    end

    yline(0,'--k');

    title('Digital Symbol Eye AFTER Real 4-ADC RX-FFE');
    xlabel('Symbol Time [UI]');
    ylabel('Equalized Symbol Level');
    xlim([-1 1]);
    grid on;
end

%% =========================================================
% 17) OPTIONAL: ORIGINAL CHANNEL + REAL-ADC FFE RESPONSE
% =========================================================

if isfile(channel_file)

    try

        S = sparameters(channel_file);

        f_sparam = S.Frequencies;
        Smat = S.Parameters;

        S21 = squeeze(Smat(2,1,:));
        S23 = squeeze(Smat(2,3,:));
        S41 = squeeze(Smat(4,1,:));
        S43 = squeeze(Smat(4,3,:));

        Sdd21 = ...
            0.5*(S21-S23-S41+S43);

        Sdd21_db = ...
            20*log10(abs(Sdd21)+eps);

        channel_abs_db = ...
            interp1( ...
                f_sparam, ...
                Sdd21_db, ...
                fffe, ...
                'linear', ...
                'extrap');

        combined_abs_db = ...
            channel_abs_db + ...
            Hffe_relative_db;

        figure( ...
            'Name','Channel + Real ADC RX-FFE', ...
            'Color','w');

        plot( ...
            fffe/1e9, ...
            channel_abs_db, ...
            'LineWidth',1.5);

        hold on;

        plot( ...
            fffe/1e9, ...
            Hffe_relative_db, ...
            'LineWidth',1.5);

        plot( ...
            fffe/1e9, ...
            combined_abs_db, ...
            'LineWidth',1.5);

        xline(16,'--','16 GHz');
        yline(0,'--');

        xlabel('Frequency [GHz]');
        ylabel('Magnitude [dB]');

        title('Measured Channel + RX-FFE Using Actual ADC Data');

        legend( ...
            'Measured Residual Channel', ...
            'Digital RX-FFE', ...
            'Combined', ...
            'Location','best');

        xlim([0 16]);
        grid on;

    catch ME

        warning( ...
            'Could not generate optional channel plot: %s', ...
            ME.message);
    end
end

%% =========================================================
% 18) FINAL SUMMARY
% =========================================================

fprintf('\n========================================================\n');
fprintf(' FINAL REAL 4-ADC + RX-FFE SUMMARY\n');
fprintf('========================================================\n');
fprintf('Data rate                 = %.1f GT/s\n',data_rate/1e9);
fprintf('ADC lanes                 = %d\n',num_adc_lanes);
fprintf('ADC rate per lane         = %.1f GS/s\n',adc_lane_rate/1e9);
fprintf('Interleaved ADC rate      = %.3f GS/s\n', ...
    measured_aggregate_rate/1e9);
fprintf('ADC resolution            = %d bits\n',adc_resolution_bits);
fprintf('ADC latency               = %d cycles = %.1f ps\n', ...
    adc_latency_cycles, ...
    adc_latency_cycles*adc_lane_period*1e12);
fprintf('PRBS alignment corr.      = %.4f\n',best_corr);
fprintf('RX-FFE taps               = %d\n',num_ffe_taps);
fprintf('RX-FFE boost near 16 GHz  = %.2f dB\n',ffe_boost_16_db);
fprintf('BER before RX-FFE         = %.3e\n',ber_before_ffe);
fprintf('BER after RX-FFE          = %.3e\n',ber_after_ffe);
fprintf('Digital eye after FFE     = %.2f %%\n', ...
    actual_eye_opening_pct);
fprintf('Selected alpha            = %.3e\n',best_alpha);
fprintf('Selected lambda           = %.3e\n',best_lambda);
fprintf('========================================================\n');

%% =========================================================
% LOCAL FUNCTION: PRBS31
% =========================================================

function bits = generate_prbs31(N)

    shift_register = ones(31,1);
    bits = zeros(N,1);

    for i = 1:N

        bits(i) = shift_register(end);

        feedback = ...
            xor( ...
                shift_register(end), ...
                shift_register(28));

        shift_register = ...
            [feedback; ...
             shift_register(1:end-1)];
    end
end

%% =========================================================
% LOCAL FUNCTION: PERCENTILE
% =========================================================

function v = local_percentile(x,p)

    x = sort(x(:));
    n = length(x);

    if n == 0
        v = NaN;
        return;
    end

    pos = ...
        1 + (p/100)*(n-1);

    lo = floor(pos);
    hi = ceil(pos);

    if lo == hi
        v = x(lo);
    else
        v = ...
            x(lo) + ...
            (pos-lo)*(x(hi)-x(lo));
    end
end

%% =========================================================
% LOCAL FUNCTION: CLEAN TWO-COLUMN FILE
% =========================================================

function A = clean2(A)

    if size(A,2) < 2
        error('Input file does not contain two numeric columns.');
    end

    A = A(:,1:2);
    A = A(all(isfinite(A),2),:);
    A = sortrows(A,1);
end

%% =========================================================
% LOCAL FUNCTION: CLOCK THRESHOLD CROSSING
% =========================================================

function tc = crossing_times(t,v,idx,thr)

    tc = zeros(length(idx),1);

    for k = 1:length(idx)

        i = idx(k);

        t1 = t(i);
        t2 = t(i+1);

        v1 = v(i);
        v2 = v(i+1);

        if v2 == v1
            tc(k) = t2;
        else
            tc(k) = ...
                t1 + ...
                (thr-v1)*(t2-t1)/(v2-v1);
        end
    end
end
