%% =========================================================
% REAL S-PARAMETER RX-FFE VALIDATION
%
% RX-FFE validation flow:
%
% PRBS31
%   -> TX pulse shaping
%   -> REAL S4P residual channel (~15-17 dB at 16 GHz)
%   -> ideal 4-bit ADC
%   -> 41-tap digital RX-FFE
%
% Channel model:
% The S-parameter channel represents the RESIDUAL channel
% after the assumed CTLE compensation.
%
% 
%
% Main objectives:
% 1) Verify the real channel loss at 16 GHz
% 2) Generate degraded eye before RX-FFE
% 3) Train digital RX-FFE
% 4) Measure RX-FFE equalization boost in dB
% 5) Compare BER and eye before/after RX-FFE
%% =========================================================

clear;
clc;
close all;

%% =========================================================
% 0) SYSTEM PARAMETERS
% =========================================================

num_bits = 30000;

% PCIe Gen5-like PAM2 / NRZ signaling
data_rate = 32e9;                  % 32 GT/s
samples_per_symbol = 16;
sample_rate = data_rate * samples_per_symbol;

data_nyquist_hz = data_rate / 2;   % 16 GHz

% A high SNR is used so the channel ISI remains the main impairment.
% This makes it easier to evaluate the equalization performed by the RX-FFE.
signal_to_noise_ratio_db = 60;

% ADC
adc_resolution_bits = 4;
adc_total_levels = 2^adc_resolution_bits;
adc_mid_threshold = adc_total_levels / 2;

% RX-FFE
num_ffe_taps = 41;
training_data_fraction = 0.8;

% Plot parameters
max_eye_traces = 600;

% Real channel file
channel_file = 'Cable_BKP_24dB_0p575m_thru1.s4p';

rng(1);


%% =========================================================
% 1) PRBS31 DATA GENERATION
% =========================================================

binary_data = generate_prbs31(num_bits);

% PAM2 / NRZ mapping:
% bit 0 -> -1
% bit 1 -> +1
symbols_ref = 2 * binary_data - 1;

symbols_tx = symbols_ref;

bit_plot_start = 500;
bit_plot_end = bit_plot_start + 79;

figure('Name','PRBS31 Input','Color','w');

subplot(2,1,1);
stairs(bit_plot_start:bit_plot_end, ...
       binary_data(bit_plot_start:bit_plot_end), ...
       'LineWidth',1.2);

title('PRBS31 Binary Data');
xlabel('Bit Index');
ylabel('Bit');
ylim([-0.2 1.2]);
grid on;

subplot(2,1,2);
stairs(bit_plot_start:bit_plot_end, ...
       symbols_ref(bit_plot_start:bit_plot_end), ...
       'LineWidth',1.2);

title('PAM2 / NRZ Symbols');
xlabel('Symbol Index');
ylabel('Symbol Level');
ylim([-1.2 1.2]);
grid on;


%% =========================================================
% 2) TX PULSE SHAPING
% =========================================================

rrc_rolloff_factor = 0.35;
rrc_filter_span = 10;

rrc_impulse_response = ...
    rcosdesign(rrc_rolloff_factor, ...
               rrc_filter_span, ...
               samples_per_symbol, ...
               'sqrt');

tx_waveform = upfirdn(symbols_tx, ...
                      rrc_impulse_response, ...
                      samples_per_symbol);

tx_waveform = tx_waveform / max(abs(tx_waveform));

figure('Name','TX Waveform','Color','w');

plot(tx_waveform(1200:1700),'LineWidth',1);

title('Transmitted PAM2 Waveform');
xlabel('Sample Index');
ylabel('Normalized Amplitude');
ylim([-1.2 1.2]);
grid on;


%% =========================================================
% 3) READ REAL 4-PORT S-PARAMETER CHANNEL
% =========================================================

fprintf('\n========================================================\n');
fprintf('           REAL S-PARAMETER CHANNEL\n');
fprintf('========================================================\n');

S = sparameters(channel_file);

f_sparam = S.Frequencies;
Smat = S.Parameters;

% ---------------------------------------------------------
% PORT MAPPING USED FOR THIS S4P FILE
%
% Differential TX pair: ports 1 and 3
% Differential RX pair: ports 2 and 4
%
% Differential through transfer:
%
% Sdd21 = 0.5 * (S21 - S23 - S41 + S43)
% ---------------------------------------------------------

S21 = squeeze(Smat(2,1,:));
S23 = squeeze(Smat(2,3,:));

S41 = squeeze(Smat(4,1,:));
S43 = squeeze(Smat(4,3,:));

Sdd21 = 0.5 * (S21 - S23 - S41 + S43);

Sdd21_db = 20*log10(abs(Sdd21) + eps);


%% =========================================================
% 3.1 CHECK LOSS AT 16 GHz
% =========================================================

[~, idx16] = min(abs(f_sparam - data_nyquist_hz));

actual_frequency_16 = f_sparam(idx16);
channel_gain_16_db = Sdd21_db(idx16);
channel_loss_16_db = -channel_gain_16_db;

fprintf('Frequency checked              = %.3f GHz\n', ...
        actual_frequency_16/1e9);

fprintf('Differential Sdd21 at 16 GHz   = %.2f dB\n', ...
        channel_gain_16_db);

fprintf('Channel insertion loss         = %.2f dB\n', ...
        channel_loss_16_db);


%% =========================================================
% 3.2 PLOT REAL CHANNEL RESPONSE
% =========================================================

figure('Name','Real S-Parameter Channel','Color','w');

plot(f_sparam/1e9, ...
     Sdd21_db, ...
     'LineWidth',1.5);

hold on;

plot(actual_frequency_16/1e9, ...
     channel_gain_16_db, ...
     'o', ...
     'MarkerSize',8, ...
     'LineWidth',1.5);

xline(16,'--','16 GHz');

title('Measured Differential Channel Response');
xlabel('Frequency [GHz]');
ylabel('|S_{dd21}| [dB]');

xlim([0 32]);

grid on;


%% =========================================================
% 4) APPLY REAL S-PARAMETER CHANNEL TO TX WAVEFORM
% =========================================================

% Apply both components of the measured differential response:
%   magnitude of Sdd21
%   phase of Sdd21
%
% Using both magnitude and phase includes the measured channel loss, delay,
% and dispersion in the time-domain waveform.

signal_length = length(tx_waveform);

% Zero padding prevents circular wrap-around.
Nfft_channel = 2^nextpow2(2 * signal_length);

positive_frequency_axis = ...
    (0:Nfft_channel/2).' * sample_rate/Nfft_channel;


% ---------------------------------------------------------
% Interpolate the real and imaginary parts separately.
% ---------------------------------------------------------

H_real = interp1(f_sparam, ...
                 real(Sdd21), ...
                 positive_frequency_axis, ...
                 'pchip', ...
                 0);

H_imag = interp1(f_sparam, ...
                 imag(Sdd21), ...
                 positive_frequency_axis, ...
                 'pchip', ...
                 0);

H_channel_positive = H_real + 1j*H_imag;


% ---------------------------------------------------------
% Build the negative-frequency half using conjugate symmetry
% so time-domain output remains real.
% ---------------------------------------------------------

H_channel_full = ...
    [H_channel_positive; ...
     conj(H_channel_positive(end-1:-1:2))];


% Transform the transmitted waveform to the frequency domain
tx_fft = fft(tx_waveform, Nfft_channel);

% Apply the measured channel transfer function
rx_fft = tx_fft .* H_channel_full;

rx_after_channel_full = real(ifft(rx_fft));

% Keep the original signal-length portion of the output
rx_after_channel = ...
    rx_after_channel_full(1:signal_length);


%% =========================================================
% 4.1 CHANNEL IMPULSE RESPONSE AND DELAY
% =========================================================

channel_impulse_response = real(ifft(H_channel_full));

% Search the first half of the impulse response for the main channel peak.
[~, channel_peak_index] = ...
    max(abs(channel_impulse_response(1:Nfft_channel/2)));

channel_delay_samples = channel_peak_index - 1;

fprintf('Estimated channel delay        = %d samples\n', ...
        channel_delay_samples);

fprintf('Estimated channel delay        = %.2f UI\n', ...
        channel_delay_samples/samples_per_symbol);


figure('Name','Real Channel Impulse Response','Color','w');

plot(channel_impulse_response(1:min(4000,...
     length(channel_impulse_response))), ...
     'LineWidth',1);

title('Impulse Response from Measured S-Parameters');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;


%% =========================================================
% 5) ADD SMALL NOISE
% =========================================================

signal_power = mean(abs(rx_after_channel).^2);

noise_power = ...
    signal_power / ...
    (10^(signal_to_noise_ratio_db/10));

noise = sqrt(noise_power) * ...
        randn(size(rx_after_channel));

rx_noisy = rx_after_channel + noise;


%% =========================================================
% 6) IDEAL RECEIVER GAIN BEFORE ADC
% =========================================================
%
% This gain stage is only used to match the waveform to the ADC input range.
%
% It applies a constant amplitude scaling to the complete waveform.
% No frequency-dependent equalization is introduced here.

peak_rx = max(abs(rx_noisy));

if peak_rx > 0
    rx_adc_input = 0.90 * rx_noisy / peak_rx;
else
    rx_adc_input = rx_noisy;
end


figure('Name','Signal After Real Channel','Color','w');

rx_display_start = channel_delay_samples + 1200;
rx_display_end   = rx_display_start + 500;

plot(rx_display_start:rx_display_end, ...
     rx_adc_input(rx_display_start:rx_display_end), ...
     'LineWidth',1);

title('Signal After Real S-Parameter Channel');
xlabel('Sample Index');
ylabel('Normalized ADC Input');
ylim([-1.1 1.1]);
grid on;

title('Signal After Real S-Parameter Channel');
xlabel('Sample Index');
ylabel('Normalized ADC Input');
ylim([-1.1 1.1]);
grid on;


%% =========================================================
% 7) SAMPLING PHASE SEARCH
% =========================================================
%
% The initial sampling position includes the known RRC group delay
% and the delay introduced by the measured S-parameter channel.

rrc_group_delay = ...
    rrc_filter_span * samples_per_symbol / 2;

base_sampling_offset = ...
    round(rrc_group_delay + ...
          channel_delay_samples + 1);

fprintf('\nInitial sampling offset        = %d samples\n', ...
        base_sampling_offset);


best_phase = 0;
best_phase_score = -inf;

best_sampling_indices = [];
best_adc_digital_codes = [];
best_ffe_input_signal = [];

% Search over a limited symbol range to align the sampled data with the PRBS reference
max_symbol_shift_for_phase = 300;


for phase = 0:samples_per_symbol-1

    first_sample = base_sampling_offset + phase;

    test_sampling_indices = ...
        first_sample + ...
        (0:num_bits-1)*samples_per_symbol;

    test_sampling_indices = ...
        test_sampling_indices( ...
        test_sampling_indices > 0 & ...
        test_sampling_indices <= length(rx_adc_input));

    if length(test_sampling_indices) < 1000
        continue;
    end


    %% -----------------------------------------------------
    % Ideal 4-bit ADC quantization
    %% -----------------------------------------------------

    test_sampled_voltages = ...
        rx_adc_input(test_sampling_indices);

    test_adc_codes = ...
        round(adc_mid_threshold + ...
        test_sampled_voltages * ...
        (adc_mid_threshold - 1));

    test_adc_codes = ...
        min(max(test_adc_codes,0), ...
        adc_total_levels-1);

    % Center the ADC codes around zero before digital equalization
    test_ffe_input = ...
        double(test_adc_codes) - ...
        adc_mid_threshold;


    %% -----------------------------------------------------
    % Find the coarse PRBS alignment for this sampling phase
    %% -----------------------------------------------------

    best_corr_this_phase = -inf;

    for shift = 0:max_symbol_shift_for_phase

        L = min(length(test_ffe_input), ...
                length(symbols_ref)-shift);

        if L > 1000

            x = test_ffe_input(1:L);
            r = symbols_ref(shift+1:shift+L);

            x = x(:) - mean(x(:));
            r = r(:) - mean(r(:));

            denominator = ...
                sqrt(sum(abs(x).^2) * ...
                     sum(abs(r).^2)) + eps;

            corr_val = ...
                abs(sum(x.*r)) / denominator;

            if corr_val > best_corr_this_phase
                best_corr_this_phase = corr_val;
            end
        end
    end


    if best_corr_this_phase > best_phase_score

        best_phase_score = best_corr_this_phase;
        best_phase = phase;

        best_sampling_indices = ...
            test_sampling_indices;

        best_adc_digital_codes = ...
            test_adc_codes;

        best_ffe_input_signal = ...
            test_ffe_input;
    end
end


sampling_indices = best_sampling_indices;
adc_digital_codes = best_adc_digital_codes;
ffe_input_signal = best_ffe_input_signal;

fprintf('Best sampling phase            = %d / %d samples per UI\n', ...
        best_phase, samples_per_symbol);


%% =========================================================
% 8) DISPLAY ADC OUTPUT
% =========================================================

sampled_voltages = rx_adc_input(sampling_indices);

figure('Name','Ideal 4-bit ADC','Color','w');

subplot(2,1,1);

segment_start = sampling_indices(200);

segment_end = min(segment_start + ...
                  20*samples_per_symbol, ...
                  length(rx_adc_input));

plot(segment_start:segment_end, ...
     rx_adc_input(segment_start:segment_end), ...
     'LineWidth',1);

hold on;

valid_samples = ...
    sampling_indices( ...
    sampling_indices >= segment_start & ...
    sampling_indices <= segment_end);

plot(valid_samples, ...
     rx_adc_input(valid_samples), ...
     'o');

title('Selected Sampling Points');
xlabel('Sample Index');
ylabel('ADC Input');
grid on;


subplot(2,1,2);

stairs(adc_digital_codes(1:min(120,...
       length(adc_digital_codes))), ...
       'LineWidth',1.2);

title('Ideal 4-bit ADC Output');
xlabel('Symbol Index');
ylabel('ADC Code');

ylim([-1 adc_total_levels]);

grid on;


%% =========================================================
% 9) ALIGN ADC SAMPLES WITH PRBS REFERENCE
% =========================================================

max_symbol_shift = 300;

best_shift = 0;
best_corr = -inf;


for shift = 0:max_symbol_shift

    L = min(length(ffe_input_signal), ...
            length(symbols_ref)-shift);

    if L > 100

        test_input = ...
            ffe_input_signal(1:L);

        test_ref = ...
            symbols_ref(shift+1:shift+L);

        test_input = ...
            test_input(:) - mean(test_input(:));

        test_ref = ...
            test_ref(:) - mean(test_ref(:));

        denominator = ...
            sqrt(sum(abs(test_input).^2) * ...
                 sum(abs(test_ref).^2)) + eps;

        corr_val = ...
            abs(sum(test_input.*test_ref)) / ...
            denominator;

        if corr_val > best_corr

            best_corr = corr_val;
            best_shift = shift;

        end
    end
end


num_avail = ...
    min(length(ffe_input_signal), ...
        length(symbols_ref)-best_shift);

ffe_input_signal = ...
    ffe_input_signal(1:num_avail);

ref_bits = ...
    symbols_ref(best_shift+1: ...
                best_shift+num_avail);

fprintf('Best reference shift           = %d symbols\n', ...
        best_shift);


%% =========================================================
% 10) BUILD 41-TAP RX-FFE MATRIX
% =========================================================

ffe_center_tap = ...
    ceil(num_ffe_taps/2);

padded_in = ...
    [zeros(ffe_center_tap-1,1); ...
     ffe_input_signal(:); ...
     zeros(num_ffe_taps-ffe_center_tap,1)];

U_matrix = ...
    zeros(num_avail,num_ffe_taps);


for i = 1:num_avail

    U_matrix(i,:) = ...
        padded_in(i:i+num_ffe_taps-1).';

end


num_train = ...
    floor(training_data_fraction*num_avail);


%% =========================================================
% 11) REGULARIZED LEAST-SQUARES FFE TRAINING
%     WITH FFE BOOST MONITORING
% =========================================================

% Candidate regularization strengths.
% Each value is scaled relative to the energy in U'*U so that
% the regularization remains consistent with the current input data.
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

U_train = U_matrix(1:num_train,:);
d_train = ref_bits(1:num_train);

val_start = num_train + 1;

if val_start > num_avail - 1000
    val_start = floor(0.6*num_avail);
end

val_idx = val_start:num_avail;

% ---------------------------------------------------------
% Scale lambda relative to the input correlation matrix
% ---------------------------------------------------------
R = U_train.' * U_train;

regularization_scale = ...
    trace(R) / num_ffe_taps;

% ---------------------------------------------------------
% Target FFE boost range used during candidate selection
% ---------------------------------------------------------
min_desired_boost_db = 6;
max_desired_boost_db = 12;

ffe_sample_rate = data_rate;

best_score = -inf;
best_weights = [];
best_output = [];

best_alpha = NaN;
best_lambda = NaN;

best_eye_opening_pct = -inf;
best_boost_db = NaN;
best_ber = inf;

fprintf('\n========================================================\n');
fprintf('       RX-FFE REGULARIZATION / BOOST SEARCH\n');
fprintf('========================================================\n');

for ai = 1:length(alpha_list)

    alpha = alpha_list(ai);

    lambda_reg = ...
        alpha * regularization_scale;

    % -----------------------------------------------------
    % Solve the regularized least-squares problem
    % -----------------------------------------------------
    w_try = ...
        (R + lambda_reg*eye(num_ffe_taps)) \ ...
        (U_train.' * d_train);

    y_try = U_matrix * w_try;

    % -----------------------------------------------------
    % Correct a possible polarity inversion
    % -----------------------------------------------------
    if mean(y_try(ref_bits == 1)) < ...
       mean(y_try(ref_bits == -1))

        y_try = -y_try;
        w_try = -w_try;
    end

    % -----------------------------------------------------
    % Normalize the equalized output level
    % -----------------------------------------------------
    y_val = y_try(val_idx);
    ref_val = ref_bits(val_idx);

    pos_val = y_val(ref_val > 0);
    neg_val = y_val(ref_val < 0);

    level_sep = ...
        mean(pos_val) - mean(neg_val);

    if level_sep > 0
        gain = 2 / level_sep;
    else
        gain = 1;
    end

    y_try = y_try * gain;
    w_try = w_try * gain;

    % -----------------------------------------------------
    % Evaluate BER and digital eye opening
    % -----------------------------------------------------
    y_val = y_try(val_idx);

    pos_val = y_val(ref_val > 0);
    neg_val = y_val(ref_val < 0);

    p5_pos = ...
        local_percentile(pos_val,5);

    p95_neg = ...
        local_percentile(neg_val,95);

    eye_opening = ...
        p5_pos - p95_neg;

    eye_opening_pct = ...
        100 * eye_opening / 2;

    decisions_val = ...
        y_val >= 0;

    bits_val = ...
        ref_val > 0;

    ber_val = ...
        sum(decisions_val(:) ~= bits_val(:)) / ...
        length(bits_val);

    % -----------------------------------------------------
    % Measure the FFE boost at the 16-GHz Nyquist frequency
    % -----------------------------------------------------
    [H_try, f_try] = ...
        freqz(w_try,1,8192,ffe_sample_rate);

    H_try_db = ...
        20*log10(abs(H_try)+eps);

    % Relative FFE response:
    % high-frequency gain referenced to the DC gain
    H_try_relative_db = ...
        H_try_db - H_try_db(1);

    [~, idx_try_16] = ...
        min(abs(f_try-data_nyquist_hz));

    boost_16_db = ...
        H_try_relative_db(idx_try_16);

    fprintf(['alpha = %-8.1e | ' ...
             'lambda = %-10.3e | ' ...
             'BER = %.3e | ' ...
             'eye = %6.2f%% | ' ...
             'FFE boost = %6.2f dB\n'], ...
             alpha, ...
             lambda_reg, ...
             ber_val, ...
             eye_opening_pct, ...
             boost_16_db);

    % -----------------------------------------------------
    %
    % FFE boost within the selected 6-12 dB range.
    %
    % Within this range:
    %   1) Lower BER has the highest priority.
    %   2) A larger digital eye opening is preferred.
    %
    % The tap coefficients are still obtained from the regularized
    % least-squares solution rather than being forced to a fixed boost.
    % -----------------------------------------------------

    boost_is_acceptable = ...
        boost_16_db >= min_desired_boost_db && ...
        boost_16_db <= max_desired_boost_db;

    if boost_is_acceptable

        score = ...
            eye_opening_pct - ...
            1000*ber_val;

        if score > best_score

            best_score = score;

            best_weights = w_try;
            best_output = y_try;

            best_alpha = alpha;
            best_lambda = lambda_reg;

            best_eye_opening_pct = ...
                eye_opening_pct;

            best_boost_db = ...
                boost_16_db;

            best_ber = ...
                ber_val;
        end
    end
end


%% ---------------------------------------------------------
% If no candidate falls within the 6-12 dB range,
% keep the best unconstrained solution instead of
% forcing the equalizer to meet the target range.
%% ---------------------------------------------------------

if isempty(best_weights)

    warning(['No tested regularization value produced an FFE ' ...
             'boost between 6 and 12 dB.']);

    fprintf(['Keeping the minimum-error solution instead. ' ...
             'Do NOT force the result artificially.\n']);

    best_total_score = -inf;

    for ai = 1:length(alpha_list)

        alpha = alpha_list(ai);
        lambda_reg = alpha * regularization_scale;

        w_try = ...
            (R + lambda_reg*eye(num_ffe_taps)) \ ...
            (U_train.'*d_train);

        y_try = U_matrix*w_try;

        if mean(y_try(ref_bits==1)) < ...
           mean(y_try(ref_bits==-1))

            y_try = -y_try;
            w_try = -w_try;
        end

        y_val = y_try(val_idx);
        ref_val = ref_bits(val_idx);

        pos_val = y_val(ref_val>0);
        neg_val = y_val(ref_val<0);

        level_sep = ...
            mean(pos_val)-mean(neg_val);

        if level_sep > 0
            gain = 2/level_sep;
        else
            gain = 1;
        end

        y_try = y_try*gain;
        w_try = w_try*gain;

        y_val = y_try(val_idx);

        pos_val = y_val(ref_val>0);
        neg_val = y_val(ref_val<0);

        eye_opening_pct = ...
            100 * ...
            (local_percentile(pos_val,5) - ...
             local_percentile(neg_val,95)) / 2;

        ber_val = ...
            sum((y_val>=0) ~= (ref_val>0)) / ...
            length(ref_val);

        score = ...
            eye_opening_pct - ...
            1000*ber_val;

        if score > best_total_score

            best_total_score = score;

            best_weights = w_try;
            best_output = y_try;

            best_alpha = alpha;
            best_lambda = lambda_reg;

            best_eye_opening_pct = ...
                eye_opening_pct;

            best_ber = ber_val;

            [H_try,f_try] = ...
                freqz(w_try,1,8192,ffe_sample_rate);

            H_try_db = ...
                20*log10(abs(H_try)+eps);

            H_try_relative_db = ...
                H_try_db-H_try_db(1);

            [~,idx_try_16] = ...
                min(abs(f_try-data_nyquist_hz));

            best_boost_db = ...
                H_try_relative_db(idx_try_16);
        end
    end
end


weights = best_weights;
equalized_output = best_output;

fprintf('\n========================================================\n');
fprintf('                SELECTED RX-FFE\n');
fprintf('========================================================\n');

fprintf('Selected alpha                  = %.3e\n', ...
        best_alpha);

fprintf('Selected lambda                 = %.3e\n', ...
        best_lambda);

fprintf('Validation BER                  = %.3e\n', ...
        best_ber);

fprintf('Validation eye opening          = %.2f %%\n', ...
        best_eye_opening_pct);

fprintf('Selected FFE boost at 16 GHz    = %.2f dB\n', ...
        best_boost_db);

fprintf('========================================================\n');

%% =========================================================
% 12) ALIGN FFE OUTPUT
% =========================================================

max_decision_delay = 50;

best_decision_delay = 0;
best_errors = inf;


for d = 0:max_decision_delay

    L = min(length(equalized_output)-d, ...
            length(ref_bits));

    if L > 100

        test_output = ...
            equalized_output(d+1:d+L);

        test_ref = ...
            ref_bits(1:L);

        test_decisions = ...
            test_output >= 0;

        test_bits = ...
            test_ref > 0;

        test_errors = ...
            sum(test_decisions(:) ~= ...
                test_bits(:));


        if test_errors < best_errors

            best_errors = ...
                test_errors;

            best_decision_delay = d;

        end
    end
end


fprintf('FFE decision delay             = %d symbols\n', ...
        best_decision_delay);


L_final = ...
    min(length(equalized_output)- ...
        best_decision_delay, ...
        length(ref_bits));


equalized_output = ...
    equalized_output( ...
    best_decision_delay+1: ...
    best_decision_delay+L_final);


ref_bits = ref_bits(1:L_final);

ffe_input_signal = ...
    ffe_input_signal(1:L_final);


if mean(equalized_output(ref_bits==1)) < ...
   mean(equalized_output(ref_bits==-1))

    equalized_output = ...
        -equalized_output;

    weights = -weights;

end


%% =========================================================
% 13) BER BEFORE AND AFTER RX-FFE
% =========================================================

reference_bits = ref_bits > 0;


% Before RX-FFE
decided_before = ...
    ffe_input_signal >= 0;

errors_before = ...
    sum(decided_before(:) ~= ...
        reference_bits(:));

ber_before_ffe = ...
    errors_before / ...
    length(reference_bits);


% After RX-FFE
decided_after = ...
    equalized_output >= 0;

errors_after = ...
    sum(decided_after(:) ~= ...
        reference_bits(:));

ber_after_ffe = ...
    errors_after / ...
    length(reference_bits);


fprintf('\n========================================================\n');
fprintf('                    BER RESULTS\n');
fprintf('========================================================\n');

fprintf('BER before RX-FFE = %.3e (%d errors)\n', ...
        ber_before_ffe, ...
        errors_before);

fprintf('BER after RX-FFE  = %.3e (%d errors)\n', ...
        ber_after_ffe, ...
        errors_after);


%% =========================================================
% 14) DIGITAL EYE OPENING AFTER FFE
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


fprintf('\nDigital eye opening after FFE = %.4f\n', ...
        actual_eye_opening);

fprintf('Digital eye opening percent   = %.2f %%\n', ...
        actual_eye_opening_pct);


%% =========================================================
% 15) RX-FFE FREQUENCY RESPONSE
% =========================================================
%
% This section evaluates the frequency-domain behavior of the selected RX-FFE.
%
% The FFE boost is measured at the 16-GHz Nyquist frequency
% relative to the low-frequency response of the equalizer.

ffe_sample_rate = data_rate;

[Hffe, fffe] = ...
    freqz(weights, ...
          1, ...
          8192, ...
          ffe_sample_rate);


Hffe_db = ...
    20*log10(abs(Hffe)+eps);


% Normalize the FFE response to its low-frequency (DC) gain
Hffe_relative_db = ...
    Hffe_db-Hffe_db(1);


[~, idx_ffe_16] = ...
    min(abs(fffe-data_nyquist_hz));


ffe_boost_16_db = ...
    Hffe_relative_db(idx_ffe_16);


fprintf('\n========================================================\n');
fprintf('               FFE EQUALIZATION RESULT\n');
fprintf('========================================================\n');

fprintf('RX-FFE boost near 16 GHz       = %.2f dB\n', ...
        ffe_boost_16_db);


figure('Name','RX-FFE Frequency Response','Color','w');

plot(fffe/1e9, ...
     Hffe_relative_db, ...
     'LineWidth',1.5);

hold on;

plot(fffe(idx_ffe_16)/1e9, ...
     ffe_boost_16_db, ...
     'o', ...
     'MarkerSize',8, ...
     'LineWidth',1.5);

xline(16,'--','16 GHz');

title('Digital RX-FFE Equalization Response');
xlabel('Frequency [GHz]');
ylabel('Relative FFE Gain [dB]');

xlim([0 16]);

grid on;


%% =========================================================
% 16) CHANNEL + FFE FREQUENCY-DOMAIN VIEW
%     USING ABSOLUTE CHANNEL LOSS AT 16 GHz
% =========================================================

% Interpolate the measured channel response onto the FFE frequency axis
channel_abs_db = ...
    interp1(f_sparam, ...
            Sdd21_db, ...
            fffe, ...
            'linear');

% Plot the FFE response relative to its low-frequency gain
ffe_relative_db = Hffe_relative_db;

% Combined response:
% measured channel loss plus relative FFE compensation
combined_abs_db = ...
    channel_abs_db + ...
    ffe_relative_db;

figure('Name','Channel and RX-FFE Equalization','Color','w');

plot(fffe/1e9, ...
     channel_abs_db, ...
     'LineWidth',1.5);

hold on;

plot(fffe/1e9, ...
     ffe_relative_db, ...
     'LineWidth',1.5);

plot(fffe/1e9, ...
     combined_abs_db, ...
     'LineWidth',1.5);

%% ---------------------------------------------------------
% MARK RESPONSE VALUES AT 16 GHz
%% ---------------------------------------------------------

[~, idx16_plot] = ...
    min(abs(fffe - data_nyquist_hz));

channel_16_plot = ...
    channel_abs_db(idx16_plot);

ffe_16_plot = ...
    ffe_relative_db(idx16_plot);

combined_16_plot = ...
    combined_abs_db(idx16_plot);

% Mark the 16-GHz values
plot(16, ...
     channel_16_plot, ...
     'o', ...
     'MarkerSize',8, ...
     'LineWidth',1.5);

plot(16, ...
     ffe_16_plot, ...
     'o', ...
     'MarkerSize',8, ...
     'LineWidth',1.5);

plot(16, ...
     combined_16_plot, ...
     'o', ...
     'MarkerSize',8, ...
     'LineWidth',1.5);

% Add numerical labels for the 16-GHz values
text(13.9, ...
     channel_16_plot-0.8, ...
     sprintf('Channel = %.2f dB', channel_16_plot), ...
     'FontSize',10);

text(12.8, ...
     ffe_16_plot+0.8, ...
     sprintf('RX-FFE = +%.2f dB', ffe_16_plot), ...
     'FontSize',10);

text(13.5, ...
     combined_16_plot+0.8, ...
     sprintf('Residual = %.2f dB', combined_16_plot), ...
     'FontSize',10);

xline(16,'--','16 GHz');
yline(0,'--');

title('Channel Loss and Digital RX-FFE Compensation');

xlabel('Frequency [GHz]');
ylabel('Magnitude [dB]');

legend('Measured Residual Channel', ...
       'Digital RX-FFE Compensation', ...
       'Residual After RX-FFE', ...
       'Location','best');

xlim([0 16]);

grid on;


%% =========================================================
% 17) EYE BEFORE RX-FFE
% =========================================================

samples_per_eye = ...
    2*samples_per_symbol;

time_axis_before = ...
    linspace(-1,1,...
    samples_per_eye+1);


% Center the eye traces around the selected sampling location
first_eye_start = ...
    sampling_indices(100) - ...
    samples_per_symbol;


num_traces_before = ...
    min(max_eye_traces, ...
    floor((length(rx_adc_input) - ...
    first_eye_start - ...
    samples_per_eye) / ...
    samples_per_symbol));


matrix_before = ...
    zeros(samples_per_eye+1, ...
          num_traces_before);


idx = first_eye_start;


for t = 1:num_traces_before

    if idx > 0 && ...
       idx+samples_per_eye <= ...
       length(rx_adc_input)

        matrix_before(:,t) = ...
            rx_adc_input( ...
            idx:idx+samples_per_eye);

    end

    idx = idx + ...
          samples_per_symbol;

end


%% =========================================================
% 18) BEFORE / AFTER FFE VISUALIZATION
% =========================================================

figure('Name','Eye Before and After RX-FFE', ...
       'Color','w', ...
       'Position',[100 100 850 650]);


%% ---------------------------------------------------------
% Before RX-FFE
%% ---------------------------------------------------------

subplot(2,1,1);

plot(time_axis_before, ...
     matrix_before, ...
     'LineWidth',0.5);

title('Eye BEFORE Digital RX-FFE');

xlabel('Time [UI]');
ylabel('Normalized Amplitude');

ylim([-1.2 1.2]);

grid on;


%% ---------------------------------------------------------
% After RX-FFE
%% ---------------------------------------------------------
%
% The FFE output contains one sample per symbol.
% The plot is therefore a digital symbol-eye representation
% rather than a reconstructed oversampled analog waveform.

subplot(2,1,2);

hold on;

num_after_eye_traces = ...
    min(600, ...
        length(equalized_output)-2);


after_start = ...
    max(1, ...
        floor(0.8*length(equalized_output)) - ...
        num_after_eye_traces);


after_end = ...
    min(after_start + ...
        num_after_eye_traces, ...
        length(equalized_output)-2);


for k = after_start:after_end

    plot([-1 0 1], ...
         [equalized_output(k), ...
          equalized_output(k+1), ...
          equalized_output(k+2)], ...
         'LineWidth',0.5);

end


yline(0,'--k');

title('Digital Symbol Eye AFTER RX-FFE');

xlabel('Symbol Time [UI]');
ylabel('Equalized Symbol Level');

ylim([-1.5 1.5]);

grid on;


%% =========================================================
% 19) FFE TAP COEFFICIENTS
% =========================================================

figure('Name','RX-FFE Tap Coefficients','Color','w');

stem(1:num_ffe_taps, ...
     weights, ...
     'filled');

title('41-Tap Digital RX-FFE Coefficients');

xlabel('Tap Number');
ylabel('Coefficient');

grid on;


%% =========================================================
% 20) FINAL SUMMARY
% =========================================================

fprintf('\n========================================================\n');
fprintf('       FINAL MENTOR-REQUESTED FFE VALIDATION\n');
fprintf('========================================================\n');

fprintf('Data rate                       = %.1f GT/s\n', ...
        data_rate/1e9);

fprintf('Nyquist frequency               = %.1f GHz\n', ...
        data_nyquist_hz/1e9);

fprintf('Real channel loss at 16 GHz     = %.2f dB\n', ...
        channel_loss_16_db);

fprintf('ADC resolution                  = %d bits\n', ...
        adc_resolution_bits);

fprintf('RX-FFE taps                     = %d\n', ...
        num_ffe_taps);

fprintf('RX-FFE boost near 16 GHz        = %.2f dB\n', ...
        ffe_boost_16_db);

fprintf('BER before RX-FFE               = %.3e\n', ...
        ber_before_ffe);

fprintf('BER after RX-FFE                = %.3e\n', ...
        ber_after_ffe);

fprintf('Digital eye opening after FFE   = %.2f %%\n', ...
        actual_eye_opening_pct);

fprintf('Selected alpha                  = %.3e\n', ...
        best_alpha);

fprintf('Selected lambda                 = %.3e\n', ...
        best_lambda);

fprintf('========================================================\n');


%% =========================================================
% LOCAL FUNCTION: PRBS31
% =========================================================

function bits = generate_prbs31(N)

    shift_register = ones(31,1);
    bits = zeros(N,1);

    for i = 1:N

        bits(i) = ...
            shift_register(end);

        feedback = ...
            xor(shift_register(end), ...
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