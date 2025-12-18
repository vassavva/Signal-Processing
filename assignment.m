%Task 1
clc;
close all;
clear all;

[x,fs] = audioread('Vasiliki Savva.wav');   %extract signal from audio
x = x(:);
Ts = 1/fs;                                  %sampling period
N = length(x); 
t = 0:Ts:(N-1)*Ts;                          %generate discrete time values (nTs)

figure; plot(t,x);                          %plot time-domain
xlabel('time(sec)'); ylabel('Amplitude'); title('Time Domain');


x_fft = fft(x);                             %calculate FFT                             
k = (0:1:N-1)';
f = k*(fs/N);                               %normalise frequency

x_norm= abs(x_fft)/N;
x_db=20*log10(x_norm);
figure; plot(f,x_db);              
xlim([0 fs/2]);
xlabel('Discrete Frequency (Hz)'); ylabel('Magnitude');
title('Spectrum');


samples_to_show = min(1000, N);


n_window = (0:N-1)';

han_w = 0.5*(1-cos(2*pi*k/(N-1)));              %generate N point hanning window
rec_w = ones(1, N)';                                 % Changed to row vector
ham_w = (0.54 - 0.46 * cos(2 * pi * (0:N-1) / (N - 1)))';



xh = x.*han_w;                                  %apply window to signal
x_fft_h = fft(xh);
xh_norm= abs(x_fft_h)/N;
xh_db=20*log10(xh_norm);

xr = x.*rec_w;                                  %apply window to signal
x_fft_r = fft(xr);
xr_norm= abs(x_fft_r)/N;
xr_db=20*log10(xr_norm);

xm = x.*ham_w;                                  %apply window to signal
x_fft_m = fft(xm);
xm_norm= abs(x_fft_m)/N;
xm_db=20*log10(xm_norm);

figure;


plot(f,xr_db,'b'); 
hold on;
plot(f,xm_db, 'g');
plot(f,xh_db,'r'); 
hold off;
xlabel('Discrete Frequency (Hz)'); ylabel('Magnitude');
title('')
legend('Hamming', 'Rectagular', 'Hanning');
grid on
xlim([0 fs/2])



fc = 18000;              % FROM THE GRAPH ESTIMATE 18K 
BW = 8000;              % use an 8kHz bandwidth
fmin = fc - BW/2;        
fmax = fc + BW/2    ;



%Task 2

dF = 2000;          % transition width

% Filter lengths
Nham = ceil(3.3*fs/dF);  if mod(Nham,2)==0, Nham=Nham+1; end
Nblc = ceil(5.5*fs/dF);  if mod(Nblc,2)==0, Nblc=Nblc+1; end

% Normalized frequencies
Fc1 = fmin/fs;
Fc2 = fmax/fs;

% Ideal impulse responses

% Hamming
n_h = 0:Nham-1;
a_h = (Nham-1)/2;
h_ideal_ham = 2*Fc2*sinc(2*Fc2*(n_h - a_h)) - 2*Fc1*sinc(2*Fc1*(n_h - a_h));

% Blackman
n_b = 0:Nblc-1;
a_b = (Nblc-1)/2;
h_ideal_blk = 2*Fc2*sinc(2*Fc2*(n_b - a_b)) - 2*Fc1*sinc(2*Fc1*(n_b - a_b));

% Windows (manual)

w_ham = 0.54 - 0.46*cos(2*pi*n_h/(Nham-1));
w_blk = 0.42 - 0.5*cos(2*pi*n_b/(Nblc-1)) + 0.08*cos(4*pi*n_b/(Nblc-1));

% Apply windows
h_ham = h_ideal_ham .* w_ham;
h_blk = h_ideal_blk .* w_blk;

% Verify frequency responses

figure;
subplot(2,1,1); freqz(h_ham,1,4096,fs); title('Hamming Window FIR Bandpass');
subplot(2,1,2); freqz(h_blk,1,4096,fs); title('Blackman Window FIR Bandpass');

% Manual FIR convolution

Lx = length(x);

% Preallocate output
y_ham = zeros(1, Lx);
y_blk = zeros(1, Lx);

% Hamming filter
for n = 1:Lx
    acc = 0;
    for k = 1:Nham
        if (n - k + 1) > 0
            acc = acc + h_ham(k)*x(n - k + 1);
        end
    end
    y_ham(n) = acc;
end

% Blackman filter
for n = 1:Lx
    acc = 0;
    for k = 1:Nblc
        if (n - k + 1) > 0
            acc = acc + h_blk(k)*x(n - k + 1);
        end
    end
    y_blk(n) = acc;
end
% Time-domain comparison

t = (0:Lx-1)/fs;

figure;
subplot(3,1,1); plot(t,x); title('Original AM Signal'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,2); plot(t,y_ham); title('Filtered AM Signal (Hamming)'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(3,1,3); plot(t,y_blk); title('Filtered AM Signal (Blackman)'); xlabel('Time (s)'); ylabel('Amplitude');

% Frequency-domain comparison

Nfft = 8192;
f_fft = linspace(-fs/2, fs/2, Nfft);

X  = fftshift(abs(fft(x, Nfft)));
Yh = fftshift(abs(fft(y_ham, Nfft)));
Yb = fftshift(abs(fft(y_blk, Nfft)));

figure;
subplot(3,1,1); plot(f_fft,X); title('Original AM Signal Spectrum'); xlabel('Hz');
subplot(3,1,2); plot(f_fft,Yh); title('Hamming Fsssiltered Spectrum'); xlabel('Hz');
subplot(3,1,3); plot(f_fft,Yb); title('Blackman Filtered Spectrum'); xlabel('Hz');


% Task 3 – Carrier Recovery & Mixing (Corrected)

x_bp = y_blk;           % Blackman FIR output

% --- Square law for carrier recovery
x_sq = x_bp.^2;

% Apply Hamming window to reduce spectral leakage
h_win = hamming(length(x_sq))';
x_sq_windowed = x_sq .* h_win;

% FFT with high resolution
Nfft = 2^16;
Xsq = fft(x_sq_windowed, Nfft);

% Frequency vector (0 to fs)
f_sq = (0:Nfft-1)*(fs/Nfft);

% Expected carrier from Task 1
expected_fc = 18000;   % Hz

% Search range around 2*fc
search_range_low  = 2*expected_fc - 5000;  % Hz
search_range_high = 2*expected_fc + 5000;  % Hz
idx_band = (f_sq >= search_range_low) & (f_sq <= search_range_high);

% Find the peak in squared signal spectrum
[~, idx_max] = max(abs(Xsq(idx_band)));

% Map peak back to frequency
f_2fc = f_sq(idx_band);
f_2fc = f_2fc(idx_max);

% Compute estimated carrier frequency
fc_est = f_2fc / 2;
fprintf('Estimated carrier frequency fc = %.0f Hz\n', fc_est);

% --- Generate local carrier
t_bp = (0:length(x_bp)-1)/fs;
phi = 0;  % initial phase
local_carrier = cos(2*pi*fc_est*t_bp + phi);

% --- Mixing (multiply bandpass signal with local carrier)
x_mix = x_bp .* local_carrier;

% --- Time domain plots
figure;
subplot(2,1,1); plot(t_bp, x_bp); title('Bandpass Filtered AM Signal'); xlabel('Time (s)'); ylabel('Amplitude');
subplot(2,1,2); plot(t_bp, x_mix); title('After Mixing with Local Carrier'); xlabel('Time (s)'); ylabel('Amplitude');

% --- Frequency domain plot
Xmix = fftshift(abs(fft(x_mix, Nfft)));
f_fft = linspace(0, fs, Nfft);
figure;
plot(f_sq, 20*log10(abs(Xsq)+eps)); hold on;
xline(2*expected_fc, 'r--', 'LineWidth', 2);  % Expected 2fc
xline(f_2fc, 'g--', 'LineWidth', 2);         % Detected 2fc
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title('Squared Signal Spectrum for Carrier Detection');
grid on;
legend('Spectrum','Expected 2fc','Detected 2fc');

% Task 4 – IIR Lowpass Filter (Butterworth)



% Assume x_mix is the mixed signal from Task 3 (Blackman FIR output multiplied by local carrier)
signal = x_mix;   % input to IIR lowpass

% Design Butterworth filter
order = 4;
fc_lp = 4000;       % cutoff frequency in Hz
Wn = fc_lp/(fs/2);  % normalized cutoff (0–1, Nyquist)


[b, a] = butter(order, Wn, 'low');
% Verify frequency response

figure;
freqz(b, a, 4096, fs);
title('4th-order Butterworth Lowpass Filter');

% Manual filtering via difference equation

Lx = length(signal);
y_iir = zeros(1,Lx);

% Apply difference equation: y[n] = b0*x[n] + b1*x[n-1] + ... - a1*y[n-1] - ...
for n = 1:Lx
    acc_b = 0;  % numerator sum
    for k = 1:length(b)
        if (n - k + 1) > 0
            acc_b = acc_b + b(k)*signal(n - k + 1);
        end
    end
    
    acc_a = 0;  % denominator sum
    for k = 2:length(a)
        if (n - k + 1) > 0
            acc_a = acc_a + a(k)*y_iir(n - k + 1);
        end
    end
    
    y_iir(n) = acc_b - acc_a;  % difference equation
end

% Plot time-domain

t_signal = (0:Lx-1)/fs;
figure;
subplot(2,1,1);
plot(t_signal, signal);
title('Input to IIR Lowpass Filter (Mixed Signal)');
xlabel('Time (s)'); ylabel('Amplitude');

subplot(2,1,2);
plot(t_signal, y_iir);
title('Output of IIR Lowpass Filter (Manual)');
xlabel('Time (s)'); ylabel('Amplitude');

% Frequency-domain plot

Y_iir = fftshift(abs(fft(y_iir, Nfft)));

figure;
plot(f_sq, Y_iir);
title('Spectrum After IIR Lowpass Filtering');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Task 5 – Improved Phase Adjustment and Playback


signal_lp = y_iir;  % Output of manual IIR lowpass (Task 4)
y_final = zeros(size(signal_lp));

% --- Stage 1: Coarse Phase Search ---
phi_coarse = linspace(0, pi, 50);  % coarse 50-point grid
rms_coarse = zeros(size(phi_coarse));

for p = 1:length(phi_coarse)
    phi = phi_coarse(p);
    local_carrier = cos(2*pi*fc_est*t_bp + phi);
    mixed_signal = y_blk .* local_carrier;  % y_blk = filtered bandpass from Task 2
    
    % Apply same IIR filter
   % --- Manual IIR filtering (inline) ---
Lx = length(mixed_signal);
y_temp = zeros(1,Lx);

for n = 1:Lx
    acc_b = 0;
    for k = 1:length(b)
        if (n - k + 1) > 0
            acc_b = acc_b + b(k)*mixed_signal(n - k + 1);
        end
    end
    
    acc_a = 0;
    for k = 2:length(a)
        if (n - k + 1) > 0
            acc_a = acc_a + a(k)*y_temp(n - k + 1);
        end
    end
    
    y_temp(n) = acc_b - acc_a;
end

    
    rms_coarse(p) = sqrt(mean(y_temp.^2));  % RMS metric
end

[~, idx_best_coarse] = max(rms_coarse);
phi_best_coarse = phi_coarse(idx_best_coarse);

% --- Stage 2: Fine Phase Search Around Coarse Optimum ---
phi_fine = linspace(max(0, phi_best_coarse-pi/20), min(pi, phi_best_coarse+pi/20), 41);  % +/- 9° around coarse
rms_fine = zeros(size(phi_fine));

for p = 1:length(phi_fine)
    phi = phi_fine(p);
    local_carrier = cos(2*pi*fc_est*t_bp + phi);
    mixed_signal = y_blk .* local_carrier;
    % --- Manual IIR filtering (inline) ---
Lx = length(mixed_signal);
y_temp = zeros(1,Lx);

for n = 1:Lx
    acc_b = 0;
    for k = 1:length(b)
        if (n - k + 1) > 0
            acc_b = acc_b + b(k)*mixed_signal(n - k + 1);
        end
    end
    
    acc_a = 0;
    for k = 2:length(a)
        if (n - k + 1) > 0
            acc_a = acc_a + a(k)*y_temp(n - k + 1);
        end
    end
    
    y_temp(n) = acc_b - acc_a;
end

    rms_fine(p) = sqrt(mean(y_temp.^2));
end

[~, idx_best_fine] = max(rms_fine);
phi_best = phi_fine(idx_best_fine);

% --- Apply best phase and filter ---
local_carrier_best = cos(2*pi*fc_est*t_bp + phi_best);
mixed_signal_best = y_blk .* local_carrier_best;
% --- Apply best phase and filter ---
Lx = length(mixed_signal_best);
y_final = zeros(1,Lx);

for n = 1:Lx
    acc_b = 0;
    for k = 1:length(b)
        if (n - k + 1) > 0
            acc_b = acc_b + b(k)*mixed_signal_best(n - k + 1);
        end
    end
    
    acc_a = 0;
    for k = 2:length(a)
        if (n - k + 1) > 0
            acc_a = acc_a + a(k)*y_final(n - k + 1);
        end
    end
    
    y_final(n) = acc_b - acc_a;
end


fprintf('Optimal phase selected: %.4f rad (%.2f degrees)\n', phi_best, phi_best*180/pi);

% --- Optional: Simple SNR check ---
signal_band = 300:3400;  % speech band
fft_y = abs(fft(y_final));
Nfft = length(fft_y);
freq_axis = (0:Nfft-1)/Nfft*fs;
speech_idx = freq_axis >= signal_band(1) & freq_axis <= signal_band(2);
noise_idx = freq_axis > 4000 & freq_axis <= fs/2;  % above speech band
SNR_est = 10*log10(sum(fft_y(speech_idx).^2)/sum(fft_y(noise_idx).^2));
fprintf('Estimated SNR: %.2f dB\n', SNR_est);

% --- Normalize and Playback ---
x_playback = y_final / max(abs(y_final)) * 0.9;
sound(x_playback, fs);
pause(length(x_playback)/fs + 0.5);
fprintf('Audio playback complete.\n');

% --- Save audio file ---
audiowrite('demodulated_message_improved.wav', x_playback, fs);
fprintf('Audio saved as demodulated_message_improved.wav\n');
