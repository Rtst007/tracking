% arbitrary waveform sig proc 

function Phase_Coded_Pulse_Compression

close all
clear
clc

%% User Inputs

% target range (m), relative to range window start
tgt.rng = [40, 43, 150];

% target radial velocity (m/s); negative for in-bound
tgt.vel = [0, 000, 000];

% target RCS (dBsm);
tgt.rcs = [0, 0, -5];

% minimum detectable signal (dBsm)
tgt.mds = -25;

% maximum sampling rate (Hz)
sp.max_sampling_rate = 3e9;

% chip duration (sec)
sp.tChip = 1e-8;

% number of registers for generating pseudo-random noise code
sp.nRegisters = 10;

% range window length (m)
sp.rw_length = 200;

% carrier frequency (Hz)
sp.fo = 9e9;

% radial velocity (m/s) for replica doppler compensation 
replica_vel = [0, 000, 000];


%% Signal Processing

% initialize the required signal processing variables
sp = initialize_prn_sp(sp);

% generate the replica(s)
H_f = arrayfun(@(x) get_prn_replica_spectrum(sp, x), replica_vel, 'uniformOutput', 0);

% generate the (noisy) receive signal spectrum
X_f = get_prn_return_signal_spectrum(tgt, sp);
% 
% perform pulse compression (all replicas)
y_t = cellfun(@(x) abs(ifft(X_f .* x)), H_f, 'uniformOutput', 0);


%% Plotting

figure;%('units', 'normalized', 'position', [0.1, 0.10, 0.3, 0.8]);

% rcs vs. range, no doppler compensation
% subplot(3,1,1)
plot(sp.rng, 20*log10(y_t{1}(1 : sp.rw_samps)), 'b'); grid on;
%ylim([-60, 0]); 
ylabel('RCS (dBsm)');
title('Pulse compression: uncompensated replica');

% rcs vs. range, with 3 Km/s doppler compensation
% subplot(3,1,2)
% plot(sp.rng, 20*log10(y_t{2}(1 : sp.rw_samps)), 'b'); grid on;
% ylim([-60, 0]); ylabel('RCS (dBsm)');
% title('Pulse compression: -3 km/s range rate compensated replica');
% 
% % rcs vs. range, with 6 Km/s doppler compensation
% subplot(3,1,3)
% plot(sp.rng, 20*log10(y_t{3}(1 : sp.rw_samps)), 'b'); grid on;
% ylim([-60, 0]); ylabel('RCS (dBsm)'); xlabel('Range Window (m)');
% title('Pulse compression: -6 km/s range rate compensated replica');

end




%% Subfunctions


% initialization of required signal processing parameters
function sp = initialize_prn_sp(sp)
    
    %{
        required inputs:
    
            max_sampling_rate
            tChip
            nRegisters
            rw_length
            fo
    %}
    
    % constant: speed of light in vacuum (m/sec)
    c = 299792458.0;
    
    % up-sampling factor
    upsamp_fac = floor(sp.tChip * sp.max_sampling_rate);
    
	% pseudo-random noise code sequence
    sequence = get_prn_mls(sp.nRegisters);
    
    % sequence phase angle array
    theta(1,:) = pi * sequence;
    
    % up-sampled sequence phase angle array
    theta = reshape(repmat(theta, upsamp_fac, 1), 1, []);
    
    % sampling rate (Hz)
    sampling_rate = upsamp_fac / sp.tChip;
    
    % sampling time (sec)
    tSamp = 1 / sampling_rate;
    
    % number of chips in pulse
    nChips = 2 ^ sp.nRegisters - 1;
    
    % number of samples in transmitted pulse
    tx_samps = upsamp_fac * nChips;
    
    % maximum time delay for target in the range window
    max_delay = 2 * sp.rw_length / c;
    
    % number of receive samples
    rx_samps = ceil(max_delay / tSamp) + tx_samps - 1;
    
    % receive sample time
    rx_t = 0 : tSamp : tSamp * (rx_samps - 1);
    
    % FFT size
    nfft = pow2(nextpow2(rx_samps));
    
    % range window samples
    rw_samps = ceil(max_delay / tSamp);
    
    % range in window
    rng = 0.5 * c * (0 : tSamp : tSamp * (rw_samps - 1));
    
    % display some data to the screen
    fprintf('\n')
    fprintf('     range resolution (m): %.4f\n', 0.5 * c * sp.tChip);
    fprintf('           bandwidth (Hz): %.4e\n', 1 / sp.tChip);
    fprintf('          pulse width (s): %.4e\n\n', tSamp * tx_samps);
    
    fprintf(' number of chips in pulse: %d\n', nChips);
    fprintf('       up-sampling factor: %d\n', upsamp_fac);
    fprintf('         transmit samples: %d\n', tx_samps);
    fprintf('          receive samples: %d\n', rx_samps);
    fprintf('                 FFT size: %d\n\n', nfft);
    
    % required outputs
    sp.c = c;
    sp.upsamp_fac = upsamp_fac;
    sp.theta = theta;
    sp.tSamp = tSamp;
    sp.tx_samps = tx_samps;
    sp.rx_samps = rx_samps;
    sp.rx_t = rx_t;
    sp.nfft = nfft;
    sp.rw_samps = rw_samps;
    sp.rng = rng;
    
end


% generate a pseudo-random noise maximal length sequence
function sequence = get_prn_mls(N)
    
    % limit register size
    max_N = 10;
    if N < 2 || N > max_N
        error('Number of registers outside limits');
    end
    
    % get the feedback taps
    taps = get_mls_polynomials(max_N);
    
    % initialize the linear feedback shift register
    register = mod(1:N, 2);
    
    % initialize the sequence array
    n_sequence = 2^N - 1;
    sequence = zeros(1, n_sequence);
    
    % generate the sequence
    for i = 1 : n_sequence
        
        % calculate the feedback using modulo-2 addition on the taps
        feedback = mod(sum(register(taps{N})), 2);
        
        % shift register contents one register to the right
        register = circshift(register, 1);
        
        % store feedback into the first bit
        register(1) = feedback;
        
        % update the output sequence with the contents of the last register
        sequence(i) = register(N);
    end
    
end


% table of maximal-length polynomials for shift register length (up to 24)
% source:  https://en.wikipedia.org/wiki/Linear-feedback_shift_register
function out = get_mls_polynomials(N)

    fp = cell(24,1);
    
    % feedback polynomial           period (2^N - 1)
    fp{2} = [1, 2];                 % 3
    fp{3} = [2, 3];                 % 7
    fp{4} = [3, 4];                 % 15
    fp{5} = [3, 5];                 % 31
    fp{6} = [5, 6];                 % 63
    fp{7} = [6, 7];                 % 127
    fp{8} = [4, 5, 6, 8];           % 255
    fp{9} = [5, 9];                 % 512
    fp{10} = [7, 10];               % 1,023
    fp{11} = [9, 11];               % 2,047
    fp{12} = [4, 10, 11, 12];       % 4,095
    fp{13} = [8, 11, 12, 13];       % 8,191
    fp{14} = [2, 12, 13, 14];       % 16,383
    fp{15} = [14, 15];              % 32,767
    fp{16} = [4, 13, 15, 16];       % 65,535
    fp{17} = [14, 17];              % 131,071
    fp{18} = [11, 18];              % 262,143
    fp{19} = [14, 17, 18, 19];      % 524,287
    fp{20} = [17, 20];              % 1,048,575
    fp{21} = [19, 21];              % 2,097,151
    fp{22} = [21, 22];              % 4,194,303
    fp{23} = [18, 23];              % 8,388,607
    fp{24} = [17, 22, 23, 24];      % 16,777,215
    
    out = fp(1:N, 1);

end


% compute doppler compensated replica spectrum for PRN waveform
function H_f = get_prn_replica_spectrum(sp, vel)
    
    % compute Doppler-induced pulse dilation parameters
	D = dilation_parameters(sp, vel);
    
    % compute the sampling times for the dilated signal
    t = 0 : sp.tSamp : sp.tSamp * floor(D.pw / sp.tSamp);
    
    % interpolate the coded waveform phase angle at the sample times from
    % the time dilated signal
    phi = interp1(0 : D.tSamp : D.pw, sp.theta, t, 'previous', 'extrap');
    
    % compute the optimum filter impulse response for the dilated signal
    h_t = D.amp_scale * exp(1j * phi) .* exp(1j * 2 * pi * D.fd * t);
    
    % compute the windowing weights
    nSamps = numel(h_t);
    %weights = 0.54 - 0.46*cos(2*pi*(1 : nSamps)/(nSamps-1));
    weights = ones(1, nSamps);
    
    % generate the optimum filter transfer function for the dilated signal
    H_f = conj(fft(weights .* h_t / sp.nfft, sp.nfft));
    
end


% % compute Doppler-induced pulse dilation parameters
function D = dilation_parameters(sp, vel)

    % speed of light in vacuum (m/sec)
    c = 299792458.0;
    
    % temporal dilation factor
    beta = (c - vel) / (c + vel);
    
    % Doppler frequency
    D.fd = sp.fo * (beta - 1);
    
    % amplitude scaling factor
    D.amp_scale = sqrt(beta);
    
    % dilated chip width
    D.tChip = sp.tChip / beta;
    
    % dilated sampling interval
    D.tSamp = D.tChip / sp.upsamp_fac;
    
    % dilated pulse width
    D.pw = D.tSamp * (sp.tx_samps - 1);
    
end


% generate the target return signal for the transmitted PRN waveform
function X_f = get_prn_return_signal_spectrum(tgt, sp)
    
    % initialize the input signal array (i.e., the target return)
    x_t = zeros(1, sp.rx_samps);
    
    % loop through the scatterers...
    for i = 1 : numel(tgt.rng)
        
        % compute Doppler-induced pulse dilation parameters
        D = dilation_parameters(sp, tgt.vel(i));
        
        % compute the time delay to the target in the range window
        tDelay = 2 * tgt.rng(i) / sp.c;
        
        % determine the index for first ADC sample time for this scatterer
        start_idx = find(sp.rx_t >= tDelay, 1, 'first');
        
        % compute the time delta between the target return time delay and
        % the first ADC sample time
        delta_t = sp.rx_t(start_idx) -  tDelay;
        
        % determine the index for last ADC sample time for this scatterer
        stop_idx = find(sp.rx_t < tDelay + D.pw - delta_t, 1, 'last');
        
        % specify the indices for all ADC sample times for this scatterer
        tIdx = start_idx : stop_idx;
        
        % get the array of ADC sample times for this scatterer
        sample_times = sp.rx_t(tIdx);
        
        % compute the time array required to generate the received signal
        t = delta_t + (sample_times - sample_times(1));
        
        % interpolate the coded waveform phase angle at the sample times
        % from the time delayed dilated received signal
        phi = interp1(0 : D.tSamp : D.pw, sp.theta, t, 'previous');
        
        % compute the return signal for this scatterer
        r_t = D.amp_scale * exp(1j * phi) .* exp(1j * 2 * pi * D.fd * t);
        
        % update the composite return signal
        x_t(tIdx) = x_t(tIdx) + (sp.nfft / numel(tIdx) * 10^(tgt.rcs(i) / 20) * r_t);
        
    end
    
    % generate receiver noise and add it to the input signal
    n_t = sqrt(sp.rx_samps/2) * 10^(tgt.mds/20) * complex(randn(1, sp.rx_samps), randn(1, sp.rx_samps));
    x_t = x_t + n_t;
    
    % compute the received signal spectrum
    X_f = fft(x_t, sp.nfft);
    
end
