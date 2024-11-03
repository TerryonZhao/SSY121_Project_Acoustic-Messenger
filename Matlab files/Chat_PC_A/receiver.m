
%%%% RECEIVER %%%%

function [audio_recorder] = receiver(fc)
    fs = 15000; % Sampling frequency
    audio_recorder = audiorecorder(fs,24,1); % Create the recorder
    
    % Attach callback function
    time_value = 0.5; % How often the function should be called in seconds
    set(audio_recorder,'TimerPeriod',time_value,'TimerFcn',@audioTimerFcn); % Attach a function that should be called every second, the function that is called is specified below.
    
    % ADD USER DATA FOR CALLBACK FUNCTION (DO NOT CHANGE THE NAMES OF THESE VARIABLES!)
    audio_recorder.UserData.receive_complete = 0;   % This is a flag that the while loop in the GUI will check
    audio_recorder.UserData.pack  = [];             % Allocate for data package
    audio_recorder.UserData.pwr_spect = [];         % Allocate for PSD
    audio_recorder.UserData.const = [];             % Allocate for constellation
    audio_recorder.UserData.eyed  = [];             % Allocate for eye diagram
    audio_recorder.UserData.fc = fc;                % Save fc so it can be accessed from callback
    
    record(audio_recorder); % Start recording
end
    
    
% CALLBACK FUNCTION
function audioTimerFcn(recObj, event, handles)
    
    %-----------------------------------------------------------
    % Set up:
    %-----------------------------------------------------------
    disp('Callback triggered')
    
    fc = recObj.UserData.fc;                                    % Carrier frequency
    wc = 2*pi*fc;                                               % Angular frequency
    B_max = 400;                                                % Maximum signal bandwidth (LPF cut-off)
    fs = 15000;                                                 % Sampling frequency
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 - 1i)]/sqrt(2);    % QPSK (Gray coded)
    Rs = 150;                                                   % Symbol rate [symb/s]
    Ts = 1/Rs;                                                  % Symbol time [s/symb]
    fsfd = fs/Rs;                                               % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]
    span = 6;                                                   % Set span = 6
    [pulse, t] = rtrcpuls(1,Ts,fs,span);                        % Create root raised-cosine pulse
    barker = [1,1,1,-1,-1,1,-1]; 
    
    disp('Finish setting')
    
    %-----------------------------------------------------------
    % Receive signal:
    %-----------------------------------------------------------
    s_rx = getaudiodata(recObj)';
    s_rx = s_rx./max(abs(s_rx));
    
    disp('Signal received')
    
    %-----------------------------------------------------------
    % Down-convert signal:
    %-----------------------------------------------------------
    t = (1:length(s_rx))*(1/fs);
    s_rx = s_rx*sqrt(2).*exp(1i*wc*t-1);
    I = real(s_rx);
    Q = imag(s_rx);
    
    disp('Demodulation done')
    
    %-----------------------------------------------------------
    % LPF:
    %-----------------------------------------------------------
    I_filt = lowpass(I,B_max,fs);
    Q_filt = lowpass(Q,B_max,fs);
    S = I_filt - 1i*Q_filt;
    
    disp('Signal passed LPF successfully')
    
    %-----------------------------------------------------------
    % Check preamble and Time sync:
    %-----------------------------------------------------------
    PA = barker + 1i*barker;                    % Preamble is sent in both I and Q channels
    PA_train = conv(pulse,upsample(PA,fsfd));   % Pulse shape
    corr = matched_filter(S,PA_train);          % Correlate with preamble
    
    [peak,idx] = max(abs(corr));
    P = length(PA);
    delay_hat = idx - fsfd*(P+1);
    phi_hat = mod(angle(corr(idx))*(180/pi),360);
    
    % Add a threshold to avoid noise
    % Recognize noise from signal & Check completeness
    PA_threshold = 4;
    
    if (PA_threshold > max(abs(corr)) || ((length(S) - idx) < fsfd * 216))
        disp("Can't find signal...")
        return % End this callback function
    end
    
    %-----------------------------------------------------------
    % MF:
    %-----------------------------------------------------------
    D = 216; % Symbol number
    data_start = delay_hat + 1 + fsfd*(P+1);
    data_ind = data_start:fsfd:(data_start+(D-1)*fsfd);
    
    S_mf = matched_filter(S,pulse);
    rx_vec = S_mf(data_ind); % Get sample points

    disp('Signal passed MF successfully')
    
    %-----------------------------------------------------------
    % Phase compensation
    %-----------------------------------------------------------
    rx_vec = rx_vec*exp(-1j*phi_hat/180*pi);
    
    %-----------------------------------------------------------
    % Minimum Eucledian distance detector
    %-----------------------------------------------------------
    rx_vec = rx_vec/max(abs(rx_vec));
    metric = abs(repmat(rx_vec.',1,4) - repmat(const,length(rx_vec),1)).^2; % Compute the distance to each possible symbol
    [tmp, m_hat] = min(metric, [], 2);                                      % Find the closest for each received symbol
    m_hat = m_hat'-1;                                                       % Get the index of the symbol in the constellation
    m_hatb = de2bi(m_hat, 2, 'left-msb')';                                  % Make symbols into bits
    a_hat = m_hatb(:)';                                                     % Write as a vector
    
    %-----------------------------------------------------------
    % Parameters needed to be set for further compute
    %-----------------------------------------------------------
    bits = a_hat;
    x = rx_vec;
    pulse_train = S_mf;
    
    %------------------------------------------------------------------------------
    % SAVE DATA FOR THE GUI
    %------------------------------------------------------------------------------
    
    % Save the estimated bits
    recObj.UserData.pack = bits;
    
    % Save the sampled symbols
    recObj.UserData.const = x;
    
    % Provide the matched filter output for the eye diagram 
    recObj.UserData.eyed.r = pulse_train((length(pulse)-1)/2:end-(length(pulse)-1)/2+1);
    recObj.UserData.eyed.fsfd = fsfd;
    
    % Compute the PSD and save it. 
    % !!!! NOTE !!!! the PSD should be computed on the BASE BAND signal BEFORE matched filtering
    [pxx, f] = pwelch(S,1024,768,1024, fs);     % Note that pwr_spect.f will be normalized frequencies
    f = fftshift(f);                            % Shift to be centered around fs
    f(1:length(f)/2) = f(1:length(f)/2) - fs;   % Center to be around zero
    p = fftshift(10*log10(pxx/max(pxx)));       % Shift, normalize and convert PSD to dB
    recObj.UserData.pwr_spect.f = f;
    recObj.UserData.pwr_spect.p = p;
    
    % In order to make the GUI look at the data, we need to set the
    % receive_complete flag equal to 1:
    recObj.UserData.receive_complete = 1; 
end
    