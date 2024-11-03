% Generate pulse train
fc = 2000;                                              % Carrier frequency
fs = 15000;                                             % sampling frequency
N = 432;                                                % number of bits
const = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % Constellation 1 - QPSK/4-QAM
M = length(const);                                      % Number of symbols in the constellation
bpsymb = log2(M);                                       % Number of bits per symbol
Rs = 750;                                               % Symbol rate [symb/s]
Ts = 1/Rs;                                              % Symbol time [s/symb]
fsfd = fs/Rs;                                           % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]
bits = randsrc(1,N,[0 1]);                              % Information bits
m_buffer = buffer(bits, bpsymb)';                       % Group bits into bits per symbol
m = bi2de(m_buffer, 'left-msb')'+1;                     % Bits to symbol index
D = const(m);                                           % Look up symbols using the indices
span = 6;                                               % Set span = 6
t_vec = -span*Ts: 1/fs :span*Ts;                        % create time vector for one sinc pulse
[pulse, t] = rtrcpuls(1,Ts,fs,span);                    % Create root raised-cosine pulse

% Generate preamble
barker = [1,1,1,-1,-1,1,-1];                % Barker code with N = 7
PA = barker + 1i*barker;                    % Send sequence in both I and Q channels
PA_train = conv(pulse,upsample(PA,fsfd));

% Create Tx signal
symbols = [PA,D];
pulse_train = conv(pulse,upsample(symbols,fsfd)); % Upsample and pulse-shape

% Modulate pulse-train onto carrier wave
wc = 2*pi*fc;
t = (1:length(pulse_train))*(1/fs);
s_tx = real(pulse_train*sqrt(2).*exp(1i*wc*t)); % I(t)*cos(wc*t)-Q(t)*sin(wc*t)
s_tx = s_tx/max(abs(s_tx)); % Avoid clipping
figure(1)
plot(t,s_tx)
title('Transmitted signal')

% Add random delay
delay = randi(50,1)
s_tx = [zeros(1,delay),s_tx];

% Phase 
phi=randi(360,1)*1  % random phase between 0 and 360 degrees
s_tx=s_tx*exp(j*phi/180*pi);

% Noise
SNR_dB = 15; % Desired SNR
s_tx = awgn(s_tx,SNR_dB); 
figure(2)
plot(1:length(s_tx),s_tx)
title('Transmitted signal w/ added Gaussian noise')

% Down-convert signal and check preamble
t = (1:length(s_tx))*(1/fs);
s_rx = s_tx*sqrt(2).*exp(1i*wc*t-1);
I = real(s_rx);
Q = imag(s_rx);
I_filt = lowpass(I,0.5*fc,fs);
Q_filt = lowpass(Q,0.5*fc,fs);
S = I_filt - 1i*Q_filt;
corr = matched_filter(S,PA_train);  % Correlate with preamble
figure(3)
plot(real(corr))
title('PA correlation')

% Detect PA and estimate delay/phase
[peak,idx] = max(abs(corr));
P = length(PA);
delay_hat = idx - fsfd*(P+1);
phi_hat = mod(angle(corr(idx))*(180/pi),360);


% Matched filter 
D = length(D);
data_start = delay_hat + 1 + fsfd*(P+1);
data_ind = data_start:fsfd:(data_start+(D-1)*fsfd);
%fulse_data_ind = (data_start+50):fsfd:(data_start+50+(D-1)*fsfd); %test
%constellation
S_mf = matched_filter(S,pulse);
rx_vec = S_mf(data_ind); % Get sample points
rx_vec = rx_vec./mean(abs(rx_vec));
figure(4)
plot(real(S_mf))
title('Matched filter output (real)')

% Phase compensation
rx_vec = rx_vec*exp(-1j*phi_hat/180*pi);

% Plot received symbols
scatterplot(rx_vec)

% Minimum Eucledian distance detector
rx_vec = rx_vec/max(abs(rx_vec));
metric = abs(repmat(rx_vec.',1,4) - repmat(const,length(rx_vec),1)).^2; % Compute the distance to each possible symbol
[tmp, m_hat] = min(metric, [], 2); % Find the closest for each received symbol
m_hat = m_hat'-1; % Get the index of the symbol in the constellation
SER = sum(m-1 ~= m_hat) % Count symbol errors
m_hat = de2bi(m_hat, 2, 'left-msb')'; % Make symbols into bits
a_hat = m_hat(:)'; % Write as a vector
BER = sum(bits ~= a_hat) % Count number of bit errors