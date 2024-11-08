%%%%%%%%%%%%%%%%%%%%%
%%%%% Set up %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
fs = 15000;                                              % sampling frequency
fc = 2000;
N = 432;                                                % number of bits
const = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % Constellation 1 - QPSK/4-QAM
M = length(const);                                      % Number of symbols in the constellation
bpsymb = log2(M);                                       % Number of bits per symbol
Rs = 500;                                               % Symbol rate [symb/s]
Ts = 1/Rs;                                              % Symbol time [s/symb]
fsfd = fs/Rs;                                           % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]

%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Transmitter%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% bits to symbol
bits = randsrc(1,N,[0 1]);                              % Information bits(generate 432 bits)
m_buffer = buffer(bits, bpsymb)';                       % Group bits into bits per symbol
m = bi2de(m_buffer, 'left-msb')'+1;                     % Bits to symbol index
x = const(m);                                           % Look up symbols using the indices
x_i = real(x);
x_q = imag(x);
% Plot
figure(1);
plot (x_i, x_q,' *');
xlabel ('Real');
ylabel ('Imag' );
title('Consetellation(Tx)')
%x = awgn(x,15);                                          % add artificial noise

%%% RC Pulse
span = 6;                                               % Set span = 6
t_vec = -span*Ts: 1/fs :span*Ts;                        % create time vector for one sinc pulse
alpha = 0.5;
pulse = (sin(pi*t_vec/Ts)./(pi*t_vec/Ts)) .* (cos(pi*alpha*t_vec/Ts)./(1-(2*alpha*t_vec/Ts).^2));

% 处理 t_vec = 0 的情况（避免 0 除错误）
pulse(t_vec == 0) = 1;  % 解决 sin(0)/0 的问题
pulse(abs(2*alpha*t_vec/Ts) == 1) = (alpha/2) * sin(pi/(2*alpha));  % 处理分母为0的点

x_i_upsample = upsample(x_i, fsfd);   
x_q_upsample = upsample(x_q, fsfd);   
Sig_I = conv(x_i_upsample,pulse, 'same');
Sig_Q = conv(x_q_upsample,pulse, 'same');

% plot
figure(2);
subplot(2,3,1);plot(Sig_I(1:fsfd*20),'LineWidth',2);title(['Baseband Signal I']);
subplot(2,3,4); plot(Sig_Q(1:fsfd*20),'LineWidth',2);title(['Baseband Signal Q']);

%%% modulate
t = (1 : length(Sig_Q)) / fs;
carrier1 = cos(2 * pi * fc * t);    % 同相载波
carrier2 = sin(2 * pi * fc * t);    % 正交载波
Pass_Sig_I = Sig_I .* carrier1;
Pass_Sig_Q = Sig_Q .* -carrier2;
Sig_qpsk = Pass_Sig_I+Pass_Sig_Q;
subplot(2,3,2); plot(Pass_Sig_I(1:fsfd*20));grid on; title(['Passband I'])
subplot(2,3,5); plot(Pass_Sig_Q(1:fsfd*20));grid on; title(['Passband Q'])
subplot(2,3,3); plot(Sig_qpsk(1:fsfd*20));grid on; title(['Modulated Signal'])

%%% AWGN channel
EbN0 = 12;
snr = EbN0 - 10 * log10(0.5 * fsfd) - 10 * log10(M);
Sig_passband = awgn(Sig_qpsk, snr, 'measured');
subplot(2,3,6);plot(Sig_passband(1:fsfd*20));grid on; title(['Passband Signal(AWGN)'])


%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Receiver%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%


