
%%%% TRANSMITTER %%%%

% pack = message to be transmitted (consists of 432 bits from the GUI, always!)
% fc = carrier frequency
function transmitter(pack,fc)
    
    %-----------------------------------------------------------
    % Set up:
    %-----------------------------------------------------------
    wc = 2*pi*fc;                                               % Angular frequency
    fs = 15000;                                                 % Sampling frequency
    const = [(1 + 1i) (1 - 1i) (-1 + 1i) (-1 - 1i)]/sqrt(2);    % QPSK (Gray coded)
    M = length(const);                                          % Number of symbols in the constellation
    bpsymb = log2(M);                                           % Number of bits per symbol
    Rs = 150;                                                   % Symbol rate [symb/s]
    Ts = 1/Rs;                                                  % Symbol time [s/symb]
    fsfd = fs/Rs;                                               % Number of samples per symbol (choose fs such that fsfd is an integer) [samples/symb]
    span = 6;                                                   % Set span = 6
    [pulse, t] = rtrcpuls(1,Ts,fs,span);                        % Create root raised-cosine pulse
    
    %-----------------------------------------------------------
    % Label:
    %-----------------------------------------------------------
    data = pack;
    m_buffer = buffer(data, bpsymb)';   % Group bits into bits per symbol
    m = bi2de(m_buffer, 'left-msb')'+1; % Bits to symbol index
    D = const(m);                       % Look up symbols using the indices
    
    %-----------------------------------------------------------
    % Preamble:
    %-----------------------------------------------------------
    barker = [1,1,1,-1,-1,1,-1];                % Barker code with N = 7
    PA = barker + 1i*barker;                    % Send sequence in both I and Q channels
    
    %-----------------------------------------------------------
    % Create Tx signal
    %-----------------------------------------------------------
    symbols = [PA,D];
    pulse_train = conv(pulse,upsample(symbols,fsfd)); % Upsample and pulse-shape
    
    %-----------------------------------------------------------
    % Modulate pulse-train onto carrier wave
    %-----------------------------------------------------------
    t = (1:length(pulse_train))*(1/fs);
    s_tx = real(pulse_train*sqrt(2).*exp(1i*wc*t)); % I(t)*cos(wc*t)-Q(t)*sin(wc*t)
    s_tx = s_tx/max(abs(s_tx)); % Avoid clipping
    figure(1)
    plot(t,s_tx)
    title('Transmitted signal')
    
    %-----------------------------------------------------------
    % Audioplayer
    %-----------------------------------------------------------
    player = audioplayer(s_tx,fs);
    playblocking(player);
end