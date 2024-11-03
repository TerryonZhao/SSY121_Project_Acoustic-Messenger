% Apply matched filter to complex input signal. Input parameters
% are the received signal and the complex envelope of the 
% transmitted signal.
function y = matched_filter(x_r, x_t)
    h_m = fliplr(conj(x_t));
    y = conv(x_r,h_m);
end