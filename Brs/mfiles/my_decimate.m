function [y,fs] = my_decimate(x,decf,fs,filter_order)
%%MY_DECIMATE  [y,fs] = my_decimate(x,decf,fs,filter_order)

if nargin<=3,
    filter_order = 4;
end
decf = round(decf);

if nargin>=3,
    fs = fs/decf;
else
    fs = nan;
end


% low-pass filter
[b,a] = butter(filter_order,0.9/decf,'low');
y = filtfilt(b,a,x);

% downsample signal
y = y(1:decf:end);
