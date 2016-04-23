function [mx,mt,obp] = cpt_bba2x(x,p,fs_in,fs_out)
%--------------------------------------------
%
% (c) Computational Physiology Toolbox
%
% Beat-to-Beat Average
%
%--------------------------------------------
% 
% 18 June 2015
%
% Adam Mahdi
% Dragana Nikolic

% ************************************************************************
% LPF at 20 Hz by 8th-order Butterworth zero-phase filter
[b,a] = butter(4,20/(fs_in/2),'low');
x     = filtfilt(b,a,x);
p     = filtfilt(b,a,p);

% calculate beat-to-beat average
% [obp,obpt] = wabp(p,0,1,fs_in);
[~,idia]    = cpt_abpTimes4x(p,fs_in,0.9);
obp         = idia;       % in samples
obpt        = obp/fs_in;  % in seconds 

% check if all markers are unique
[~,i,~]     = unique(obp);
obp         = obp(i);
obpt        = obpt(i);
[~,i]       = sort(obp);
obp         = obp(i);
obpt        = obpt(i);

% calculate beat-to-beat average
mx  = zeros(numel(obp)-1,1);
for i=1:numel(obp)-1,
    mx(i,1) = mean(x(obp(i):obp(i+1)));
end

% interpolate beat-to-beat average
mt = 0:1/fs_out:(length(x)-1)/fs_in;
obpt = obpt(2:end);
if obpt(1)~=0, 
    obpt = [0;obpt];
    mx = [mx(1);mx];
end
if obpt(end)~=(length(x)-1)/fs_in, 
    obpt = [obpt;(length(x)-1)/fs_in];
    mx = [mx;mx(end)];
end   
mx = interp1(obpt,mx,mt,'linear','extrap')';
