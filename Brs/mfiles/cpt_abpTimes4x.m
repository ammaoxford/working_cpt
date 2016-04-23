function [ind_systol,ind_diastol] = cpt_abpTimes4x(bp,fs_in,sep)
%--------------------------------------------
%
% Computational Physiology Toolbox
%
% Adam Mahdi, University of Oxford
%
%--------------------------------------------
%
% Sept 2015

% ************************************************************************
% LPF at 20 Hz by 8th-order Butterworth zero-phase filter
[b,a]           = butter(4,20/(fs_in/2),'low');
bp              = filtfilt(b,a,bp);

% ************************************************************************
% COMPUTE SYSTOLIC INDICES
my_mpk          = round(sep*fs_in);                         % define a separation, choose +/- 0.6 sec.
[~,ind_s]       = findpeaks(bp,'MinPeakDistance',my_mpk);   % initial peak detection
ff              = 60*fs_in/round(mean(diff(ind_s)));        % estimate HR
data            = bp';
bsline          = LPFilter(data,.7/fs_in);                  % base wander removal
peaks           = PeakDetection(data-bsline,ff/fs_in);      % peak detection
ind_systol      = find(peaks);

% ************************************************************************
% COMPUTE DIASTOLIC INDICES
wind_dias       = round(0.25*(fs_in));
ST              = ind_systol(1:end-1);
BeatsNo         = length(ST);

MinDomain  = zeros(BeatsNo,wind_dias);
for i=1:wind_dias
    MinDomain(:,i) = ST-i+1;
end

MinDomain(MinDomain<1) = 1;  % Error protection
[~,Dindex]             = min(bp(MinDomain),[],2);
ind_diastol            = MinDomain(sub2ind(size(MinDomain),(1:BeatsNo)',Dindex));
%--------------------

wind_dn    = round(0.4*(fs_in));
MaxDomain  = zeros(BeatsNo,wind_dn);
for i=1:wind_dn,
    MaxDomain(:,i)  = ST+i+1;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Dicrotic Dip (dd) & Dicrotic Notch (dn)
% my_mpk      = .6;  % define a separation
% for jj=1:BeatsNo
% [pks_dd,obpx_dd]     = findpeaks(-bp(MaxDomain(jj,:)),'MinPeakDistance',my_mpk);
% [pks_dn,obpx_dn]     = findpeaks( bp(MaxDomain(jj,:)),'MinPeakDistance',my_mpk);
%
%     if length(obpx_dn)>1 % &&length(obpx_dn)>1
%         for i=1:2        % length(obpx_dn)
%         mojtest(i,1) = pks_dn(i,1)+pks_dd(i,1); % remember pks_dd is negative
%         [~,xi]       = max(mojtest(:,1));
%     DDt_index(jj,1)  = obpx_dd(xi,1);
%     DNt_index(jj,1)  = obpx_dn(xi,1);
%         end
%     else
%     DDt_index(jj,1)  = obpx_dd(1,1);
%     DNt_index(jj,1)  = obpx_dn(1,1);
%     end
% end
% ind_dd    = MaxDomain(sub2ind(size(MaxDomain),(1:BeatsNo)',DDt_index));
% ind_dn    = MaxDomain(sub2ind(size(MaxDomain),(1:BeatsNo)',DNt_index));
% %--------------------



function peaks = PeakDetection(x,ff)
%
% peaks = PeakDetection(x,f),
% R-peak detector
%
% inputs:
% x: vector of input data
% f: approximate ECG beat-rate in Hertz, normalized by the sampling frequency
%
% output:
% peaks: vector of R-peak impulse train
%
% Notes:
% - The R-peaks are found from a peak search in windows of length N; where
% N corresponds to the R-peak period calculated from the given f. R-peaks
% with periods smaller than N/2 or greater than N are not detected.
% - The signal baseline wander is recommended to be removed before the
% R-peak detection
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.

N = length(x);
peaks = zeros(N,1);

%th = .5; % Original threshold
th = 20;  % My change (Adam)
rng = floor(th/ff);

if (abs(max(x))>abs(min(x))),
    for j = 1:N,
        % index = max(j-rng,1):min(j+rng,N);
        if (j>rng && j<N-rng),
            index = j-rng:j+rng;
        elseif (j>rng),
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        
        if max(x(index))==x(j),
            peaks(j) = 1;
        end
    end
else
    for j = 1:N,
        % index = max(j-rng,1):min(j+rng,N);
        if (j>rng && j<N-rng),
            index = j-rng:j+rng;
        elseif (j>rng),
            index = N-2*rng:N;
        else
            index = 1:2*rng;
        end
        
        if min(x(index))==x(j),
            peaks(j) = 1;
        end
    end
end


% remove fake peaks
I = find(peaks);
d = diff(I);
z = find(d<rng);
peaks(I(z))=0;

function y = LPFilter(x,fc)
%
% y = LPFilter(x,fc),
% First order zero-phase Lowpass filter
%
% inputs:
% x: vector or matrix of input data (channels x samples)
% fc: -3dB cut-off frequency normalized by the sampling frequency
%
% output:
% y: vector or matrix of filtered data (channels x samples)
%
%
% Open Source ECG Toolbox, version 1.0, November 2006
% Released under the GNU General Public License
% Copyright (C) 2006  Reza Sameni
% Sharif University of Technology, Tehran, Iran -- LIS-INPG, Grenoble, France
% reza.sameni@gmail.com

% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details. You should have received a copy of the
% GNU General Public License along with this program; if not, write to the
% Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
% MA  02110-1301, USA.


k = .7;         % cut-off value
alpha = (1-k*cos(2*pi*fc)-sqrt(2*k*(1-cos(2*pi*fc))-k^2*sin(2*pi*fc)^2))/(1-k);
y = zeros(size(x));
for i = 1:size(x,1),
    y(i,:) = filtfilt(1-alpha,[1 -alpha],x(i,:));
end
