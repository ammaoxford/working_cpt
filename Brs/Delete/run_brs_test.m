clear all; close all; clc 


%% ADD PATHS & READ DATA

addpath('mfiles/')
addpath('data/mat_Y')
datanames_Y;
load(strcat('mat_Y/',dat_Y{7}))


% initialisation
no1 = 3;                % number of successive points
dsys_thresh = 1;        % in mmHg
drri_thresh = 0.005;	% in seconds



%% CALCULATE BRS
fs_in = round(1/mean(diff(t)));       % current sampling frequenc

% sys and dia indices
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);
rri = diff([-Inf;id])/fs_in;

% increase sequence of sys
incS = diff([-Inf;bp(is)])>=dsys_thresh;
incT = diff([-Inf;rri])>=drri_thresh;

% detect sequences with no1=3 or more successive increases in SP 
strS = num2str(incS');
strS(isspace(strS)) = '';
isi_b = strfind(strS,['0' repmat('1',1,no1)]);
isi_e = strfind(strS,[repmat('1',1,no1) '0'])+2;
if isi_e(1)<isi_b(1), isi_b = [1 isi_b]; end
if isi_e(end)<isi_b(end), isi_e = [isi_e length(incS)]; end

maskS = zeros(length(is),1);
for i=1:length(isi_b),
    maskS(isi_b(i):isi_e(i)) = 1;
end

% detect sequences with no1=3 or more successive increases in RRI 
strT = num2str(incT');
strT(isspace(strT)) = '';
iti_b = strfind(strT,['0' repmat('1',1,no1)]);
iti_e = strfind(strT,[repmat('1',1,no1) '0'])+2;
if iti_e(1)<iti_b(1), iti_b = [1 iti_b]; end
if iti_e(end)<iti_b(end), iti_e = [iti_e length(incT)]; end

maskT = zeros(length(is),1);
for i=1:length(iti_b),
    maskT(iti_b(i):iti_e(i)) = 1;
end

mask = maskS & maskT;

[isi_b(:) isi_e(:)]
[iti_b(:) iti_e(:)]

figure, hold on,
subplot(211),hold on
plot(t,bp)
plot(t(is(find(maskS))),bp(is(find(maskS))),'ro')
plot(t(is(mask)),bp(is(mask)),'k*')
xlabel('time [s]')
subplot(212),hold on
plot(rri)
plot(find(maskT),rri(find(maskT)),'ro')
plot(find(mask),rri(mask),'k*')
xlabel('samples')

figure, hold on, 
stem(1.4*mask,'m')
stem(1.2*maskT,'r')
stem(maskS)



