clear all;close all; clc 


%% ADD PATHS & READ DATA

addpath('mfiles/')
addpath('data/mat_Y')
datanames_Y;
load(strcat('mat_Y/',dat_Y{7}))

% t  = t(1:3000);
% bp = bp(1:3000);


%% CALCULATE BRS
fs_in = round(1/mean(diff(t)));       % current sampling frequenc

% sys and dia indices
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);

rri = diff([-Inf;is])/fs_in;


no1 = 3;                % number of successive points
dsys_thresh = 1;        % in mmHg
drri_thresh = 0.005;	% in seconds

% increase sequence of sys
incS = diff([-Inf;bp(is)])>=dsys_thresh;
incT = diff([-Inf;rri])>=drri_thresh;


% detect sequences with no1=3 or more successive increases in SP 
strS = num2str(incS');
strS(isspace(strS)) = '';
isi_b = strfind(strS,['0' repmat('1',1,no1)]);
isi_e = strfind(strS,[repmat('1',1,no1) '0'])+2;
if length(isi_e)<length(isi_b), isi_e(end+1) = length(incS); end

maskS = zeros(length(is),1);
for i=1:length(isi_b),
    maskS(isi_b(i):isi_e(i)) = 1;
end

% detect sequences with no1=3 or more successive increases in RRI 
strT = num2str(incT');
strT(isspace(strT)) = '';
iti_b = strfind(strT,['0' repmat('1',1,no1)]);
iti_e = strfind(strT,[repmat('1',1,no1) '0'])+2;
if length(iti_e)<length(iti_b), iti_e(end+1) = length(incT); end

maskT = zeros(length(is),1);
for i=1:length(iti_b),
    maskT(iti_b(i):iti_e(i)) = 1;
end


[isi_b(:) isi_e(:)]
[iti_b(:) iti_e(:)]

figure, hold on,
subplot(211),hold on
plot(t,bp)
plot(t(is(find(maskS))),bp(is(find(maskS))),'ro')
subplot(212),hold on
plot(rri)
plot(find(maskT),rri(find(maskT)),'ro')

figure, hold on, 
stem(1.4*maskS.*maskT,'m')
stem(1.2*maskT,'r')
stem(maskS)


return


% minimum correlation value between systolic pressure and RRI 0.8

% monotonic sequence of size >=4.
incx = inc1;
incy = false(1,length(incx))';
for j=1:length(incx)-2,
    if incx(j)==1,
        if incx(j+1)==1,
            if incx(j+2)==1,
               % incy(j-1) = 1;
               incy(j)   = 1;
               incy(j+1) = 1;
               incy(j+2) = 1;
            end
        end
    end
end

inc_seq1 = is(inc1);
inc_seq2 = is(inc2);    
inc_seq_sum = union(inc_seq1,inc_seq2);
inc_seqy = is(incy); 

% choose type of sequence
inc_seq = inc_seqy;

% calculate something
per = diff([t(id);0;0]);



%% PLOTTING

figure, hold on
plot(t,bp)
plot(t(id),bp(id),'ro')
% plot(t(is),bp(is),'r+')
plot(t(inc_seq),bp(inc_seq),'k+')

figure, hold on
subplot(211)
plot(bp(is(1:end-4)))
subplot(212)
plot(per(1:end-2))


figure, hold on
scatter(bp(is(incy(1:end-2))),per(incy(2:end-1)))
xlim([60 200])
ylim([0.3 1])


return
%% PSD

% psd for sys
Psys = bp(is(1:end-3));
tsys = t(is(1:end-3));
[Pow_sys,fre_sys]   = plomb(Psys,tsys);     % using LOMB method

yperx = diff(bp(id));
yper = yperx(1:end-1);
tper = t(id(1:end-2));
[Pow_per,fre_per]   = plomb(yper,tper);     % using LOMB method


figure
subplot(211), hold on
plot(fre_sys,Pow_sys)
plot(fre_per,Pow_per)
xlim([0.09,0.3])
subplot(212)
plot(fre_per,Pow_per)
xlim([0.09,0.3])

%% Coherence

C_t = 0.5;
[Cxy,F] = mscohere(Psys,yper);    % find the coherence function

% plot coherence
figure;
plot(F,Cxy);

% % find frequencies of good coherence, and truncate
% F_good = find((Cxy > C_t) & (F < fH/(2*pi)));
% hold on
% plot(F(F_good),Cxy(F_good),'ro');
