clear all;close all; clc 

%% ADD PATHS & READ DATA

addpath('mfiles/')
addpath('data/mat_Y')
datanames_Y;
load(strcat('mat_Y/',dat_Y{7}))

% t  = t(1:3000);
% bp = bp(1:3000);


%% EXCTRACT SBP AND RR
fs_in = round(1/mean(diff(t)));     % current sampling frequenc

% sys and dia indices
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);
%rri     = diff([-Inf;id])/fs_in; 
rri     = diff([-Inf;t(id)]);

%% PARAMETERS FOR THE SEQUENCE METHOD

no1         = 3;                    % # successive points
dsys_thresh = 1;                    % in mmHg
drri_thresh = 0.005;                % in seconds
lag_inc     = 1;                    % correction for delay RR (increase)

params = [fs_in, no1, dsys_thresh, drri_thresh, lag_inc];

%%  INCREASE in SBP and RR
incS = diff([-Inf;bp(is)])>=dsys_thresh;
incT = diff([-Inf;rri])>=drri_thresh;

% Correction for the delay between SYS and RR (increase)
%incS  = [zeros(lag_inc,1);incS];
%incT  = [incT;zeros(lag_inc,1)];

% detect sequences with no1=3 or more successive increases in SYS 
strS                = num2str(incS');
strS(isspace(strS)) = '';
isi_b               = strfind(strS,['0' repmat('1',1,no1)]);
isi_e               = strfind(strS,[repmat('1',1,no1) '0'])+no1-1;
if length(isi_e)<length(isi_b), isi_e(end+1) = length(incS); end

maskS = zeros(length(is),1);
for i=1:length(isi_b),
    maskS(isi_b(i):isi_e(i)) = 1;
end

% detect sequences with no1=3 or more successive increases in RRI 
strT                = num2str(incT');
strT(isspace(strT)) = '';
iti_b               = strfind(strT,['0' repmat('1',1,no1)]);
iti_e               = strfind(strT,[repmat('1',1,no1) '0'])+no1-1;
if length(iti_e)<length(iti_b), iti_e(end+1) = length(incT); end

maskT = zeros(length(is),1);
for i=1:length(iti_b),
    maskT(iti_b(i):iti_e(i)) = 1;
end

% Correction for the delay between SYS and RR (increase)
maskS_inclag  = [zeros(lag_inc,1);maskS];
maskT_inclag  = [maskT;zeros(lag_inc,1)];

% Simultaneous increase/decrease of BRS and RR
maskST_inclag = maskS_inclag.*maskT_inclag;

% Re-correction
maskSx = maskST_inclag(lag_inc+1:end);
maskTx = maskST_inclag(1:end-lag_inc);


%%  DECREASE in SBP and RR
incS_d = -diff([-Inf;bp(is)])>=dsys_thresh;
incT_d = -diff([-Inf;rri])>=drri_thresh;


% detect sequences with no1=3 or more successive decrease in SBP 
strS_d                  = num2str(incS_d');
strS_d(isspace(strS_d)) = '';
isi_b_d                 = strfind(strS_d,['0' repmat('1',1,no1)]);
isi_e_d                 = strfind(strS_d,[repmat('1',1,no1) '0'])+no1-1;
if length(isi_e_d)<length(isi_b_d), isi_e_d(end+1) = length(incS_d); end

%[isi_b_d' isi_e_d']

maskS_d = zeros(length(is),1);
for i=1:length(isi_b_d),
    maskS_d(isi_b_d(i):isi_e_d(i)) = 1;
end


% detect sequences with no1=3 or more successive decrese in RRI 
strT_d                  = num2str(incT_d');
strT_d(isspace(strT_d)) = '';
iti_b_d                 = strfind(strT_d,['0' repmat('1',1,no1)]);
iti_e_d                 = strfind(strT_d,[repmat('1',1,no1) '0'])+no1-1;
if length(iti_e_d)<length(iti_b_d), iti_e_d(end+1) = length(incT_d); end

maskT_d = zeros(length(is),1);
for i=1:length(iti_b_d),
    maskT_d(iti_b_d(i):iti_e_d(i)) = 1;
end

% Correction for the delay between SYS and RR (decrease)

maskS_declag  = [zeros(lag_inc,1);maskS_d];
maskT_declag  = [maskT_d;zeros(lag_inc,1)];

% Simultaneous increase/decrease of BRS and RR
maskST_declag = maskS_declag.*maskT_declag;

% Re-correction
maskSx_d = maskST_declag(lag_inc+1:end);
maskTx_d = maskST_declag(1:end-lag_inc);


%% Fitting 

dSBP     = bp(is(find(maskSx)));
dRR      = rri(find(maskTx));

dSBP_d   = bp(is(find(maskSx_d)));
dRR_d    = rri(find(maskTx_d));

xx = [dRR; dRR_d];
yy = [dSBP; dSBP_d];

p  =polyfit(xx,yy,1)
vv = polyval(p,xx);

figure, hold on
plot(xx,yy,'o')
plot(xx,vv)

%return;

%% PLOTTING
% [isi_b(:) isi_e(:)]
% [iti_b(:) iti_e(:)]


% figure, hold on, 
% %stem(1.4*maskS.*maskT,'m')
% stem(2*maskST_inclag)
% stem(1.2*maskT,'r')
% stem(maskS)


figure, hold on,
subplot(211),hold on
    plot(t,bp)
    plot(t(is(find(maskS))),bp(is(find(maskS))),'.k', 'MarkerSize',20)
    plot(t(is(find(maskS_d))),bp(is(find(maskS_d))),'ro')
subplot(212),hold on
    plot(rri)
    plot(find(maskT),rri(find(maskT)),'.k', 'MarkerSize',20)
    plot(find(maskT_d),rri(find(maskT_d)),'ro')

    
% figure, hold on, 
% stem(2*maskST_declag)
% stem(1.2*maskT_d,'r')
% stem(maskS_d)

% figure, hold on,
% subplot(211),hold on
% plot(t,bp)
% %plot(t(is(find(maskS))),bp(is(find(maskS))),'ro')
% plot(t(is(find(maskSx_d))),bp(is(find(maskSx_d))),'ro')
% subplot(212),hold on
% plot(rri)
% %plot(find(maskT),rri(find(maskT)),'ro')
% plot(find(maskTx_d),rri(find(maskTx_d)),'ro')