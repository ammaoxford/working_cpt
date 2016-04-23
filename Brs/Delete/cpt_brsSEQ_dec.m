function [maskSx, maskTx, maskS, maskT] = cpt_brsSEQ_inc(bp,is,id,rri, params)


%% PARAMETERS

fs_in       = params(1);
no1         = params(2);         % # successive points
dsys_thresh = params(3);         % in mmHg
drri_thresh = params(4);         % in seconds
lag_inc     = params(5);         % correction for delay RR (increase)



%% COMPUTE RR
% % sys and dia indices
% [is,id] = cpt_abpTimes4x(bp,fs_in,0.9);
%rri     = diff([-Inf;is])/fs_in;

       
%%  INCREASE in SBP and RR
incS = -diff([-Inf;bp(is)])>=dsys_thresh;
incT = -diff([-Inf;rri])>=drri_thresh;

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
