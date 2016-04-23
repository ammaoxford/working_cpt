function [maskSx,maskTx,maskS,maskT,indSx,indTx,indS,indT] = cpt_brsSEQ_incdec(bp,rri,is,id,params,opt)

% This function computes the following indices
% maskSx: lagged & simultaneous indices of BRS increase sequences
% maskTx: lagged & simultaneous indices of RR increase sequences
% maskS:  indices of BRS increase sequences
% maskT:  indices of RR increase sequences
% opt: '1' - increase of brs sequences; '2'- decrease of brs sequences.



%% PARAMETERS
fs_in       = params(1);
no1         = params(2);         % # successive points
dsys_thresh = params(3);         % in mmHg
drri_thresh = params(4);         % in seconds
brs_lag     = params(5);         % correction for delay RR (increase)


       
%%  INCREASE in SBP and RR

% Choosing '1' for increase and '2' for decrease
if opt==1,
    incS = diff([-Inf;bp(is)])>=dsys_thresh;
    incT = diff([-Inf;rri])>=drri_thresh;
elseif opt==2,
    incS = -diff([-Inf;bp(is)])>=dsys_thresh;
    incT = -diff([-Inf;rri])>=drri_thresh;
end

% detect sequences with no1=3 or more successive increases in SBP 
strS                = num2str(incS');
strS(isspace(strS)) = '';
isi_b               = strfind(strS,['0' repmat('1',1,no1)]);
isi_e               = strfind(strS,[repmat('1',1,no1) '0'])+no1-1;
if length(isi_e)<length(isi_b),isi_e(end+1) = length(incS); end
if length(isi_e)>length(isi_b),isi_b = [1 isi_b]; end

maskS = zeros(length(is),1);
for i=1:length(isi_b),
    maskS(isi_b(i):isi_e(i)) = 1;
end
indS = [isi_b(:) isi_e(:)];
clear isi_*

% detect sequences with no1=3 or more successive increases in RRI 
strT                = num2str(incT');
strT(isspace(strT)) = '';
iti_b               = strfind(strT,['0' repmat('1',1,no1)]);
iti_e               = strfind(strT,[repmat('1',1,no1) '0'])+no1-1;
if length(iti_e)<length(iti_b),iti_e(end+1) = length(incT); end
if length(iti_e)>length(iti_b),iti_b = [1 iti_b]; end

maskT = zeros(length(is),1);
for i=1:length(iti_b),
    maskT(iti_b(i):iti_e(i)) = 1;
end
indT = [iti_b(:) iti_e(:)];
clear iti_*

% lag between BRS and RR 
maskS_lag  = [zeros(brs_lag,1);maskS];
maskT_lag  = [maskT;zeros(brs_lag,1)];

% condition for simultaneous increase/decrease of BRS and RR
maskST_lag = maskS_lag.*maskT_lag;

% reconstruction for the introduced lag
maskSx = maskST_lag(brs_lag+1:end);
maskTx = maskST_lag(1:end-brs_lag);

indSx = indS+brs_lag;
indTx = indS;  %------- IS THIS RIGHT?????????

return
figure,clf
ha(1)=subplot(3,1,1); hold on
stem(1.2*maskT,'k')
stem(maskS,'r')
ha(2)=subplot(3,1,2); hold on
stem(1.4*maskST_lag,'m')
stem(1.2*maskT_lag,'k')
stem(maskS_lag,'r')
ha(3)=subplot(3,1,3); hold on
stem(1.2*maskTx,'k')
stem(maskSx,'r')
linkaxes(ha,'x')
zoom(gcf,'xon')
%keyboard

