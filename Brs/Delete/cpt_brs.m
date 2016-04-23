function [brs] = cpt_brs(t,bp)


%% PARAMETERS FOR THE SEQUENCE METHOD

no1         = 3;                    % # successive points
dsys_thresh = 1;                    % in mmHg
drri_thresh = 0.005;                % in seconds
lag_inc     = 1;  


%% EXCTRACT SBP AND RR
fs_in = round(1/mean(diff(t)));     % current sampling frequenc

% sys and dia indices
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);
%rri     = diff([-Inf;is])/fs_in;
rri     = diff([-Inf;id])/fs_in;

     


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




%%  DECREASE in SBP and RR
incS_d = -diff([-Inf;bp(is)])>=dsys_thresh;
incT_d = -diff([-Inf;rri])>=drri_thresh;


% detect sequences with no1=3 or more successive increases in SYS 
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



% detect sequences with no1=3 or more successive increases in RRI 
strT_d                  = num2str(incT_d');
strT_d(isspace(strT_d)) = '';
iti_b_d                 = strfind(strT_d,['0' repmat('1',1,no1)]);
iti_e_d                 = strfind(strT_d,[repmat('1',1,no1) '0'])+no1-1;
if length(iti_e_d)<length(iti_b_d), iti_e_d(end+1) = length(incT_d); end

maskT_d = zeros(length(is),1);
for i=1:length(iti_b_d),
    maskT_d(iti_b_d(i):iti_e_d(i)) = 1;
end

%% Fitting 

% Correction for the delay between SYS and RR (increase)
maskS_inclag  = [zeros(lag_inc,1);maskS];
maskT_inclag  = [maskT;zeros(lag_inc,1)];

maskS_declag  = [zeros(lag_inc,1);maskS_d];
maskT_declag  = [maskT_d;zeros(lag_inc,1)];

% Simultaneous growth
maskST_inclag = maskS_inclag.*maskT_inclag;

maskST_declag = maskS_declag.*maskT_declag;

% Re-correction
maskSx = maskST_inclag(lag_inc+1:end);
maskTx = maskST_inclag(1:end-lag_inc);

maskSx_d = maskST_declag(lag_inc+1:end);
maskTx_d = maskST_declag(1:end-lag_inc);

d_SBP    = bp(is(find(maskSx)));
d_RR     = rri(find(maskTx));

d_SBP_d  = bp(is(find(maskSx_d)));
d_RR_d   = rri(find(maskTx_d));


xx = [d_RR; d_RR_d];
yy = [d_SBP; d_SBP_d];

% Calculating BRS
p  = polyfit(xx,yy,1);

brs = p(1);

