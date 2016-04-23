function [brs,vv,brsdata] = cpt_brs(t,bp)

% This function calculates BRS using the sequence method
% If RR are not evailable they are computed from DBP

% Dependencies:
% [maskSx,maskTx,maskS,maskT] = cpt_brsSEQ_incdec(bp,is,id,rri, params,opt)
% [is,id]                     = cpt_abpTimes4x(bp,fs_in,sep);


%% EXCTRACT SBP AND RR
fs_in   = round(1/mean(diff(t)));         % current sampling frequency
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);   % sys and dia indices
rri     = diff([-Inf;id])/fs_in;          % cardiac cycles length using DBP
%rri     = diff([-Inf;is])/fs_in;          % cardiac cycles length using SBP


%% PARAMETERS FOR THE SEQUENCE METHOD

no1         = 2;                    % # successive points
dsys_thresh = 1;                    % in mmHg
drri_thresh = 0.005;                % in seconds
brs_lag     = 1;                    % correction for delay RR (increase)
params      = [fs_in,no1,dsys_thresh,drri_thresh,brs_lag];


%% MAIN: INCREASE & DECREASE SEQUENCES
% maskS:  indices of BRS sequences
% maskSx: lagged & simult. indices of BRS sequences
% etc. '_inc'- increase; '_dec' - for decrease

[maskSx_inc,maskTx_inc,maskS_inc,maskT_inc] = cpt_brsSEQ_incdec(bp,rri,is,id,params,1);
[maskSx_dec,maskTx_dec,maskS_dec,maskT_dec] = cpt_brsSEQ_incdec(bp,rri,is,id,params,2);



%% BRS COMPUTATIONS 

dSBP_inc   = bp(is(find(maskSx_inc)));
dRR_inc    = rri(find(maskTx_inc));

dSBP_dec   = bp(is(find(maskSx_dec)));
dRR_dec    = rri(find(maskTx_dec));

% figure, hold on
% set(gcf,'Units','normalized','Position',[0.2 0.2 0.7 0.5])
% plot(bp,'k')
% plot(is,bp(is),'ko')
% plot(id,bp(id),'ko')
% % plot(is(find(maskSx_inc)),bp(is(find(maskSx_inc))),'ro','MarkerFaceColor','r')
% % plot(is(find(maskSx_dec)),bp(is(find(maskSx_dec))),'bo','MarkerFaceColor','b')
% plot(is(find(maskS_inc)),bp(is(find(maskS_inc))),'ro','MarkerFaceColor','r')
% plot(is(find(maskS_dec)),bp(is(find(maskS_dec))),'bo','MarkerFaceColor','b')
% keyboard

% union of increase and decrease sequences
yy  = [ dRR_inc;  dRR_dec]*1000;
xx  = [dSBP_inc; dSBP_dec];

%xx = xx-mean(xx);
%yy = yy-mean(yy);

% Local Mean Detrending

% fit the line (least squares sense)
% p   = polyfit(dRR_inc,dSBP_inc,1);
p   = polyfit(xx,yy,1);
brs = p(1);             % BRS (slope of the line)
vv  = polyval(p,xx);    % both coefficients

brsdata = [xx yy];
