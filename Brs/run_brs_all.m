clear all;close all; clc 


%% ADD PATHS & READ DATA

addpath('mfiles/')
addpath('data/mat_Y')
datanames_Y;
load(strcat('mat_Y/',dat_Y{8}))

% t  = t(1:3000);
% bp = bp(1:3000);


%% EXCTRACT SBP AND RR
fs_in = round(1/mean(diff(t)));     % current sampling frequenc

% sys and dia indices
[is,id] = cpt_abpTimes4x(bp,fs_in,0.9);
%rri     = diff([-Inf;id])/fs_in; 
rri     = diff([-Inf;is])/fs_in;


%% PARAMETERS FOR THE SEQUENCE METHOD

no1         = 3;                    % # successive points
dsys_thresh = 1;                    % in mmHg
drri_thresh = 0.005;                % in seconds
lag_inc     = 1;                    % correction for delay RR (increase)
params      = [fs_in,no1,dsys_thresh,drri_thresh,lag_inc];


%% MAIN COMPUTATION
% [maskSx_i,maskTx_i,maskS_i,maskT_i] = cpt_brsSEQ_inc(bp,rri,is,id,params);
% [maskSx_d,maskTx_d,maskS_d,maskT_d] = cpt_brsSEQ_dec(bp,rri,is,id,params);

[maskSx_i,maskTx_i,maskS_i,maskT_i,indSx_i,indTx_i,indS_i,indT_i] = cpt_brsSEQ_incdec(bp,rri,is,id,params,1);
[maskSx_d,maskTx_d,maskS_d,maskT_d,indSx_d,indTx_d,indS_d,indT_d] = cpt_brsSEQ_incdec(bp,rri,is,id,params,2);


figure
hold on
set(gcf,'Units','normalized','Position',[0.2 0.2 0.7 0.6])
ha(1)=subplot(211); hold on
plot(t,bp,'k')
plot(t(id),bp(id),'ko')
plot(t(is(find(maskS_i))),bp(is(find(maskS_i))),'ro','MarkerFaceColor','r')
plot(t(is(find(maskS_d))),bp(is(find(maskS_d))),'bo','MarkerFaceColor','b')
ylabel('BP [mmHg]')
ha(2)=subplot(212); hold on
tr = (id(1)-1)/fs_in+[0; cumsum(rri(2:end))];
plot(tr,rri,'k')
plot(tr(find(maskT_i)),rri(find(maskT_i)),'ro','MarkerFaceColor','r')
plot(tr(find(maskT_d)),rri(find(maskT_d)),'bo','MarkerFaceColor','b')
ylabel('RRI [s]')
xlabel('time [s]')
linkaxes(ha,'x')
zoom(gcf,'xon')

%return
     
     

%% Fitting 

% method 1
dSBP_i   = bp(is(find(maskSx_i)));
dRR_i    = rri(find(maskTx_i));

dSBP_d   = bp(is(find(maskSx_d)));
dRR_d    = rri(find(maskTx_d));

xx = [ dRR_i;  dRR_d];
yy = [dSBP_i; dSBP_d];
% xx = [ dRR_i];
% yy = [dSBP_i];

p  = polyfit(xx,yy,1);
vv = polyval(p,xx);


figure, hold on
set(gcf,'Units','normalized','Position',[0.3 0.2 0.4 0.6])
plot(dRR_i,dSBP_i,'ro')
plot(dRR_d,dSBP_d,'bo')
% plot(xx,yy,'ro')
plot(xx,vv)

p(1)

% method 2
for i=1:size(indSx_i,1),
    pp_i(i,:) = polyfit(rri(indSx_i(i,1):indSx_i(i,2)),bp(is(indSx_i(i,1):indSx_i(i,2))),1);
end
for i=1:size(indSx_d,1),
    pp_d(i,:) = polyfit(rri(indSx_d(i,1):indSx_d(i,2)),bp(is(indSx_d(i,1):indSx_d(i,2))),1);
end

[mean(pp_i(:,1)) mean(pp_d(:,1)) mean([pp_i(:,1); pp_d(:,1)])]

%return


%% PLOTTING
% [isi_b(:) isi_e(:)]
% [iti_b(:) iti_e(:)]

% figure, hold on 
% % stem(1.4*maskS.*maskT,'m')
% stem(2*maskST_inclag)
% stem(1.2*maskT,'r')
% stem(maskS)


figure, hold on,
set(gcf,'Units','normalized','Position',[0.2 0.2 0.7 0.6])
subplot(211), hold on
plot(t,bp,'k')
plot(t(is(find(maskS_i))),bp(is(find(maskS_i))),'ro','MarkerFaceColor','r')
plot(t(is(find(maskS_d))),bp(is(find(maskS_d))),'bo','MarkerFaceColor','b')
ylabel('BP')
subplot(212), hold on
plot(rri,'k')
plot(find(maskT_i),rri(find(maskT_i)),'ro','MarkerFaceColor','r')
plot(find(maskT_d),rri(find(maskT_d)),'bo','MarkerFaceColor','b')
ylabel('RRI')

    
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