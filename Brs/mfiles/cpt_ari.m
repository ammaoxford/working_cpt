function [ari_i,err] = cpt_ari(bp,fv,fs)
%
% Computational Physiology Toolbox
%
% Adam Mahdi, University of Oxford
%
%--------------------------------------------
% INPUT
% bp,fv
%
% OUTPUT
% ari_i:      ARI index
% err:        vector of errors of the 10 Tiecks models
%--------------------------------------------

% ARI CONSTANTS
mBP = mean(bp);
mFV = mean(fv);
vec_const = [fs,mBP,mFV];


% Define ARI levels
ari_ref(1,:)  = [Inf    0      0];       % ARI=0
ari_ref(2,:)  = [2      1.6    0.20];    % ARI=1
ari_ref(3,:)  = [2      1.5    0.40];    % ARI=2
ari_ref(4,:)  = [2      1.15   0.60];    % ARI=3
ari_ref(5,:)  = [2      0.9    0.80];    % ARI=4
ari_ref(6,:)  = [1.9    0.75   0.90];    % ARI=5
ari_ref(7,:)  = [1.6    0.65   0.94];    % ARI=6
ari_ref(8,:)  = [1.2    0.55   0.96];    % ARI=7
ari_ref(9,:)  = [0.87   0.52   0.97];    % ARI=8
ari_ref(10,:) = [0.65   0.50   0.98];    % ARI=9

% Preallocation
vec_err = zeros(1,10);
t0 = 11;
for j=1:10,
    fv_sim = cpt_arisimul(ari_ref(j,:),bp,vec_const);
    vec_err(j) = norm(fv_sim(t0:end)-fv(t0:end),2)/mean(fv);   
    
    % % rescale simulated FV signal
    % SIGMA2  = nanmean(fv_sim(t0:end).^2);         % mean-square in each fv_sim
    % SIGMA12 = nanmean(fv_sim(t0:end).*fv(t0:end));
    % gain(j) = SIGMA12./SIGMA2;
    % fv_sim_tmp = fv_sim(t0:end)*gain(j);
    % nrmse(j) = (mean((fv(t0:end)-fv_sim_tmp).^2)/mean(fv(t0:end).^2)).^0.5;
end

% COMPUTE THE ARI index

% Integer ARI
% [err_i,ari_i] = min(vec_err);
% ari_i=ari_i-1;

% Interpolated ARI (cubic spline)
sx         = 0:9;
sxx        = 0:.1:9;
%sy         = nrmse;
sy         = vec_err;
yy         = spline(sx,sy,sxx);
[~,ari_vi] = min(yy);
ari_i      = sxx(ari_vi);
err        = yy(ari_vi);


%% ARI SIMULATION
function vsim = cpt_arisimul(parTDK,bp,vec_const)

%--------------------------------------------
% INPUT
% parTDK    = [T D K]    
% vec_const = [f mBP mFV];
%
% OUTPUT
% vsim - ARI response vector of the size bp
%--------------------------------------------

% Critical Closing Pressure
CCP = 12;

% PARAMETERS: ari_pars = [T D K]
T   = parTDK(1);
D   = parTDK(2);
K   = parTDK(3);

% CONSTANTS:  vec_const = [f mBP mFV];
f   = vec_const(1);
mBP = vec_const(2);
mFV = vec_const(3);

% NORMALIZED PRESSURE
dP  = (bp-mBP)./(mBP-CCP);

% INITIALIZE VECTOR
tsize = length(bp);
% x1 = zeros(1,tsize);
% x2 = zeros(1,tsize);

% INITIAL CONDITION
x1(1) = 2*D*dP(1);
x2(1) = dP(1);

% SIMULATE
for j=1:tsize-1,
    x2(j+1) = x2(j) + (x1(j)-2*D*x2(j))./(f*T);
    x1(j+1) = x1(j) + (dP(j)-x2(j))./(f*T);
end

% OUTPUT
vsim = mFV*(1+dP-K.*x2');
