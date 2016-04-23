clear all; close all; clc


addpath('mfiles/')
fs1 = 50;

% read all files from young dir
% *************************************************************************
dir1 = 'data/mat_Y/';
d = dir([dir1 '*.mat']);
filenames = {d.name};

for i=1:numel(filenames),
    A = load(strcat(dir1,filenames{i}));
    
    fs_in   = round(1/mean(diff(A.t)));         % current sampling frequenc
    if fs_in~=fs1,
        A.bp = my_decimate(A.bp,fs_in/fs1,fs1);
        A.t  = my_decimate(A.t,fs_in/fs1,fs1);
    end
    
    % compute BRS
    [OUT.brs1(i),~,~] = cpt_brs(A.t,A.bp);
    clear A
end

% read all files from normotensive dir
% *************************************************************************
dir2 = 'data/mat_N/';
d = dir([dir2 '*.mat']);
filenames = {d.name};

for i=1:numel(filenames),
    A = load(strcat(dir2,filenames{i}));
    
    fs_in   = round(1/mean(diff(A.t)));         % current sampling frequenc
    if fs_in~=fs1,
        A.bp = my_decimate(A.bp,fs_in/fs1,fs1);
        A.t  = my_decimate(A.t,fs_in/fs1,fs1);
    end
    
    % compute BRS
    [OUT.brs2(i),~,~] = cpt_brs(A.t,A.bp);
    clear A
end


% read all files from normotensive dir
% *************************************************************************
dir3 = 'data/mat_H/';
d = dir([dir3 '*.mat']);
filenames = {d.name};

for i=1:numel(filenames),
    A = load(strcat(dir3,filenames{i}));
    
    fs_in   = round(1/mean(diff(A.t)));         % current sampling frequenc
    if fs_in~=fs1,
        A.bp = my_decimate(A.bp,fs_in/fs1,fs1);
        A.t  = my_decimate(A.t,fs_in/fs1,fs1);
    end
    
    % compute BRS
    [OUT.brs3(i),~,~] = cpt_brs(A.t,A.bp);
    clear A
end


OUT.brs = [OUT.brs1 OUT.brs2 OUT.brs3];

figure, clf, hold on
plot(OUT.brs1(1:2:end),OUT.brs1(2:2:end),'ro')
plot(OUT.brs2(1:2:end),OUT.brs2(2:2:end),'k*')
plot(OUT.brs3(1:2:end),OUT.brs3(2:2:end),'m+')
plot(-40:40,-40:40,'k-')
xlim([-10,15])
ylim([-10,15])
xlabel('trial 1')
ylabel('trial 2')
axis square

figure 
% Young
subplot(131), hold on
    plot(OUT.brs1(1:2:end),'ro')
    plot(OUT.brs1(2:2:end),'ko')
    plot((OUT.brs1(1:2:end)+OUT.brs1(2:2:end))/2,'ko-')
    %xlim([-100,100])
    ylim([-10,30])
    %xlabel('trial 1')
    %ylabel('trial 2')
    %axis square
subplot(132), hold on
% Normotensive
    plot(OUT.brs2(1:2:end),'ro')
    plot(OUT.brs2(2:2:end),'ko')
     plot((OUT.brs2(1:2:end)+OUT.brs2(2:2:end))/2,'ko-')
    %ylim([-10,100])
    ylim([-10,30])
    %xlabel('trial 1')
    %ylabel('trial 2')
    %axis square
subplot(133), hold on
% Normotensive
    plot(OUT.brs3(1:2:end),'ro')
    plot(OUT.brs3(2:2:end),'ko')
    plot((OUT.brs3(1:2:end)+OUT.brs3(2:2:end))/2,'ko-')
    %ylim([-10,100])
    ylim([-10,30])
    %xlabel('trial 1')
    %ylabel('trial 2')
    %axis square
