%% NES Strata
% Isabel Honda
% Strata coordinates provided by Harvey Walsh from the NOAA Northeast
% Fisheries Science Center

%% Initialization
clc; clear; close all
strata47 = load('/Users/isabelhonda/Dropbox (MIT)/Research/Data/NES_polygons47.mat');
strata46 = load('/Users/isabelhonda/Dropbox (MIT)/Research/Data/NES_polygons46.mat');
strata26 = load('/Users/isabelhonda/Dropbox (MIT)/Research/Data/NES_polygons26.mat');

%% Plotting 47 Strata

figure(1)
clf
set(gcf,'color','w');
xlim([-78 -64])
ylim([35 45])
clear h
hold on
for i=1:47
    if i<=13 %MAB
        L(1) = plot(strata47.shape(i),'FaceColor','red');
    elseif i>13 && i<=25 %SNE
        L(2) = plot(strata47.shape(i),'FaceColor','blue');
    elseif i>25 && i<=32 %GB
        L(3) = plot(strata47.shape(i),'FaceColor','green');
    else
        L(4) = plot(strata47.shape(i),'FaceColor','magenta'); %cyan or yellow instead?
        
    end    
    [xn,yn] = centroid(strata47.shape(i));
    if i==6 || i==17 || i==13 || i==36 || i==46
        xn=xn-0.15;
    elseif i==12 || i==21 || i==40 || i==41
        xn=xn-0.25;
    elseif i== 20 || i==31 || i==32
        yn=yn+0.15;
    elseif i==22 || i==24 
        yn=yn-0.1;
    elseif i==33 || i==35
        yn=yn-0.25;
    elseif i==38
        xn=xn-0.9;
    elseif i==39
        xn=xn-0.2;
    elseif i==45
        xn=xn-0.5;
        yn=yn-0.1;
    end
    
    text(xn,yn,num2str(i));
    numPos_old(i,:) = [xn yn];
end
hleg = legend(L, {'MAB','SNE','GB','GoM'},'Location','Southeast')
htitle = get(hleg,'Title');
set(htitle,'String','Regions')
xlabel('Longitude (º)')
ylabel('Latitude (º)')
axis square


%% Plotting 46 Strata

figure(1)
clf
set(gcf,'color','w');
xlim([-78 -64])
ylim([35 45])
clear h
hold on
for i=1:46
    if i<=13 %MAB
        L(1) = plot(strata46.shape(i),'FaceColor','red');
    elseif i>13 && i<=25 %SNE
        L(2) = plot(strata46.shape(i),'FaceColor','blue');
    elseif i>25 && i<=32 %GB
        L(3) = plot(strata46.shape(i),'FaceColor','green');
    else
        L(4) = plot(strata46.shape(i),'FaceColor','magenta'); %cyan or yellow instead?
        
    end    
    [xn,yn] = centroid(strata46.shape(i));
    if i==6 || i==17 || i==13 || i==36 || i==46
        xn=xn-0.15;
    elseif i==12 || i==21 || i==40 || i==41
        xn=xn-0.25;
    elseif i== 20 || i==31 || i==32
        yn=yn+0.15;
    elseif i==22 || i==24 
        yn=yn-0.1;
    elseif i==33 || i==35
        yn=yn-0.25;
    elseif i==38
        xn=xn-0.9;
    elseif i==39
        xn=xn-0.2;
    elseif i==45
        xn=xn-0.5;
        yn=yn-0.1;
    end
    
    text(xn,yn,num2str(i));
    numPos_old(i,:) = [xn yn];
end
hleg = legend(L, {'MAB','SNE','GB','GoM'},'Location','Southeast')
htitle = get(hleg,'Title');
set(htitle,'String','Regions')
xlabel('Longitude (º)')
ylabel('Latitude (º)')
axis square



%% Plotting 26 Strata

figure(3)
set(gcf,'color','w');
clf
hold on
for i=1:26
        if i<=6
            h(1) = plot(strata26.shape(i),'FaceColor','red')
        elseif i>6 && i<=12
            h(2) = plot(strata26.shape(i),'FaceColor','blue')
        elseif i>12 && i<=16
            h(3) = plot(strata26.shape(i),'FaceColor','green')
        else
            h(4) = plot(strata26.shape(i),'FaceColor','magenta')      
        end
        [xn,yn] = centroid(strata26.shape(i));
    text(xn,yn,num2str(i));
end
xlabel('Longitude (º)')
ylabel('Latitude (º)')
hleg = legend(h,'MAB','SNE','GB','GoM','location','se');
htitle = get(hleg,'Title');
set(htitle,'String','Regions')
axis square


