%% Copepod Variance Ratio Calculations
% Isabel Honda
% January/February 2023

%% Copepod Community Synchrony - GAM Generated Timeseries

clc; clear; close all
cd('/Users/isabelhonda/Dropbox (MIT)/Research/VR/')

locs = ["MVCO","L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11"];
comm_VRs_entire = zeros(length(locs),2);
comm_VRs_int = zeros(length(locs),2);
comm_VRs_seas = zeros(length(locs),2);

figure(1)
clf
set(gcf,'color','w');

for i=1:length(locs)
    calfin = readtable(sprintf('GAM_transectStations/calfin_%s_cubedRootNm3.csv',locs(i)));
    pseudo = readtable(sprintf('GAM_transectStations/pseudo_%s_cubedRootNm3.csv',locs(i)));
    ctyp = readtable(sprintf('GAM_transectStations/ctyp_%s_cubedRootNm3.csv',locs(i)));
    comm_loc_entire = [calfin.GAMest pseudo.GAMest ctyp.GAMest]; % 
    comm_VRs_entire(i) = variance_ratio(comm_loc_entire);
    [lb,ub] = VR_sigTest(comm_loc_entire);
    if comm_VRs_entire(i,1) > ub || comm_VRs_entire(i,1) < lb
        comm_VRs_entire(i,2) = 1;
    end

    avgYear_calfin = accumarray(calfin.year,calfin.GAMest,[],@mean);
    avgYear_pseudo = accumarray(pseudo.year,pseudo.GAMest,[],@mean);
    avgYear_ctyp = accumarray(ctyp.year,ctyp.GAMest,[],@mean);
    comm_loc_int = [nonzeros(avgYear_calfin) nonzeros(avgYear_pseudo) nonzeros(avgYear_ctyp)]; %  
    comm_VRs_int(i) = variance_ratio(comm_loc_int);
    [lb,ub] = VR_sigTest(comm_loc_int);
    if comm_VRs_int(i,1) > ub || comm_VRs_int(i,1) < lb
        comm_VRs_int(i,2) = 1;
    end

    seas_calfin = accumarray(calfin.julian,calfin.GAMest,[],@mean);
    seas_pseudo = accumarray(pseudo.julian,pseudo.GAMest,[],@mean);
    seas_ctyp = accumarray(ctyp.julian,ctyp.GAMest,[],@mean);
    comm_loc_seas = [seas_calfin seas_pseudo seas_ctyp]; %  
    comm_VRs_seas(i) = variance_ratio(comm_loc_seas);
    [lb,ub] = VR_sigTest(comm_loc_seas);
    if comm_VRs_seas(i,1) > ub || comm_VRs_seas(i,1) < lb
        comm_VRs_seas(i,2) = 1;
    end

    if i==1
        subplot(3,3,1)
        plot(calfin.year + calfin.julian/365,comm_loc_entire,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Entire Timeseries at MVCO')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Year')

        subplot(3,3,4)
        plot(1977:2019,comm_loc_int,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Interannual Variability at MVCO')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Year')
        
        subplot(3,3,7)
        plot(1:365,comm_loc_seas,'LineWidth',2)
        xlim([1 365])
        legend('Calfin','Pseudo','Ctyp')
        title('Seasonal Variability at MVCO')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Day of Year')

    elseif i==12
        subplot(3,3,2)
        plot(calfin.year + calfin.julian/365,comm_loc_entire,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Entire Timeseries at L11')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Year')

        subplot(3,3,5)
        plot(1977:2019,comm_loc_int,'LineWidth',2)
        xlim([1977 2019])
        legend('Calfin','Pseudo','Ctyp')
        title('Interannual Variability at L11')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Year')
        
        subplot(3,3,8)
        plot(1:365,comm_loc_seas,'LineWidth',2)
        xlim([1 365])
        legend('Calfin','Pseudo','Ctyp')
        title('Seasonal Variability at L11')
        ylabel('N^{1/3} m^{-1}')
        xlabel('Day of Year')
    end
end


subplot(3,3,3)
plot(1:length(locs),comm_VRs_entire(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_entire(i,2) == 1
        h(1)=plot(i,comm_VRs_entire(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Entire Timeseries VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')

subplot(3,3,6)
plot(1:length(locs),comm_VRs_int(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_int(i,2) == 1
        h(1)=plot(i,comm_VRs_int(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Interannual Variability VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')

subplot(3,3,9)
plot(1:length(locs),comm_VRs_seas(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(locs)
    if comm_VRs_seas(i,2) == 1
        h(1)=plot(i,comm_VRs_seas(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
xlabel('Station')
ylabel('Variance Ratio')
title('Seasonal Variability VR')
xlim([1 length(locs)])
legend(h,'Statistically Significant (95% CI)','Location','southeast')


%% EcoMon Rectangles around Transect stations - Calfin, Pseudo, Ctyp

clc; clear; close all

locs = ["MVCO","L1","L2","L3","L4","L5","L6","L7","L8","L9","L10","L11"];
stations = readtable("~/Dropbox (MIT)/Research/Data/transectLocs.csv");
relStations = stations([26,1:11],:);
oldStrats = load('/Users/isabelhonda/Dropbox (MIT)/Research/Data/NES_polygons46.mat');

figure(2)
clf
set(gcf,'color','w');
subplot(2,3,3)
for i=18:25
    plot(oldStrats.shape(i),'EdgeColor','black','FaceColor','white');  
    hold on
    [xn,yn] = centroid(oldStrats.shape(i));
    if i==24
        xn = xn+0.2;
    end
    %text(xn,yn-0.1,num2str(i));
end
xlabel('Longitude')
ylabel('Latitude')
title('Stations')
axis square

VRs = zeros(length(locs),2);
for i=1:length(locs)
    subplot(2,3,3)
    v = rand(3,1); 
    curPolyVert = readtable(sprintf('EcoMon_transectRectangle/%s_coordsRectangle_0.5deg.csv',locs(i)));
    curPoly = polyshape(curPolyVert.V2,curPolyVert.V1);  
    plot(curPoly,'FaceColor',v)
    plot(relStations.longitude(i),relStations.latitude(i),'*','Color',v)
    text(relStations.longitude(i) + 0.1,relStations.latitude(i),relStations.name(i));

    curLocDat = readtable(sprintf('EcoMon_transectRectangle/%s_dataRectangle_0.5deg.csv',locs(i)));
    VRs(i) = variance_ratio([curLocDat.calfin_100m3 curLocDat.pseudo_100m3 curLocDat.ctyp_100m3]);
    [lb,ub] = VR_sigTest([curLocDat.calfin_100m3 curLocDat.pseudo_100m3, curLocDat.ctyp_100m3]);
    if VRs(i,1) > ub || VRs(i,1) < lb
        VRs(i,2) = 1;
    end

    if i==1
        subplot(2,3,[1 2])
        clear h
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.calfin_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq)
        hold on
        h(1) = plot(curLocDat.year + curLocDat.month/12, curLocDat.calfin_100m3,'.','color',"#0072BD")
        
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.pseudo_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq,'color',"#D95319")
        h(2) = plot(curLocDat.year + curLocDat.month/12, curLocDat.pseudo_100m3,'.','color',"#D95319")
        
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.ctyp_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq,'color',"#EDB120")
        h(3) = plot(curLocDat.year + curLocDat.month/12, curLocDat.ctyp_100m3,'.','color',"#EDB120")
    
        xlim([1977 2020])
        title({locs(i),sprintf('VR = %f',VRs(i,1))})
        legend(h,'Calfin','Pseudo','Ctyp')
        ylabel('N^{1/3} m^{-1}')

    elseif i==12
        subplot(2,3,[4 5])
        clear h
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.calfin_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq)
        hold on
        h(1) = plot(curLocDat.year + curLocDat.month/12, curLocDat.calfin_100m3,'.','color',"#0072BD")
        
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.pseudo_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq,'color',"#D95319")
        h(2) = plot(curLocDat.year + curLocDat.month/12, curLocDat.pseudo_100m3,'.','color',"#D95319")
        
        vq = interp1(curLocDat.year + curLocDat.month/12,curLocDat.ctyp_100m3,1977:1/12:2020,'linear');
        plot(1977:1/12:2020, vq,'color',"#EDB120")
        h(3) = plot(curLocDat.year + curLocDat.month/12, curLocDat.ctyp_100m3,'.','color',"#EDB120")
    
        xlim([1977 2020])
        title({locs(i),sprintf('VR = %f',VRs(i,1))})
        legend(h,'Calfin','Pseudo','Ctyp')
        ylabel('N^{1/3} m^{-1}')

    end

end

subplot(2,3,6)
plot(1:length(locs),VRs(:,1),'LineWidth',2,'Color','#a87d60')
hold on
for i=1:length(VRs)
    if VRs(i,2) == 1
        h(1)=plot(i,VRs(i,1),'*','Color','#a87d60','Linewidth',2)
    end   
end
set(gca, 'Xtick',1:length(locs),'XTickLabel',locs);
ylabel('Variance Ratio')
title('Variance Ratios')
%xlim([1 length(titles)])
legend(h,'Statistically Significant (95% CI)','Location','NW')


%%

function VR = variance_ratio(yy)
    % yy should be in the form of a matrix with different timeseries in each column
    % For example, columns of yy may contain the same species at different locations, or
    % different species at the same locations

    varsum = sum(nanvar(yy),'omitnan');
    covsum = sum(triu(nancov(yy),1),'all');
    VR=(varsum+2*covsum)/varsum;
end


% Significance test function
function [lb_95,ub_95] = VR_sigTest(yy)
    % Adapted from Ji's "Ji_idealized.m" reshuffling method
    % May result in spurious compensatory/synchrony significance results as in Solow & Dupisea (2007) 
    
    %VR = variance_ratio(yy);

    [m,n]=size(yy);
    vrs_all=zeros(1000,1);
    yys=zeros(m,n);

    for i=1:1000     %permute 1000 times
        for nn=1:n           
            %random permute
            mp=randperm(m);
            yys(:,nn)=yy(mp,nn);            
        end
        
        % Calculating VR for randomly reshuffled yy timeseries       
        vrs_all(i)=variance_ratio(yys);    
    end

    mu=nanmean(vrs_all);
    ss=nanstd(vrs_all);

    % for 95% confidence bound    
    lb_95=mu-1.96*ss;   %lower bound
    ub_95=mu+1.96*ss;   %uppder bound
    
    % for 90% confidence bound    
    lb_90=mu-1.64*ss;   %lower bound
    ub_90=mu+1.64*ss;   %uppder bound
    
    % VR must be outside this limit to be significant    
end


