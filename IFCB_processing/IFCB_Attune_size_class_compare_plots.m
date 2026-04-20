%clist = {'EN608' 'EN617' 'EN627' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661' 'EN668'}; 
%clist = {'AR28' 'AR31' 'AR34' 'AR39' 'AR61B' 'AR62' 'AR63' 'AT46'};

%clist = {'EN608' 'EN617' 'EN627' 'EN644' 'EN649' 'EN655' 'EN657' 'EN661' 'EN668' 'AR28' 'AR31' 'AR34' 'AR39' 'AR61B' 'AR62' 'AR63'};% 'AT46'}; 

p = '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\summary_files\';
clist = dir([p '*_IFCB_Attune_merge_sum.mat']);
clist = {clist.name}';
%clist = {'TN368' 'RB1904' 'AR29'}; %p = '\\sosiknas1\IFCB_products\SPIROPA\summary\';
%clist = {'EN695' 'HRS2303' 'EN706'};
clist = {'AR99'};
%clist = {'20220713_EN685_IFCB_Attune_merge_sum.mat'};
outdir =  '\\sosiknas1\Lab_data\Attune\cruise_data\IFCB_Attune_merge\plots_Scatter_TimeSeries\';
plotEachSample = 1;
attunebase = '\\sosiknas1\Lab_data\Attune\cruise_data\';
for count = 1:length(clist)
    %cruiseStr = clist{count}; 
    f = dir([p '*' clist{count} '*']);
    fout = regexprep(f.name, '_IFCB_Attune_merge_sum.mat', '');
    if ~exist(outdir, 'dir')
        mkdir(outdir)
    end

    load([ p char(f.name)])
    binrange = [binedges(1:end-1); binedges(2:end); 1:length(binedges)-1];
    %%
    figure
    tl = tiledlayout(2,3);
    for c = 1:5
        nexttile
        b = c+2;
        %plot([ifcbHcount.diatom(:,b)+ifcbHcount.dino(:,b)+ifcbHcount.nano(:,b)+ifcbHcount.ciliate(:,b)+ifcbHcount.otherPhyto(:,b)]./ifcbHmeta.ml_analyzed, attuneHcount.Euk(:,b)./attuneml', '.')
        plot(ifcbHcount.protist_tricho(:,b)./ifcbHmeta.ml_analyzed, attuneHcount.Euk(:,b)./attuneml', '.')
        line(xlim, xlim)
        title([num2str(binrange(1,b)) '-' num2str(binrange(2,b)) ' \mum'])
    end
    ylabel(tl, 'Attune, ml^{-1}')
    xlabel(tl, 'IFCB, ml^{-1}')
    title(tl, regexprep(f.name, '_IFCB_Attune_merge_sum.mat', ''),'interpreter', 'none')
    print([outdir fout 'SizeRangeScatter' datestr(now,'ddmmmyyyy') '.png'], '-dpng')

    %%
    figure
    tl = tiledlayout(5,1);
    dt = datetime(ifcbHmeta.sample_time, 'InputFormat', 'yyyy-MM-dd hh:mm:ss+00:00');
    ifcbHmeta.dt = datetime(ifcbHmeta.sample_time, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
    ifcbHmeta.mdate = datenum(ifcbHmeta.dt);
    [~,ind] = sort(dt);

    for c = 1:5
        nexttile
        b = c+2;
%        plot(dt(ind), [ifcbHcount.diatom(ind,b)+ifcbHcount.dino(ind,b)+ifcbHcount.nano(ind,b)+ifcbHcount.ciliate(ind,b)]./ifcbHmeta.ml_analyzed(ind), '-')
        plot(dt(ind), ifcbHcount.protist_tricho(ind,b)./ifcbHmeta.ml_analyzed(ind), '-')
        hold on
        plot(dt(ind), attuneHcount.Euk(ind,b)./attuneml(ind)', '-')
        title([num2str(binrange(1,b)) '-' num2str(binrange(2,b)) ' \mum'])
    end
    ylabel(tl, 'Concentration, ml^{-1}')
    %xlabel(tl, 'IFCB, ml^{-1}')
    title(tl, regexprep(f.name, '_IFCB_Attune_merge_sum.mat', ''),'interpreter', 'none')
    orient tall
    legend(tl.Children(end), 'IFCB', 'Attune', 'Location', 'northeast')
    print([outdir fout 'SizeRangeTimeSeries' datestr(now,'ddmmmyyyy') '.png'], '-dpng')

%%
   % ifcbHcount_sumConc = (ifcbHcount.diatom+ifcbHcount.dino+ifcbHcount.nano+ifcbHcount.ciliate+ifcbHcount.otherPhyto)./ifcbHmeta.ml_analyzed;
   % ifcbHcount_sumConc2 = ifcbHcount.particle./ifcbHmeta.ml_analyzed;
    ifcbHcount_sumConc = ifcbHcount.protist_tricho./ifcbHmeta.ml_analyzed;
    ifcbHcount_sumConc2 = (ifcbHcount.protist_tricho+ifcbHcount.Detritus)./ifcbHmeta.ml_analyzed; 
    if exist('attunemlSyn', 'var') %case for alternating runs (e.g., EN688)
        attuneHcount_sumConc = attuneHcount.Euk./attunemlEuk' + attuneHcount.Syn./attunemlSyn';
        attuneml = attunemlEuk; %used only for value displayed in plot legend
    else
        attuneHcount_sumConc = (attuneHcount.Euk+attuneHcount.Syn)./attuneml';
    end
    diamCenters = mean([binedges(1:end-1); binedges(2:end)]);
    if plotEachSample
        attunebase2 = [attunebase fout '\bead_calibrated\class\'];
        attunebase_out = [attunebase2 'sizeDistPlots\'];
        if ~exist(attunebase_out, 'dir')
            mkdir(attunebase_out)
        end
        for countSample = 1:size(ifcbHmeta,1)
            figure(99), clf
            semilogy(diamCenters,ifcbHcount_sumConc(countSample,:), 'b.-')
            hold on 
            semilogy(diamCenters,ifcbHcount_sumConc2(countSample,:), 'b-.')
            semilogy(diamCenters,attuneHcount_sumConc(countSample,:), 'g.-')
            title(ifcbHmeta.pid(countSample), 'interpreter', 'none')
            legend(['IFCB ' num2str(round(ifcbHmeta.ml_analyzed(countSample),2)) ' ml'], 'IFCB w/detritus', ['Attune +/-12min ' num2str(attuneml(countSample)) ' ml'])
            ylabel('Cells ml^{-1} per bin')
            xlabel('ESD (\mum)')
            grid
            set(gca, 'xtick',0:5:50, 'xlim',[0 50])
        %    pause
            print([attunebase_out ifcbHmeta.pid{countSample} '.png'], '-dpng')
        end
    end
    clear attuneml* 
end