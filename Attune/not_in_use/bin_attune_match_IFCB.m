load('\\sosiknas1\Lab_data\Attune\EN608\Summary\Attune')
IFCB = load('\\sosiknas1\IFCB_products\NESLTER_transect\summary\21Aug2018\IFCB_biovolume_size_classes_manual_21Aug2018');

maxdf = 20/60/24; %window size in time, 20 minutes as days
%Attune_match_IFCB = struct('Count', NaN(size(IFCB.matdate)), 'Biovol', NaN(size(IFCB.matdate)))
%f1 = fieldnames(Attune_match_IFCB);
%initialize
temp = NaN(size(IFCB.matdate));
f1 = fieldnames(Attune.Count);
for indf = 1:length(f1)
    Attune_match_IFCB.Count.(f1{indf}) = temp;
end
f2 = fieldnames(Attune.Biovol);
for indf = 1:length(f2)
    Attune_match_IFCB.Biovol.(f2{indf}) = temp;
end
Attune_match_IFCB.vol_analyzed = temp;
Attune_match_IFCB.matdate = temp;

for ifile = 1:length(IFCB.filelist)
    absdf = abs(Attune.FCSfileinfo.matdate_start-IFCB.matdate(ifile));
    ipts = find(absdf <= maxdf);
    if ~isempty(ipts)    
        for indf = 1:length(f1)
            Attune_match_IFCB.Count.(f1{indf})(ifile) = nansum(Attune.Count.(f1{indf})(ipts));
        end
        for indf = 1:length(f2)
            Attune_match_IFCB.Biovol.(f2{indf})(ifile) =nansum(Attune.Biovol.(f2{indf})(ipts));
        end
        Attune_match_IFCB.vol_analyzed(ifile) = nansum(Attune.FCSfileinfo.vol_analyzed(ipts))/1e6;
        Attune_match_IFCB.matdate(ifile) = nanmean(Attune.FCSfileinfo.matdate_start(ipts));
    end
end


figure
plot(IFCB.N10_20_phyto./IFCB.ml_analyzed, Attune_match_IFCB.Count.tentwen./Attune_match_IFCB.vol_analyzed, '.')
line(xlim, xlim)
xlabel('IFCB cells ml^{-1}, 10-20 \mum')
ylabel('Attune cells ml^{-1}, 10-20 \mum')

figure
plot(IFCB.matdate, IFCB.N10_20_phyto./IFCB.ml_analyzed)
hold on
plot(Attune_match_IFCB.matdate, Attune_match_IFCB.Count.tentwen./Attune_match_IFCB.vol_analyzed,'.-')
datetick
legend('IFCB', 'Attune')
ylabel('Cells ml^{-1}, 10-20 \mum')
