%load voldists table
function Plot_Voldists_function(basepath, outpath)

load([outpath '\AttuneVolTable.mat'])
AttuneVolTable = AttuneVolTable(AttuneVolTable.QC_flag ==1, :);

  % Get cruise name from first filename
    cruisename =  split(AttuneVolTable.Filename{1}, '_'); 
    cruisename = cruisename{2}; 

figure
subplot(4,2,1)
s = pcolor(AttuneVolTable.StartDate, [1:57], AttuneVolTable.SynDist'./sum(AttuneVolTable.SynDist'));
caxis([0 .07])
set(s, 'EdgeColor', 'none')
yticks([1 18 26 42 50])
yticklabels(syn_volbins([1 18 26 42 50]))
ylabel('Cell Volume (\mum^{3})')
title([cruisename '- Syn Scaled'])
ax = gca;

subplot(4,2,3)
s = pcolor(AttuneVolTable.StartDate, [1:57], AttuneVolTable.SynDist');
yticks([1 18 26 42 50])
yticklabels(syn_volbins([1 18 26 42 50]))
set(s, 'EdgeColor', 'none')
ylabel('Cell Volume (\mum^{3})')

title([cruisename '- Syn'])

set(gca, 'XLim', ax.XLim)

subplot(4,2,2)
a = pcolor(AttuneVolTable.StartDate, [1:66], AttuneVolTable.EukDist'./sum(AttuneVolTable.EukDist'))
caxis([0 .07])
set(a, 'EdgeColor', 'none')
yticks([1 17 32 42 57])
yticklabels(euk_volbins([1 17 32 42 57]))
ylabel('Cell Volume (\mum^{3})')
title([cruisename '- Euks Scaled'])
set(gca, 'XLim', ax.XLim)

subplot(4,2,4)
a = pcolor(AttuneVolTable.StartDate, [1:66], AttuneVolTable.EukDist')
set(a, 'EdgeColor', 'none')
yticks([1 17 32 42 57])
yticklabels(euk_volbins([1 17 32 42 57]))
title([cruisename '- Euks'])
ylabel('Cell Volume (\mum^{3})')

set(gca, 'XLim', ax.XLim)

subplot(4,2,5)
%plot(AttuneVolTable.StartDate, AttuneVolTable.rad_sw, '.k')
if ismember('rad1_sw',AttuneVolTable.Properties.VariableNames)
    plot(AttuneVolTable.StartDate, AttuneVolTable.rad1_sw, '.k')
else
    plot(AttuneVolTable.StartDate, AttuneVolTable.rad_sw, '.k')
end
ylabel('Rad SW')
set(gca, 'XLim', ax.XLim)


subplot(4,2,6)
plot(AttuneVolTable.StartDate, AttuneVolTable.temperature, '.')
ylabel('Temperature')
yyaxis right
plot(AttuneVolTable.StartDate, AttuneVolTable.salinity, '.')
ylabel('Salinity')
set(gca, 'XLim', ax.XLim)


%%
%also sanity check on attune table

subplot(4,2,7)
plot(AttuneVolTable.StartDate, AttuneVolTable.Syn_count./AttuneVolTable.VolAnalyzed_ml, '.')
ylabel('Syn (cells ml^{-1})')
set(gca, 'XLim', ax.XLim)
title('Syn Conc')
subplot(4,2,8)
plot(AttuneVolTable.StartDate, AttuneVolTable.Euk_without_PE_leq10um_count./AttuneVolTable.VolAnalyzed_ml, '.')
ylabel('Euk (cells ml^{-1})')
title('Euk Conc')
set(gca, 'XLim', ax.XLim)


set(gcf, 'Position', [2233 -373.6667 774 918.6667])

print([basepath cruisename '_Summary'], '-djpeg')
