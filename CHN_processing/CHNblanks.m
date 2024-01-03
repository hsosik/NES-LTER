

CHNBLANKS = readtable('\\sosiknas1\Lab_data\LTER\CHN\data\CHN_BLANKS.xlsx');


%get rid of entries that don't have data. aka they don't have a date run
CHNBLANKS = CHNBLANKS(~isnan(CHNBLANKS.date_combusted) & ~isnan(CHNBLANKS.date_run),:);

%make BDL == 0
CHNBLANKS.umolC(CHNBLANKS.umolC=='BDL') = '0';
CHNBLANKS.umolN(CHNBLANKS.umolN=='BDL') = '0';

%combine date combusted and date run to a unique 16digit number that's
%unique to each blank that will be applied

%blanks that were not run get a value of 11111111 in date_run field
CHNBLANKS.unq_date_combustrun = str2num([num2str(CHNBLANKS.date_combusted) num2str(CHNBLANKS.date_run)]);

% blanks2use = table('size',[length(unq_chn_date) 7],'VariableNames',{'date_combusted','date_run','CHN_blank','date_run','nitrogen','carbon','unq_date_combustrun'});
blanks2use = [];

unq_chn_date = unique(CHNBLANKS.unq_date_combustrun);
placement = 1;
for count = 1:length(unq_chn_date)
    ind = find(CHNBLANKS.unq_date_combustrun == unq_chn_date(count));
    unqblank = categorical(unique(CHNBLANKS.CHN_blank(ind)));
    for count2 = 1:length(unqblank)
        ind2 = find(CHNBLANKS.unq_date_combustrun == unq_chn_date(count) & CHNBLANKS.CHN_blank == unqblank(count2));
        blanks2use(placement).date_combusted = CHNBLANKS.date_combusted(ind2(1));
        blanks2use(placement).date_run = CHNBLANKS.date_run(ind2(1));
        blanks2use(placement).CHN_blank = CHNBLANKS.CHN_blank(ind2(1));
        blanks2use(placement).date_run = CHNBLANKS.date_run(ind2(1));
        blanks2use(placement).nitrogen = nanmean(CHNBLANKS.umolN(ind2));
        blanks2use(placement).carbon = nanmean(CHNBLANKS.umolC(ind2));
        blanks2use(placement).unq_date_combustrun = unq_chn_date(count);
        blanks2use(placement).blank_estimated = CHNBLANKS.values_estimated(ind2(1));
        placement = placement +1;
    end
end
blanks2use = struct2table(blanks2use);

clear CHNBLANKS* ans count ind unq_chn_date

save \\sosiknas1\Lab_data\LTER\CHN\CHNblanks2use blanks2use