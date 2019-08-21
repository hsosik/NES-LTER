%Comparison script to see differences in model runs from
% old vs. new processing
% two component vs. one component


%% Old vs. New Processing:


for year2do=2010 %2003:2018
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    if ismember(year2do,2016:2018) %do not have July 2016 data for most of these years
        eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/output_Jan2019/mvco_14par_dmn_' num2str(year2do) '.mat'])
        modelres_old = modelresults;
    else
        eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/output_July2016/mvco_14par_dmn_' num2str(year2do) '.mat'])
        modelres_old = modelresults;
    end
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/output_June2019/mvco_14par_dmn_' num2str(year2do) '.mat'])
    modelres_new = modelresults;
    
    [~,io,in] = intersect(modelres_old(:,1),modelres_new(:,1)); %days may have changed with new processing
    
    subplot(4,4,year2do-2002)
    line([0 2],[0 2],'color',[0.5 0.5 0.5],'linewidth',2), hold on
    plot(modelres_old(io,17),modelres_new(in,17),'.')
    xlabel('Old processing')
    ylabel('New processing')   
    title(num2str(year2do))
    axis([0 2 0 2])
end



%% One component vs. Two Oooooo....


for year2do=2003:2018
    
    switch year2do
        case 2003
            filelabel='May';
        case 2004
            filelabel='Apr';
        case 2005
            filelabel='Apr';
        case 2006
            filelabel='May';
        case 2007
            filelabel='Mar';
        otherwise
            filelabel='Jan';
    end
    
    eval(['load /Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/output_June2019/mvco_14par_dmn_' num2str(year2do) '.mat'])
    modelresults_two=modelresults;
    
    filename=['/Volumes/Lab_data/MVCO/FCB/MVCO_' filelabel num2str(year2do) '/model/onecomp_output_June2019/mvco_7par_dmn_' num2str(year2do) '.mat'];
    load(filename)
    if exist(['modelresults_one' num2str(year2do)],'var')
        eval(['modelresults_one=modelresults_one' num2str(year2do) ';'])
        disp(num2str(year2do))
    end
    
    if size(modelresults_two,1)~=size(modelresults_one,1)
        disp('Uh-oh...sizes of modelresults are not the same?')
    end
    
    subplot(4,4,year2do-2002)
    plot(modelresults(:,17),modelresults_one(:,10),'.')
    xlabel('Two components')
    ylabel('One component')
    line([0 2],[0 2],'color',[0.4 0.4 0.4])
    title(num2str(year2do))
    
    clearvars -except year2do
end

% HMmmmmm...would seem that two components are important! Woot!


