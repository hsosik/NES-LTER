%% Starting from Vol Table for a single cruise
%let's make a folder with input .mat files
% for each day we want to apply our division rate model 

%this script requires manual checks to see if light data looks ok and dawn
%is correctly identified. 

function Setup_Cruise_Model_Days_Syn(AttuneVolTable)

%save them all here 
outpath = '\\sosiknas1\Lab_data\Attune\cruise_data\Division-rate-model\';

% Get cruise name from first filename
cruisename =  split(AttuneVolTable.Filename{1}, '_'); 
cruisename = cruisename{2}; 

    %first make a directory for the cruise
    mkdir([outpath cruisename])
    
    %cut down data to just this cruise
    %subset = TallT(TallT.cruisenames == cruiselist(c), :); 
    subset = sortrows(AttuneVolTable, 2); %make sure in chronological order
    
    %find start of days
    nighttime = find(subset.rad_sw <= 10);
    dawnind = nighttime((nighttime(2:end) - nighttime(1:end-1)) > 5); % look for periods of darkness where next period of darkness isn't for a while 

    figure
    plot(subset.StartDate, subset.rad_sw)
    hold on 
    scatter(subset.StartDate(dawnind), dawnind.*0+50, 50, 'filled')
    
    check = input('look ok? y/n', 's');
    if strcmp(check, 'n')
        keyboard
    end
    
    for dnum = 1:length(dawnind) %now deal with each day individually
        daystarttime = subset.StartDate(dawnind(dnum)) ;
        dayenddtime = datenum(daystarttime)+1.0417; 
        daytable = subset((subset.StartDate >= daystarttime & datenum(subset.StartDate) < dayenddtime), :) ;
        
        t_reldawn = datenum(daytable.StartDate) - datenum(daystarttime); 
        hournum = floor(t_reldawn.*24)+1;
        
        %format light data
        Edata = [t_reldawn.*24 daytable.rad_sw];
        Edata(Edata(:,2) < 0, 2) = 0; %remove negative light values 
        
        synvolbins = 2.^[-5:0.125:2];
        N_dist = nan(57, 25);
        cellsperml = nan(1,25);
        
        %format syn data inputs
        func2 = @(x) sum(x, 1); 
        [Groupnum, id] = findgroups(hournum); 
        
        cellsperml = nan(1, 25);
        N_dist = nan(57, 25);
        cellsperml(id) = splitapply(func2, daytable.Syn_count, Groupnum)./splitapply(func2, daytable.VolAnalyzed_ml, Groupnum); 
        N_dist(:, id) = splitapply(func2, daytable.SynDist, Groupnum)';
        Vhists = N_dist./sum(N_dist); 
        
        save([outpath cruisename filesep cruisename 'day' num2str(dnum) 'input.mat'], 'daystarttime', 'daytable', 'synvolbins', 'N_dist', 'cellsperml', 'Vhists', 'Edata')
    end
    
    
end

