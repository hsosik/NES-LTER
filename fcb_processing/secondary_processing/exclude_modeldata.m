function [days2redo]=exclude_modeldata(year2do)
%division rate estimates to exclude for now or rerun to make sure:

%FOR JUNE 2019 MODEL RUNS!

switch year2do
    
    case 2003
        days2redo=[{'731750'} {'may not have found global max'};
            {'731751'} {'too high growth rate? may not have found global max'};
            {'731755'} {'too high growth rate? may not have found global max'}];
        
    case 2004
        days2redo=[{'732112'} {'bad ssc data to give to model'};
            {'732118'} {'may not have found global max; strange fit from model'};
            {'732119'} {'may not have found global max; strange fit from model'};
            {'732120'} {'may not have found global max; strange fit from model'}];
        
    case 2005
        %really bad SSC smearing starting in Oct 22:
        days2redo=[{'732607'} {'bad ssc data to give to model'};
            {'732610'} {'bad ssc data to give to model'}];
        
    case 2006
        
        days2redo=[{'732837'} {'Tail interference - may not have found global max'};
            {'732838'} {'Tail interference - may not have found global max'};
            {'732840'} {'SSC weirdness / not good model fit'};
            {'732841'} {'SSC weirdness / not good model fit'};
            {'732845'} {'SSC weirdness / not good model fit'};
            {'732846'} {'SSC weirdness / classification / not good model fit'};
            {'732847'} {'SSC weirdness / classification / not good model fit'};
            {'732922'} {'may not have found global max'}];
        
    case 2007
        
        days2redo=[{'733136'} {'may not have found global max'};
            {'733382'} {'bad ssc blip; not okay for model'}
            {'733387'} {'may not have found global max'};
            {'733398'} {'too high winter growth rate - may not have found global max'}];
        
    case 2008
        days2redo=[{'733457'} {'too high winter growth rate?; not found global max'};
            {'733531'} {'too high spring growth rate?; not found global max'};
            {'733536'} {'too high spring growth rate?; not found global max'};
            {'733702'} {'Classification interference in middle of day?'}];
        
    case 2009
        
        days2redo=[{'733801'} {'global max not found?'};
            {'733802'} {'global max not found?'};
            {'733805'} {'global max not found?'};
            {'733806'} {'global max not found?'};
            {'733926'} {'global max not found?'};
            {'734128'} {'global max not found? Not good model fit'};
            {'734132'} {'Odd SSC jump? Not good for model'};
            {'733827'} {'Classification interference in middle of day?'};
            {'733828'} {'Classification interference in middle of day?'}];
        
    case 2010
        
        days2redo=[{'734261'} {'May not have found global max?'};
            {'734372'} {'May not have found global max?'};
            {'734403'} {'May not have found global max?'};
            {'734456'} {'Strange SSC pattern for model to handle'};
            {'734457'} {'Strange SSC pattern for model to handle'};
            {'734272'} {'mu diff with old data; found global max?'};
            {'734336'} {'mu diff with old data; found global max?'};
            {'734345'} {'mu diff with old data; found global max?'};
            {'734372'} {'mu diff with old data; found global max?'};
            {'734373'} {'mu diff with old data; found global max?'};
            {'734378'} {'mu diff with old data; found global max?'};
            {'734379'} {'mu diff with old data; found global max?'};
            {'734380'} {'mu diff with old data; found global max?'};
            {'734381'} {'mu diff with old data; found global max?'};
            {'734382'} {'mu diff with old data; found global max?'};
            {'734383'} {'mu diff with old data; found global max?'};
            {'734384'} {'mu diff with old data; found global max?'};
            {'734389'} {'mu diff with old data; found global max?'};
            {'734390'} {'mu diff with old data; found global max?'};
            {'734391'} {'mu diff with old data; found global max?'};
            {'734392'} {'mu diff with old data; found global max?'};
            {'734400'} {'mu diff with old data; found global max?'};
            {'734403'} {'mu diff with old data; found global max?'};
            {'734407'} {'mu diff with old data; found global max?'};
            {'734411'} {'mu diff with old data; found global max?'};
            {'734428'} {'mu diff with old data; found global max?'};
            {'734432'} {'mu diff with old data; found global max?'};
            {'734456'} {'mu diff with old data; found global max?'};
            {'734457'} {'mu diff with old data; found global max?'};
            {'734472'}  {'mu diff with old data; found global max?'}];
        
    case 2011
        
        days2redo=[{'734646'} {'Did not find global max?'};
            {'734692'} {'Did not find global max?'};
            [cellstr(num2str((734711:734722)')) repmat({'bad SSC'},12,1)]];
        
    case 2012
        
        days2redo=[{'734989'} {'Found global max?'};
            {'735146'} {'Found global max?'};
            {'735149'} {'Found global max?'}];
        
    case 2013
        
        days2redo=[{'735247'} {'Found global max?'};
            {'735265'} {'Found global max?'};
            {'735318'} {'Noise interfering'};
            {'735324'} {'Noise interfering'};
            {'735325'} {'Noise interfering'};
            {'735328'} {'Noise interfering'};
            {'735331'} {'Noise interfering'};
            {'735539'} {'Unfair SSC pattern for model'}];
        
        %[cellstr(num2str((735318:735336)')) repmat({'too few cells plus noise!'},19,1)]];
        
    case 2014
        
        days2redo=[{'735754'} {'Did not find global max?'}];
        
    case 2015
        
        days2redo=[{'736057'} {'Has not found global max?'};
            {'736060'} {'bad SSC changes in middle of day'}
            {'736074'} {'Has not found global max?'};
            {'736086'} {'Has not found global max?'}
            {'736088'} {'Has not found global max?'};
            {'736106'} {'bad SSC in middle of day?'};
            {'736110'} {'bad SSC'};
            {'736112'} {'bad SSC'};
            {'736123'} {'Has not found global max'};
            {'736144'} {'bad SSC'};
            [cellstr(num2str((736164:736216)')) repmat({'bad SSC'},53,1)]]; %so, so very sad :(
        
    case 2016
        
        days2redo=[{'736600'} {'bad SSC'};
            {'736605'} {'bad SSC'};
            {'736606'} {'bad SSC'};
            {'736613'} {'bad SSC'};
            {'736693'} {'bad SSC / classification'}];
        
    case 2017
        
        days2redo=[{'736785'} {'Found global max?'};
            {'736992'} {'Found global max?'}];
        
    case 2018
        
        days2redo=[{'737114'} {'Found global max?'};
            [cellstr(num2str((737160:737163)')) repmat({'too few cells, bad classificaiton'},4,1)];
            [cellstr(num2str((737165:737175)')) repmat({'too few cells, bad classificaiton'},11,1)];
            {'737239'} {'Found global max?'};
            {'737258'} {'Found global max?'};
            {'737267'} {'Fair case for model?'};
            [cellstr(num2str((737271: 737280)')) repmat({'bad SSC'},10,1)];
            {'737307'} {'bad SSC at end of day'}];
        
end