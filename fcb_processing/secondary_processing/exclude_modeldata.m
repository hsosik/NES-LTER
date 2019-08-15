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
        days2redo=[{'733457'} {'too high winter growth rate?'};
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
                {'734457'} {'Strange SSC pattern for model to handle'}];
        
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
        days2redo=[{'735690'} {'too high winter growth rate?'}];
        days2exclude=[];
        
    case 2015
        days2redo=[{'736057'} {'too high growth rate? some SSC noise problems?'}];
        days2redo=[days2redo; {'736066'} {'too low growth rate? not found global max?'}];
        days2redo=[days2redo; {'736074'} {'too high growth rate? not found global max?'}];
        days2redo=[days2redo; {'736123'} {'too low growth rate? not found global max?'}];
        
        days2exclude=[{'736106'} {'bad SSC in middle of day?'}];
        days2exclude=[days2exclude; {'736060'} {'bad SSC changes in middle of day'}];
        days2exclude=[days2exclude; {'736110'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736112'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736165'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736168'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736177'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736181'} {'bad SSC'}];
        days2exclude=[days2exclude; [cellstr(num2str((736185:736202)')) repmat({'bad SSC'},18,1)]];
        days2exclude=[days2exclude; {'736208'} {'strange SSC change in middle of day'}];
        
    case 2016 
        days2exclude=[{'736545'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736546'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736552'} {'bad SSC'}];
        days2exclude=[days2exclude; {'736600'} {'bad classification?'}];
        days2exclude=[days2exclude; {'736605'} {'bad data'}];
        days2exclude=[days2exclude; {'736606'} {'bad data'}];
        days2exclude=[days2exclude; {'736613'} {'bad data'}];
        days2exclude=[days2exclude; {'736614'} {'bad data'}];
        days2exclude=[days2exclude; {'736693'} {'bad classification?'}];
        days2redo=[];
    case 2017 
        days2exclude=[{'737012'} {'bad SSC for morning'}];
        days2redo=[];
    case 2018 
        days2exclude=[{'737158'} {'too few cells and bad classification'}];
        days2exclude=[days2exclude; [cellstr(num2str((737160:737164)')) repmat({'too few cells, bad classificaiton'},5,1)]];
        days2exclude=[days2exclude; [cellstr(num2str((737167: 737175)')) repmat({'too few cells, bad classificaiton'},9,1)]];
        days2exclude=[days2exclude; [cellstr(num2str((737264: 737282)')) repmat({'bad SSC'},19,1)]];
        days2exclude=[days2exclude; {'737307'} {'bad SSC at end of day'}];

        days2redo=[];
end