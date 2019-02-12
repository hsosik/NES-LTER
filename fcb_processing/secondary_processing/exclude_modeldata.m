function [days2redo, days2exclude]=exclude_modeldata(year2do)
%division rate estimates to exclude for now or rerun to make sure:

%days commented out have been done - either incoporated or did not make a
%difference!

switch year2do
    
    case 2003
        days2redo=[{'731750'} {'may not have found global max'}];
        days2redo=[days2redo; {'731751'} {'too high growth rate? may not have found global max'}];
%         days2redo=[days2redo; {'731753'} {'may not have found global max'}];
%         days2redo=[days2redo; {'731755'} {'too high growth rate?'}];
%         days2redo=[days2redo; {'731794'} {'too high growth rate?'}];
%         days2redo=[days2redo; {'731787'} {'too high growth rate?'}];
        
        days2exclude=[];

    case 2004
%         days2redo=[{'732112'} {'may not have found global max'}];
%         days2redo=[days2redo; {'732118'} {'may not have found global max'}];
%         days2redo=[days2redo; {'732119'} {'may not have found global max'}];
%         days2redo=[days2redo; {'732164'} {'may not have found global max'}];
%         days2redo=[days2redo; {'732167'} {'too high growth rate?'}];
        days2redo=[];
        days2exclude=[];

    case 2005
%         days2redo=[{'732507'} {'Classification or tidal interference'}];
%         days2redo=[days2redo; {'732512'} {'Classification or tidal interference'}];
%         days2redo=[days2redo; {'732497'} {'Classification or tidal interference'}];
%         days2redo=[days2redo; {'732485'} {'Classification or tidal interference'}];
        days2redo=[];
        days2exclude=[];
   
    case 2006
        days2redo=[{'732832'} {'Tail interference - may not have found global max'}];
        days2redo=[days2redo; {'732837'} {'Tail interference - may not have found global max'}];
        days2redo=[days2redo; {'732838'} {'Tail interference - may not have found global max'}];
        days2redo=[days2redo; {'732922'} {'may not have found global max'}];
        
        days2exclude=[{'732845'} {'Classification interference in middle of day'}];
        days2exclude=[{'732978'} {'Tidal or classificaiton interference?'}];
    
    case 2007
        days2redo=[{'733398'} {'too high winter growth rate - may not have found global max'}];
%         days2redo=[{'733137'} {'may not have found global max'}];
%         days2redo=[days2redo; {'733136'} {'may not have found global max'}];
%         days2redo=[days2redo; {'733133'} {'may not have found global max'}];
%         days2redo=[days2redo; {'733132'} {'may not have found global max'}];
%         days2redo=[days2redo; {'733131'} {'may not have found global max'}];
%         days2redo=[days2redo; {'733387'} {'may not have found global max'}];
%         
        days2exclude=[{'733131'} {'Classification interference?'}];
    
    case 2008
        days2redo=[{'733457'} {'too high winter growth rate?'}];
        %days2redo=[];
        days2exclude=[{'733702'} {'Classification interference in middle of day?'}];
    
    case 2009
       days2redo=[{'734132'} {'too high winter growth rate?'}];
       % days2redo=[];
        days2exclude=[{'733827'} {'Classification interference in middle of day?'}];
        days2exclude=[days2exclude; {'733828'} {'Classification interference in middle of day?'}];
    
    case 2010
        days2redo=[{'734272'} {'too high growth rate? found global max?'}];
        days2redo=[days2redo; {'734456'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734457'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734472'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734408'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734388'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734389'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734390'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734393'} {'growth rate too high? - may not have found global max'}];
        days2redo=[days2redo; {'734394'} {'growth rate too high? - may not have found global max'}];
        
        days2exclude=[];
    
    case 2011
        days2redo=[{'734590'} {'too high winter growth rate? found global max?'}];
        days2redo=[days2redo; {'734594'} {'too high winter growth rate? found global max?'}];
%         days2redo=[{'734539'} {'too high growth rate? found global max?'}];
%         days2redo=[days2redo; {'734531'} {'growth rate too high - may not have found global max'}];
%        days2redo=[];  
        
        days2exclude=[{'734719'} {'bad merging / SSC'}];
        days2exclude=[days2exclude; {'734720'} {'bad merging /SSC'}];
        days2exclude=[days2exclude; {'734718'} {'bad merging /SSC'}];
        days2exclude=[days2exclude; {'734716'} {'bad merging /SSC'}];
        days2exclude=[days2exclude; {'734717'} {'bad merging /SSC'}];
        days2exclude=[days2exclude; {'734714'} {'bad merging /SSC'}];
        days2exclude=[days2exclude; {'734862'} {'bad SSC'}];
        days2exclude=[days2exclude; {'734863'} {'bad SSC'}];
    
    case 2012
        
        days2redo=[{'734948'} {'too high growth rate? found global max?'}];
        days2redo=[days2redo; {'734947'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734933'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734929'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734919'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734917'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734931'} {'growth rate too low - may not have found global max'}];
        days2redo=[days2redo; {'734935'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'734948'} {'growth rate too high - may not have found global max'}];
        days2redo=[days2redo; {'735034'} {'growth rate too low - may not have found global max'}];
        
        days2exclude=[];

    case 2013
        days2redo=[{'735265'} {'too high winter growth rate?'}];
%         days2redo=[{'735247'} {'too high growth rate? found global max?'}];
%         days2redo=[days2redo; {'735318'} {'growth rate too high - may not have found global max, bad SSC'}];
%         days2redo=[days2redo; {'735265'} {'growth rate too high - may not have found global max, bad SSC'}];
%         days2redo=[days2redo; {'735324'} {'growth rate too high - may not have found global max, bad SSC'}];
%         days2redo=[days2redo; {'735325'} {'growth rate too high - may not have found global max, bad SSC'}];
        
        days2exclude=[{'735328'} {'bad SSC'}];
        
    case 2014
        days2redo=[{'735690'} {'too high winter growth rate?'}];
%         days2redo=[{'735685'} {'too high growth rate? some SSC noise problems'}];
%         days2redo=[days2redo; {'735693'} {'growth rate too high some SSC noise problems?'}];
%         days2redo=[days2redo; {'735694'} {'growth rate too high some SSC noise problems?'}];
%         days2redo=[days2redo; {'735717'} {'growth rate too high some SSC noise problems?'}];
%         days2redo=[days2redo; {'735725'} {'growth rate too high some SSC noise problems?'}];
%         days2redo=[days2redo; {'735767'} {'growth rate too high - not found a global max?'}];
        
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
        
    case 2016 %none so far!
        days2exclude=[];
        days2redo=[];
    case 2017 %none so far!
        days2exclude=[];
        days2redo=[];
    case 2018 %none so far!
        days2exclude=[];
        days2redo=[];
end