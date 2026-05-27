    function [ class , bounds] = assign_class_EN712( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )

plot_flag = 0; 

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H %chlorophyll 
    npar_eukY = 18; %GL2-H %phycoerythrin 
    
    %parameter numbers for main syn polygon
    npar_synX = 12; %11 is FSC-H, 12 is SSC-H
    npar_synY = 18; %GL2-H %phycoerythrin 
    
    %just for initial gates
    synmaxY = 3.3e5; 
    synminX = 770 ; 
    synXcorners = [4200 15000]; 

    eukminX = 9e3; %5e3; 
    eukcorner = [5e4 800]; 
    eukmaxY = 2e4; 
    eukmaxYlower = 300; 

    gl2_noise_thresh = 13000; %basically synminY 
  

    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.2; %PE to CHL, ?1.2 with .8 offset? was 1.3
    synGL1H2BL3Hoffset = 0; %PE to CHL .4 on RB, .3 on TN, .8? AT46 -.20;
    syneukBL3H2SSCHslope = 1.12; %CHL to SSC
    syneukBL3H2SSCHoffset = -1.1; %-.6; %PE to CHL -.8 on TN
    nonsynfactorA = 15; %6 %Move diagonal lines in GL2H GL3H?
    nonsynfactorB = 3; %2.5
    
    %syn main gate
    gsyn_main_gate = [synminX gl2_noise_thresh ; synXcorners(1) gl2_noise_thresh; synXcorners(2) synmaxY; synminX synmaxY]; %[Xmin Ymin; Xmax Ymax]
    %euk gate 
    geuk_main_gate = [eukminX eukmaxYlower;  eukcorner(1) eukcorner(2); 1100000 eukmaxY; 1100000 1; eukminX 1];
    
        
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
    
    %first look in gates, then cut out extremes? or add if pileup at edges
    minX = prctile(fcsdat(in_syn,npar_synX),10)*.3; maxX = prctile(fcsdat(in_syn, npar_synX), 90)*10; 
    if isnan(minX)  %heidi hack for no points in_syn here **CHECK
        minX = 0;
    end
    minY = prctile(fcsdat(in_syn,npar_synY),10)*.3;
   
    %%compare syn mimimum Y to noise level. 
    in_noise_tail = (fcsdat(:,npar_synY)>100 & fcsdat(:,npar_synX)<100);
    minY2 = prctile(fcsdat(in_noise_tail, npar_synY), 99); 
    minY = min(minY, minY2);

    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = max([eukminX 500]); %Pretty sure its always eukminX
    minY = max([minY 100]); %not below trigger level for this cruise

    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; synmaxY; synmaxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(4,1) = 1800; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 


    % assingments for euk and syn 
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));


    %% part 3

    %look for things with low syn level phycoerythrin & low GL2/GL3 ratio?
    %& not big FCS with low phycoerythrin
    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < synmaxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorB & ~(fcsdat(:,npar_synY)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > synmaxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) > nonsynfactorB & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorA;
    
    %assign values in class vector
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;

    class(in_syn) = 2; %must be done after nonsyn

    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class((fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1))) & fcsdat(:, npar_synY)> minY) = 3; %more euks

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<minX) = 0; 
   
    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
   
    %group things with very high PE signals. 
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > synmaxY;
    class(in_nonsyn_hiPE) = 4;

    class(fcsdat(:,npar_eukX) < 200 & fcsdat(:,npar_synY) < 250) = 0; %noise %LAST
  
   
    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    

    
