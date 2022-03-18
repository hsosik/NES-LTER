function [ class , bounds] = assign_class_EN655( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )

if startdate < 7.379983992129629e+05
    phase = 1
else
    phase = 2
end

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H %chlorophyll 
    npar_eukY = 18; %GL2-H %phycoerythrin 
    
    %parameter numbers for main syn polygon
    npar_synX = 12; %11 is FSC-H, 12 is SSC-H
    npar_synY = 18; %GL2-H %phycoerythrin 
    
    %just for initial gates
    synmaxY = 2e4; 
    synminX = 200 ; 
    synXcorners = [10000 70000]; 

    eukminX = 600; 
    eukcorner = [7500 800]; 
    eukmaxY = 60000; 
    eukmaxYlower = 400; 

    gl2_noise_thresh = 1500; %basically synminY 
  

    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.3; %PE to CHL, ?1.2 with .8 offset?
    synGL1H2BL3Hoffset = -.7; %PE to CHL .4 on RB, .3 on TN, .8?
    syneukBL3H2SSCHslope = 1.1; %CHL to SSC
    syneukBL3H2SSCHoffset = -1.5; %CHL to SSC
    nonsynfactorA = 15; %6
    nonsynfactorB = 1; %2.5
    
    if phase == 2
       synXcorners = [4000 19000]; 

    end

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
    minY = prctile(fcsdat(in_syn,npar_synY),10)*.3; maxY = prctile(fcsdat(in_syn,npar_synY),90)*10;
 
    %%compare syn mimimum Y to noise level. 
    in_noise_tail = (fcsdat(:,npar_synY)>100 & fcsdat(:,npar_synX)<100);
    minY2 = prctile(fcsdat(in_noise_tail, npar_synY), 99); 
    minY = min(minY, minY2);

    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = max([eukminX 500]); %Pretty sure its always eukminX
    minY = max([minY 100]); %not below trigger level for this cruise
   

    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(4,1) = minX; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 


    % assingments for euk and syn 
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));

    %this cruise has another smaller syn population 
    gsyn_small_gate = [minX  minY; 10 1800; 10 2e4; synXcorners(2) maxY];
    in_small_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_small_gate(:,1)),log10(gsyn_small_gate(:,2))));

    %% Part 2
    %it would be really nice if we could adjust the diagonal line in the
    %Chl PE relationship to move with the data
    
    frac_coinc = sum(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset))./sum(in_syn);
    while frac_coinc > .03
        synGL1H2BL3Hoffset = synGL1H2BL3Hoffset - .1;
        frac_coinc = sum(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset))./sum(in_syn);
    end

    %% part 3

    %look for things with low syn level phycoerythrin & low GL2/GL3 ratio?
    %& not big FCS with low phycoerythrin
    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < maxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorB & ~(fcsdat(:,npar_synY)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) > nonsynfactorB & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorA;
    
    %assign values in class vector
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;
    %keyboard

    class(in_syn) = 2; %must be done after nonsyn

    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class((fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1))) & fcsdat(:, npar_synY)> minY) = 3; %more euks
    class(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1)))) = 5; %Syn&Euk coincident
    %also use FSC-W to see coincidence?

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<minX) = 0; 
   
    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    
    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat(class==6, npar_synX));
    class(class == 6 & fcsdat(:,npar_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY;
    class(in_nonsyn_hiPE) = 4;

    class(fcsdat(:,npar_eukX) < 200 & fcsdat(:,npar_synY) < 250) = 0; %noise %LAST
  
    %weird second cluster of syn 
    class(in_small_syn) = 2;

    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    

    
