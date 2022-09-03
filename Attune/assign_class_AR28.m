function [ class, bounds ] = assign_class_AR28( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )

if startdate < 737163 %Apr 13
    phase = 1; 
elseif startdate <  7.371632694444444e+05 %Apr 13 6:28
    phase = 2; 
else 
     phase = 3; 
end

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H %chlorophyll 
    npar_eukY = 19; %GL1-H %phycoerythrin 
    
    %parameter numbers for main syn polygon
    npar_synX = 12; %11 is FSC-H, 12 is SSC-H
    npar_synY = 19; %GL1-H %phycoerythrin 
   
    
    %just for initial gates
    synmaxY = 10000; 
    synminX = 30 ; 
    synXcorners = [400 900]; 

    eukminX = 700; 
    eukcorner = [40000 200]; 
    eukmaxY = 2e4; 
    eukmaxYlower = 150; 

    gl2_noise_thresh = 1300; %basically synminY 
  

    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.18; %PE to CHL
    synGL1H2BL3Hoffset = .2; %PE to CHL 
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -0.8; %PE to CHL 
    nonsynfactorA = 22; %6
    nonsynfactorB = 4; %2.5


    if phase == 1
    synmaxY = 10000; 
    synXcorners = [400 900]; 

    end
    if phase ==3
        eukcorner = [40000 500]; 
        syneukBL3H2SSCHoffset = -1.8; %PE to CHL 
        gl2_noise_thresh = 1900; %basically synminY
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
    
    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    %minY = max([minY 100]); %not below trigger level for this cruise
  
    if phase ==1 
    %look at minimum density of minimum threshholds
    if sum(in_euk & log10(fcsdat(:, npar_eukX))< 4) > 0; 
    [f, xi] = ksdensity(log10(fcsdat(in_euk & log10(fcsdat(:, npar_eukX))< 3.4, npar_eukX)));
    if sum(islocalmin(f))>0 & f(islocalmin(f))./sum(f) < .015; 
        eukminX = 10.^xi(islocalmin(f));
        eukminX = eukminX(1)
    end
    end
    else %super noisy on last day of cruise. 
        eukminX = 700; 
        
    end


    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(5,1) = minX; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 

    if sum(in_euk) ~= 0
    % assingments for euk and syn 
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
    end

    %% part 3
    in_syn(fcsdat(:, npar_synY) < 500) = 0; %really weird picking up syn not in polygon???

    %look for things with low syn level phycoerythrin & low GL2/GL3 ratio?
    %& not big FCS with low phycoerythrin
    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < maxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorB & ~(fcsdat(:,npar_synY)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) > nonsynfactorB & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorA;
    
    %assign values in class vector
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;

    class(in_syn) = 2; %must be done after nonsyn

    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class((fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1))) & fcsdat(:, npar_synY)> minY) = 3; %more euks
    class(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1)))) = 5; %Syn&Euk coincident
    %also use FSC-W to see coincidence?

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<minX) = 0; 
    class( in_nonsyn_lowPE & fcsdat(:, npar_synX)<400 ) = 0;, 

        %weird low PE group
        in_nonsyn_lowPE = (fcsdat(:,npar_eukY) > eukmaxYlower & fcsdat(:,npar_eukY) <6000 & fcsdat(:,npar_eukX) > 4000); 
        class(in_nonsyn_lowPE & class == 0) = 3; 

    %erroneous syn during flow issues
   class(class == 2 & fcsdat(:,15)< 100) = 0; 

    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    
    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat(class==6, npar_synX));
    class(class == 6 & fcsdat(:,npar_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY;
    class(in_nonsyn_hiPE) = 4;

    class(fcsdat(:,npar_eukX) < 150 & fcsdat(:,npar_synY) < 200) = 0; %noise %LAST
    class(fcsdat(:,npar_eukX) <150) = 0; 
   
    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 
end
    

    
