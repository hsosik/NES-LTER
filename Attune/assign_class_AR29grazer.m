function [ class, bounds ] = assign_class_AR29grazer( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )


    phase = 1; 

    euk_cluster = 0; %don't look for euk cluster in noisy zones 
    if startdate > 7.371676805555555e+05 & startdate < 7.371785625000000e+05
    euk_cluster = 1;
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
    synmaxY = 53000; 
    synminX = 30 ; 
    synXcorners = [400 10000]; 

    eukminX = 200; 
    eukcorner = [40000 200]; 
    eukmaxY = 2e4; 
    eukmaxYlower = 150; 

    gl2_noise_thresh = 1300; %basically synminY 
  

    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.1; %PE to CHL  EP: This line and the one below it affect the line position on plot one,  which designates coincidence. 
    synGL1H2BL3Hoffset = .2; %PE to CHL 
    syneukBL3H2SSCHslope = 1.2; %CHL to SSC
    syneukBL3H2SSCHoffset = -0.8; %PE to CHL 
    nonsynfactorA = 22; %6
    nonsynfactorB = 1; %2.5


    if phase == 1
    synmaxY = 10000; 
    eukminX = 1200; 
    synXcorners = [400 900]; 

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
    %eukminX = max([eukminX 500]); %Pretty sure its always eukminX
    minY = max([minY 100]); %not below trigger level for this cruise
    

    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(5,1) = minX; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 


    % assingments for euk and syn 
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));

    %% Part 2
    %it would be really nice if we could adjust the diagonal line in the
    %Chl PE relationship to move with the data
    
    frac_coinc = sum(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset))./sum(in_syn)
    while frac_coinc > .03 & synGL1H2BL3Hoffset > 0.1;
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

    class(in_syn) = 2; %must be done after nonsyn

    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class((fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1))) & fcsdat(:, npar_synY)> minY) = 3; %more euks
    class(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1)))) = 5; %Syn&Euk coincident
    %also use FSC-W to see coincidence?

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<minX) = 0; 
    class( in_nonsyn_lowPE & fcsdat(:, npar_synX)<400 ) = 0;, 

    if phase == 1
        %weird low PE group
        in_nonsyn_lowPE = (fcsdat(:,npar_eukY) > eukmaxYlower & fcsdat(:,npar_eukY) <6000 & fcsdat(:,npar_eukX) > 4000); 
        class(in_nonsyn_lowPE & class == 0) = 3; 
    end

    %erroneous syn during flow issues
    class(class == 2 & fcsdat(:,15)< 200) = 0; 

    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    
    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat(class==6, npar_synX));
    class(class == 6 & fcsdat(:,npar_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY;
    class(in_nonsyn_hiPE) = 4;

    if euk_cluster
        %large ssc euk cluster
        pe_euk_main_gate = [3e4 150 ; 3e4 500; 3e5 1000; 3e5 500];
        in_lowPE = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(pe_euk_main_gate(:,1)),log10(pe_euk_main_gate(:,2))));
        
          %first look in gates, then cut out extremes? or add if pileup at edges
          minX = prctile(fcsdat(in_lowPE,npar_synX),10)*.3; maxX = prctile(fcsdat(in_lowPE, npar_synX), 90)*10; 
          minY = prctile(fcsdat(in_lowPE,npar_synY),10)*.3; maxY = prctile(fcsdat(in_lowPE,npar_synY),90)*10;
        
          pe_euk_main_gate = [minX minY; minX 500; maxX maxY; maxX 500];


        class(in_lowPE) = 3;
    end

    class(fcsdat(:,npar_eukX) < eukminX & fcsdat(:,npar_synY) < 250) = 0; %noise %LAST
    class(fcsdat(:,17) < 200) = 0; %use Gl3-H to rule out noise

    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    

    
