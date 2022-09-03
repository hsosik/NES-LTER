function [ class, bounds ] = assign_class_spiropa_RB1904( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )
%function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%
%   input:  data and hdr from FCS file (fcsdat, fcshdr)
%           boolean flag for plotting of cytograms
%   output: vector of integers specifying classes These gates configured
%   for SPIROPA cruise AR29, April 2018
%

phase =1 ;

if startdate > 7.375629208333333e+05; %17-May-2019 22:06:00
    phase = 2; 
end
if startdate > 7.375658298611111e+05; % 20-May-2019 19:55:00
    phase =3; 
end
if startdate >  737570; 
    phase = 1
end


%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H
    npar_eukY = 19; %GL1-H
    
    %parameter numbers for main syn polygon
    npar_synX = 12; %SSC-H
    npar_synY = 19; %GL1-H

    synmaxY = 100000; 
    synminX = 50; 
    synXcorners = [1000 800];

    gl1_noise_thresh = 5000; 

    nonsynfactorB = 2; 
    nonsynfactorA = 7.4; 

    eukminX = 200;
    
    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.2; %PE to CHL, ?1.2 with .8 offset?
    synGL1H2BL3Hoffset = 0.4; %PE to CHL .4 on RB, .8?
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -.8; %-.6; %PE to CHL


    if phase ==  2
        synmaxY = 100000; 
        synminX = 300; 
        synXcorners = [5000 20000];
        gl1_noise_thresh = 400; 
        syneukBL3H2SSCHslope = 1.3; %CHL to SSC
        syneukBL3H2SSCHoffset = -2.2 %-.6; %PE to CHL
    end


    if phase ==  3
        synmaxY = 18000; 
        synminX = 1; 
        synXcorners = [5000 20000];
        gl1_noise_thresh = 200; 
        syneukBL3H2SSCHslope = 1.5; %CHL to SSC
        syneukBL3H2SSCHoffset = -2.2 %-.6; %PE to CHL
    end

    geuk_main_gate = [eukminX 210; 10000 310; 1100000 11000; 1100000 1; eukminX 1];

    gsyn_main_gate = [synminX gl1_noise_thresh ; synXcorners(1) gl1_noise_thresh; synXcorners(2) synmaxY; synminX synmaxY]; %[Xmin Ymin; Xmax Ymax]


    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));


    minY = prctile(fcsdat(in_syn,19),10)*.3; maxY = prctile(fcsdat(in_syn,19),90)*10;
    if phase == 3
        minY = 290;
    end; 

    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = max([eukminX 350]); %not below trigger level = 300 for this cruise
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; %[Xmin Ymin; Xmax Ymax]
    geuk_main_gate(1) = eukminX; 
    geuk_main_gate(9) = eukminX; 



    %look at minimum density of minimum threshholds
    if sum(in_euk & log10(fcsdat(:, npar_eukX))< 3.5) > 0; 
    [f, xi] = ksdensity(log10(fcsdat(in_euk & log10(fcsdat(:, npar_eukX))< 3.5, npar_eukX)));
    if sum(islocalmin(f))>0 & f(islocalmin(f))./sum(f) < .015; 
        eukminX2 = 10.^xi(islocalmin(f));
        eukminX2 = eukminX2(1);
        if eukminX2< 2000; 
            eukminX = eukminX2;
        end
    end
    geuk_main_gate(1) = eukminX; 
    geuk_main_gate(9) = eukminX;  


    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
   
   
%% PART 2 
   
    %it would be really nice if we could adjust the diagonal line in the
    %Chl PE relationship to move with the data

    subset = in_syn & log10(fcsdat(:, npar_eukX))>2 & log10(fcsdat(:, npar_eukX))<3;
    lm = fitlm(log10(fcsdat(subset, npar_eukX)), log10(fcsdat(subset, npar_eukY)));
    synGL1H2BL3Hslope = lm.Coefficients.Estimate(2) +.54;
    synGL1H2BL3Hoffset = lm.Coefficients.Estimate(1) - 2.05;


   if sum(in_euk)< 1000 %syneukcioncident only show up when euk concentrations are high
        synGL1H2BL3Hslope = 1.2;
        synGL1H2BL3Hoffset = -.8;
    end


  

        %% part 3

    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < maxY & fcsdat(:,18)./fcsdat(:,17) > nonsynfactorB & fcsdat(:,18)./fcsdat(:,17) < nonsynfactorA & ~(fcsdat(:,npar_synY)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,17) > 10000 & fcsdat(:,18)./fcsdat(:,17) > nonsynfactorA ;   

    %check that euks with PE have appropriate chl signals
    class((class == 3 | class ==4) & fcsdat(:,npar_eukX)< eukminX) = 0; 
    
    
    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class(~in_nonsyn_hiPE & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1))) & fcsdat(:, npar_synY)> minY) = 5; %more euks
    %also use FSC-W to see coincidence?

    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;

    class(in_syn) = 2; %must be done after nonsyn
    class(in_syn & (fcsdatlog(:,npar_synY)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1)))) = 5; %Syn&Euk coincident

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<synminX) = 0; 
   
    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    if phase ==3
        class(fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    end

    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat(class==6, npar_synX));
    class(class == 6 & fcsdat(:,npar_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY;
    class(in_nonsyn_hiPE) = 4;
    class(class == 3 & fcsdat(:,npar_synY)> 9e5) = 4; 

    class(fcsdat(:,npar_eukX) < 200 & fcsdat(:,npar_synY) < 200) = 0; %noise %LAST

    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    
