 function [ class , bounds] = assign_class_AE2426( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )

plot_flag = 0; 

%sampled alternating from start of cruise until 
%sampled pro from Nov 1 1714? BL = 390 or 370 until 11:15 Nov 2 (local), experiment 01Nov(19) BL hv
%= 340; If BL3 voltage is >340, it is a Pro sample and needs pro gate.
%Other gates will need to move as well on CHl channel.

%for files that have pro on AE2426 - used date range in the past, hopefully
%this works on triggersettings, it did for AR82
if endsWith(fcshdr.tr1_par, "_SSC") && fcshdr.tr1_level < 500
    phase = 2;
else
%regular shelf settings
phase =1;
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
    synmaxY = 4e5; 
    synminX = 10 ; 
    synXcorners = [7000 30000]; 

    eukminX = 5e3; 
    eukcorner = [20000 1200]; 
    eukmaxY = 4e4; 
    eukmaxYlower = 400;%300; 
    gl2_noise_thresh = 10000; %basically synminY
    
    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1.1; %PE to CHL, ?1.2 with .8 offset?
    synGL1H2BL3Hoffset = .4; %PE to CHL .4 on RB, .3 on TN, .8?
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -2.5; %-.6; %PE to CHL -.8 on TN
    nonsynfactorA = 25; %6
    nonsynfactorB = 6; %2.5


    if phase == 2
    eukminX =  5000;%2000;%5000;%5e3; 
    eukcorner = [10000 900]; %[20000 1200]; 
    eukmaxY = 4e4; 
    eukmaxYlower = 400;%300; 
    gl2_noise_thresh = 3000;%10000; %basically synminY 
    synGL1H2BL3Hoffset = -.2; %PE to CHL .4 on RB, .3 on TN, .8?
    end
  
    
    
    % next 4 lines added for case with Pro.
    prominX = 300;  
    promaxX = 3000;
    prominY = 15;
    promaxY = 400;

%     prominX =  2000;%5000;%5e3; 
%     procorner = [10000 600]; %[20000 1200]; 
%     promaxY = 4e4; 
%     promaxYlower = 600;%300; 
    
%     if phase == 2
%         eukminX = 1300; 
%     end


    %syn main gate
    gsyn_main_gate = [synminX gl2_noise_thresh ; synXcorners(1) gl2_noise_thresh; synXcorners(2) synmaxY; synminX synmaxY]; %[Xmin Ymin; Xmax Ymax]
    %euk gate 
    geuk_main_gate = [eukminX eukmaxYlower;  eukcorner(1) eukcorner(2); 1100000 eukmaxY; 1100000 1; eukminX 1];
    

    %pro main gates
    if phase == 2
    pro_main_gate = [prominX promaxY;  prominX prominY; promaxX prominY; promaxX promaxY]; %gates Pro on GL2/BL3 plot
    pro_2nd_gate = [400 0; 400 800; 8000 8000; 8000 0]; %gates Pro on GL2/SSC plot
    end
    
    
        
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
    
    %first look in gates, then cut out extremes? or add if pileup at edges
    minX = prctile(fcsdat(in_syn,npar_synX),10)*.3; maxX = prctile(fcsdat(in_syn, npar_synX), 90)*10; 
    minY = prctile(fcsdat(in_syn,npar_synY),10)*.3; maxY = prctile(fcsdat(in_syn,npar_synY),90)*10;
 
    %eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = max([eukminX 500]); %Pretty sure its always eukminX
    minY = max([minY 100]); %not below trigger level for this cruise


    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(4,1) = minX; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 


%     if phase == 1
%         geuk_main_gate = [geuk_main_gate(1,1) 100; 4000 geuk_main_gate(1,2); geuk_main_gate(2:end,:)]; 
%     else
%         geuk_main_gate = [geuk_main_gate(1,1) 100; eukminX geuk_main_gate(1,2); geuk_main_gate(2:end,:)];
%     end


    % assingments for euk, syn and pro 
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
   if exist("pro_main_gate", "var")
   in_pro_chl =inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(pro_main_gate(:,1)),log10(pro_main_gate(:,2)));
    in_pro_ssc = inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(pro_main_gate(:,1)),log10(pro_main_gate(:,2)));
      in_pro = in_pro_chl & in_pro_ssc;
   %disregard pro gating if it is spread out along the scatter channel,
   %probably mostly detritus
       if sum(in_pro)/sum(in_pro_chl) < .8
        in_pro = ~logical(fcsdatlog(:,npar_synX));
       
       end

%        in_pro = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(pro_main_gate(:,1)),log10(pro_main_gate(:,2)))...
%         & (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(pro_2nd_gate(:,1)),log10(pro_2nd_gate(:,2))));

   else 
       in_pro = ~logical(fcsdatlog(:,npar_synX));
   end

   

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

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat(:, npar_synX)<minX) = 0; 
    
    %classify pro
    class(in_pro) = 7;
   
    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog(:,npar_eukX)<fcsdatlog(:,npar_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    
    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat(class==6, npar_synX));
    class(class == 6 & fcsdat(:,npar_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE  = class == 3 & fcsdat(:,npar_synY) > 6.5e5;
    class(in_nonsyn_hiPE) = 4;

    class(fcsdat(:,npar_eukX) < 200 & fcsdat(:,npar_synY) < 250) = 0; %noise %LAST
  

   
    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    

    
