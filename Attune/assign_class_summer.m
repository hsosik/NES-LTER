function [ class, bounds] = assign_class_AR92( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate, pro_measured)
% function [ class, bounds] = assign_class_AR92( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )

fcsdat = array2table(fcsdat, 'VariableNames', {fcshdr.par.name});
%plot_flag = 0; 

%sampled alternating from start of cruise until 
%sampled pro from Nov 1 1714? BL = 390 or 370 until 11:15 Nov 2 (local), experiment 01Nov(19) BL hv
%= 340; If BL3 voltage is >340, it is a Pro sample and needs pro gate.
%Other gates will need to move as well on CHl channel.

%for files that have pro on AE2426 - used date range in the past, hopefully
%this works on triggersettings, it did for AR82
% if endsWith(fcshdr.tr1_par, "_SSC") && fcshdr.tr1_level < 500
%     phase = 2;
%     pro_measured =1;
% else
% %regular shelf settings
% phase =1;
% pro_measured =0;
% end

%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    par_eukX = 'BL3-H'; %15; %BL3-H %chlorophyll 
    par_eukY = 'GL2-H'; %18; %GL2-H %phycoerythrin 
    
    %parameter numbers for main syn polygon
    par_synX = 'SSC-H'; %12; %11 is FSC-H, 12 is SSC-H
    par_synY = 'GL2-H'; %18; %GL2-H %phycoerythrin 
    
    %just for initial gates
    synmaxY = 4e5; 
    synminX = 10 ; 
    synXcorners = [7000 30000]; 

    eukminX = 1000;% Heidi Sep 2026--this seems WAY too high for HRS2606... 3.5e3; 
    eukcorner = [20000 1200]; 
    eukmaxY = 4e4; 
    eukmaxYlower = 400;%300; 
    gl2_noise_thresh = 1000; %basically synminY
    
    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 1; %1.1; %PE to CHL, ?1.2 with .8 offset?
    synGL1H2BL3Hoffset = .3; %.4; %PE to CHL .4 on RB, .3 on TN, .8?
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -2.5; %-.6; %PE to CHL -.8 on TN
    nonsynfactorA = 25; %6
    nonsynfactorB = 6; %2.5

    if pro_measured %phase == 2
    eukminX =  8000;%2000;%5000;%5e3; 
    eukcorner = [30000 900]; %EP had at 10000...but needs to go up instead [20000 1200]; 
    eukmaxY = 4e4; 
    eukmaxYlower = 400;%300; 
    gl2_noise_thresh = 1000;%10000; %basically synminY 
    synGL1H2BL3Hoffset = -.2; %PE to CHL .4 on RB, .3 on TN, .8?
    end
  
    %syn main gate
    gsyn_main_gate = [synminX gl2_noise_thresh ; synXcorners(1) gl2_noise_thresh; synXcorners(2) synmaxY; synminX synmaxY]; %[Xmin Ymin; Xmax Ymax]
    %euk gate 
    geuk_main_gate = [eukminX eukmaxYlower;  eukcorner(1) eukcorner(2); 1100000 eukmaxY; 1100000 1; eukminX 1];

    % 
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    fcsdatlog{:,:} = real(fcsdatlog{:,:});
    a = fcsdatlog{:,:}; a(isinf(a)) = log10(1); %get rid of the -inf from raw = 0
    fcsdatlog{:,:} = a;
    in_euk = inpolygon(fcsdatlog.(par_eukX),fcsdatlog.(par_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog.(par_synX),fcsdatlog.(par_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
    
    %first look in gates, then cut out extremes? or add if pileup at edges
    %minX = prctile(fcsdat{in_syn,par_synX},10)*.3; maxX = prctile(fcsdat{in_syn, par_synX}, 90)*10; 
    %minY = prctile(fcsdat{in_syn,par_synY},10)*.3; maxY = prctile(fcsdat{in_syn,par_synY},90)*10;
    minX = prctile(fcsdat{in_syn,par_synX},5)*.5; maxX = prctile(fcsdat{in_syn, par_synX}, 90)*10; 
    minY = prctile(fcsdat{in_syn,par_synY},5)*.3; maxY = prctile(fcsdat{in_syn,par_synY},90)*10;
 
    %eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    %eukminX = prctile(fcsdat{in_euk,par_eukX},10)*.3;
    %eukminX = max([eukminX 500]); %Pretty sure its always eukminX
    b = 3:.2:4.5; %log bin edges
    [m f] = mode(discretize(log10(fcsdat.(par_eukX)(in_euk)), b,(b(1:end-1))+.05)); %log mode
%    eukminX = prctile(fcsdat{in_euk&fcsdatlog.(par_eukX)>(m-.5),par_eukX},2)/2;
    eukminX = prctile(fcsdat{in_euk&fcsdatlog.(par_eukX)>(m-.5),par_eukX},2)/5;  %get too much then refine with mahal below

    minY = max([minY 100]); %not below trigger level for this cruise

    %make new gates with adapted boundaries
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    gsyn_main_gate(1,1) = minX; 
    gsyn_main_gate(4,1) = minX; 
    geuk_main_gate(1,1) = eukminX; 
    geuk_main_gate(5,1) = eukminX; 

    % assingments for euk, syn and pro 
    in_euk = inpolygon(fcsdatlog.(par_eukX),fcsdatlog.(par_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2)));
    in_syn = (inpolygon(fcsdatlog.(par_synX),fcsdatlog.(par_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
    
    in_euk = find(in_euk); 
    eukdist = mahal(fcsdatlog{in_euk,{'SSC-H' 'BL3-A'}},fcsdatlog{in_euk,{'SSC-H' 'BL3-A'}});
    %ind = eukdist<10;
    %ind_tight = eukdist<8 & fcsdatlog.('BL3-H')(in_euk)<yt & fcsdatlog.('SSC-H')(in_euk)<xt;
    yt = prctile(fcsdatlog.('BL3-H')(in_euk(eukdist<10)),10);
    xt = prctile(fcsdatlog.('SSC-H')(in_euk(eukdist<10)),10);
    ind = fcsdatlog.('BL3-H')(in_euk)<yt+1 & fcsdatlog.('SSC-H')(in_euk)<xt+1; %bottom decades
    eukdist = mahal(fcsdatlog{in_euk,{'SSC-H' 'BL3-A'}},fcsdatlog{in_euk(ind),{'SSC-H' 'BL3-A'}});
    ind = eukdist<3 & fcsdatlog.('BL3-H')(in_euk)<yt+1 & fcsdatlog.('SSC-H')(in_euk)<xt+1; %bottom decades
    eukdist = mahal(fcsdatlog{in_euk,{'SSC-H' 'BL3-A'}},fcsdatlog{in_euk(ind),{'SSC-H' 'BL3-A'}});
    yt = prctile(fcsdatlog.('BL3-H')(in_euk(eukdist<10)),50);
    xt = prctile(fcsdatlog.('SSC-H')(in_euk(eukdist<10)),50);
   
    in_euk(eukdist>20 & fcsdatlog.('BL3-H')(in_euk)<yt ) = []; %& fcsdatlog.('SSC-H')(in_euk)<xt) = []; %not euks 

    if pro_measured %exist("pro_main_gate", "var")

        prominX = 200;
        promaxX = eukminX; %/2; %not higher than 2x bottom of Euks on chl %4000;
        prominY = 0;
        promaxY = 400;
        pro_main_gate = [prominX promaxY;  prominX prominY; promaxX prominY; promaxX promaxY]; %gates Pro on GL2/BL3 plot
        prominX2 = 200;
        promaxX2 = 4000;
        prominY2 = 0;
        promaxY2 = minY; %bottom on syn on PE % 400; 
        pro_main_gate_PEvsSSC = [prominX2 promaxY2;  prominX2 prominY2; promaxX2 prominY2; promaxX2 promaxY2]; %gates Pro on GL2/BL3 plot
        pro_2nd_gate = [400 0; 400 800; 8000 8000; 8000 0]; %gates Pro on GL2/SSC plot

        %in_pro_chl =inpolygon(fcsdatlog.(par_eukX),fcsdatlog.(par_eukY),log10(pro_main_gate(:,1)),log10(pro_main_gate(:,2)));
        %in_pro_ssc = inpolygon(fcsdatlog.(par_synX),fcsdatlog.(par_synY),log10(pro_main_gate_PEvsSSC(:,1)),log10(pro_main_gate_PEvsSSC(:,2)));
        in_pro_chl = fcsdat.(par_eukX)>prominX & fcsdat.(par_eukX)<promaxX & fcsdat.(par_eukY)>prominY & fcsdat.(par_eukY)<promaxY ;
        in_pro_ssc = fcsdat.(par_synX)>prominX2 & fcsdat.(par_synX)<promaxX2 & fcsdat.(par_synY)>prominY2 & fcsdat.(par_synY)<promaxY2;
        in_pro = in_pro_chl & in_pro_ssc;
      
        %in_pro = single(in_pro);
        %disregard pro gating if it is spread out along the scatter channel,
        %probably mostly detritus
        if sum(in_pro)/sum(in_pro_chl) < .8
            %in_pro = ~logical(fcsdatlog(:,npar_synX));
            in_pro(1:length(class), :) = logical(0);
        else
           %one more step to clean up pro cluster
           prodist = mahal(fcsdatlog{in_pro,{'SSC-H' 'BL3-A'}},fcsdatlog{in_pro,{'SSC-H' 'BL3-A'}});
           prodist = mahal(fcsdatlog{in_pro,{'SSC-H' 'BL3-A'}},fcsdatlog{in_pro(prodist<2),{'SSC-H' 'BL3-A'}});
           in_pro = find(in_pro); in_pro = in_pro(prodist<12); %3 and 8
        end
    else
        in_pro(1:length(class), :) =logical(0);
    end
  

    %% Part 2
    %it would be really nice if we could adjust the diagonal line in the
    %Chl PE relationship to move with the data
%    frac_coinc = sum(in_syn & (fcsdatlog.(par_synY)<fcsdatlog.("BL3-H")*synGL1H2BL3Hslope+synGL1H2BL3Hoffset))./sum(in_syn);
%    while frac_coinc > .03
%        synGL1H2BL3Hoffset = synGL1H2BL3Hoffset - .1;
%        frac_coinc = sum(in_syn & (fcsdatlog.(par_synY)<fcsdatlog.("BL3-H")*synGL1H2BL3Hslope+synGL1H2BL3Hoffset))./sum(in_syn);
%    end
    %okay, let's try just adjusting the intercept according to the Syn peak location
    ymed = prctile(fcsdatlog.("GL2-H")(in_syn),50);
    xmed = prctile(fcsdatlog.("BL3-H")(in_syn),50);
    synGL1H2BL3Hoffset = ymed-(xmed+.6)*synGL1H2BL3Hslope; 
    
    %% part 3

    %look for things with low syn level phycoerythrin & low GL2/GL3 ratio?
    %& not big FCS with low phycoerythrin
    in_nonsyn_lowPE = fcsdat.(par_synY) > minY & fcsdat.(par_synY) < maxY/2 & fcsdat.(par_synY)./fcsdat.("GL3-H") < nonsynfactorB & ~(fcsdat.(par_synY)<1e4 & fcsdat.("FSC-H")>1e4);
    in_nonsyn_hiPE = fcsdat.(par_synY) > maxY/2 & fcsdat.(par_synY)./fcsdat.("GL3-H") > nonsynfactorB & fcsdat.(par_synY)./fcsdat.("GL3-H") < nonsynfactorA;
    
    %assign values in class vector
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;
    %keyboard

    class(in_syn) = 2; %must be done after nonsyn

    %now use diagonal line in plot 1 to distinguish syn from euks and coincident
    class((fcsdatlog.(par_synY)<fcsdatlog.("BL3-H")*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat.("BL3-H")> min(geuk_main_gate(:,1))) & fcsdat.(par_synY)> minY) = 3; %more euks

    class(in_euk) = 1; %AFTER lowPE

    %rule out eukaryotes classed with smaller ssc than syn minimum
    class(fcsdat.(par_synX)<minX) = 0; 
        
    %classify pro
    class(in_pro) = 7;
   
    %use size to chl ratio to rule out noise between syn and euks
    class(class ~= 2 & fcsdatlog.(par_eukX)<fcsdatlog.(par_synX)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 0; %more noise
    
    %find "syn, euk coincident" mean CHL
    %look for "low PE euks" with less chlorophyll, and mark them noise
    meancoincX = nanmean(fcsdat{class==6, par_synX});
    class(class == 6 & fcsdat.(par_synX)<meancoincX) = 0; 

    %group things with very high PE signals. 
    in_nonsyn_hiPE  = class == 3 & fcsdat.(par_synY) > 6.5e5;
    class(in_nonsyn_hiPE) = 4;

    class(fcsdat.(par_eukX) < 200 & fcsdat.(par_synY) < 250) = 0; %noise %LAST
    
    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
end
    

    
