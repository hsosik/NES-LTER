function [ class, bounds ] = assign_class_spiropa_TN368( fcsdat, fcshdr, plot_flag, filename, QC_flag, startdate )
%function [ class ] = assign_class_spiropa( fcsdat, fcshdr, plot_flag )
%
%   input:  data and hdr from FCS file (fcsdat, fcshdr)
%           boolean flag for plotting of cytograms
%   output: vector of integers specifying classes These gates configured
%   for SPIROPA cruise AR29, April 2018
%

if startdate < 7.376119368055556e+05
    phase = 1;
else 
    phase = 2; 
end
if startdate > 7.376122708333334e+05 & startdate  < 7.376123333333334e+05
    phase = 3; 
end
if startdate > 7.376133041666667e+05 & startdate < 7.376133305555555e+05; 
    phase =4;
end
if startdate > 7.376136805555555e+05
    phase =5;
end
if (startdate > 7.376138583333333e+05 & startdate < 7.376139458333333e+05) | startdate > 7.376140625000000e+05 & startdate < 7.376141770833334e+05; 
    phase = 6;
end
if startdate > 7.376142361111111e+05
    phase = 6; 
end
if startdate > 7.376144861111111e+05
    phase = 7; 
end
if startdate > 7.376146805555555e+05 
    phase = 2;
end
if startdate > 7.376148173611110e+05
    phase = 6;
end
if startdate > 7.376148722222223e+05 & startdate < 7.376152916666666e+05
    %phase = 1; 
end
if startdate > 7.376152986111111e+05
    phase = 5; 
end
if startdate > 7.376155361111110e+05
    phase =0;
end


%Initialze class vector
    class = zeros(size(fcsdat,1),1);
    
    %parameter numbers for main euk polygon
    npar_eukX = 15; %BL3-H
    npar_eukY = 19; %GL1-H
    
    %parameter numbers for main syn polygon
    npar_synX = 11; %FSC-H
    npar_synY = 19; %GL1-H

    synmaxY = 4e5; 
    synminX = 0; 
    synXcorners = [800 800];

    gl1_noise_thresh = 1200; 

    nonsynfactorB = 2.5; 
    nonsynfactorA = 14; 

    eukminX = 500;
    geuk_main_gate = [eukminX 1100; 8e3 1100; 1100000 1e5; 1100000 1; eukminX 1];
    
    synGL1A2GL1Hmax = 4; %PE area to height
    synGL1H2BL3Hslope = 2.0; %PE to CHL, ?1.2 with .8 offset?
    synGL1H2BL3Hoffset = -1.9; %PE to CHL .4 on RB, .8?
    syneukBL3H2SSCHslope = 1.3; %CHL to SSC
    syneukBL3H2SSCHoffset = -.8; %-.6; %PE to CHL


    gsyn_main_gate = [synminX gl1_noise_thresh ; synXcorners(1) gl1_noise_thresh; synXcorners(2) synmaxY; synminX synmaxY]; %[Xmin Ymin; Xmax Ymax]

    
    if startdate > datenum(2019,5,17,22,40,0) & startdate < datenum(2019,5,20,19,50,0) & ~strmatch('SPIROPA_RB1904_Grazer5',filename) | strmatch('SPIROPA_RB1904_Grazer7',filename)
        m = 11;
        syneukBL3H2SSCHoffset = syneukBL3H2SSCHoffset-log10(m); %PE to CHL
    end;
    
    %find indices of cells within the gates
    fcsdatlog = log10(fcsdat); %use log10 to make sure inpolygon corresponds to view of polygon on log-log plots
    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));


    minY = prctile(fcsdat(in_syn,19),10)*.3; maxY = prctile(fcsdat(in_syn,19),90)*10;
    eukminX = prctile(fcsdat(in_euk,npar_eukX),10)*.3;
    eukminX = max([eukminX 350]); %not below trigger level = 300 for this cruise
    minY = max([minY gl1_noise_thresh]); %not below trigger level = 1100 for this cruise
    gsyn_main_gate = [0 minY ; 8e3 minY; 8e3 maxY; 0 maxY]; %[Xmin Ymin; Xmax Ymax]
    geuk_main_gate = [eukminX 1100; 8e3 1100; 1100000 1e5; 1100000 1; eukminX 1];
    if (startsWith(filename, 'SPIROPA_TN368_Grazer6')) % ||...
            %startsWith(filename, 'SPIROPA_RB1904_Grazer6') || startsWith(filename, 'SPIROPA_RB1904_Grazer7'))
        minY = 6e3;
        maxY = 8e5;
        if eukminX < 1e4
            eukminX = max([1500 eukminX]);
        else 
            eukminX = 1500;
        end
        geuk_main_gate = [eukminX 1000; 2e4 1000; 1100000 1e5; 1100000 1; eukminX 1];
        gsyn_main_gate = [0 minY ; 8e3 minY; 8e3 maxY; 0 maxY]; %[Xmin Ymin; Xmax Ymax] %maybe only for Grazer4 and after?
    elseif (startsWith(filename, 'SPIROPA_TN368_Grazer8'))% || startsWith(filename, 'SPIROPA_TN368_Grazer10'))
        minY = 1e4;
        maxY = 1.1e6;
        if eukminX < 1e4
            eukminX = max([1800 eukminX]);
        else 
            eukminX = 1800;
        end    
        geuk_main_gate = [eukminX 1000; 2e4 1000; 1100000 1e5; 1100000 1; eukminX 1];
        gsyn_main_gate = [0 minY ; 8e3 minY; 8e3 maxY; 0 maxY]; %[Xmin Ymin; Xmax Ymax] %maybe only for Grazer4 and after?
    elseif (startsWith(filename, 'SPIROPA_TN368_Grazer10') || startsWith(filename, 'SPIROPA_TN368_Grazer12')) 
        minY = 1e4;
        maxY = 1.1e6;
        if eukminX < 1e4
            eukminX = max([1800 eukminX]);
        else 
            eukminX = 1800;
        end    
        geuk_main_gate = [eukminX 1000; 2e4 1000; 1100000 1e5; 1100000 1; eukminX 1];
        gsyn_main_gate = [0 minY ; 1e4 minY; 1e4 maxY; 0 maxY]; %[Xmin Ymin; Xmax Ymax] %maybe only for Grazer4 and after?
    elseif startsWith(filename, 'SPIROPA_TN368_Grazer')
        minY = 3e4;
        maxY = 1.1e6;
        if eukminX < 1e4
            eukminX = max([1500 eukminX]);
        else 
            eukminX = 1500;
        end    
        geuk_main_gate = [eukminX 1000; 2e4 1000; 1100000 1e5; 1100000 1; eukminX 1];
        gsyn_main_gate = [0 minY ; 1e4 minY; 1e4 maxY; 0 maxY]; %[Xmin Ymin; Xmax Ymax] %maybe only for Grazer4 and after?
    end;


      %%compare syn mimimum Y to noise level. 
    in_noise_tail = (fcsdat(:,npar_synY)>100 & fcsdat(:,npar_synX)<15);
    minY2 = prctile(fcsdat(in_noise_tail, npar_synY), 99); 
    minY = min(minY, minY2);
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 


    if phase >1
    %look at minimum density of minimum threshholds
    if sum(in_syn & log10(fcsdat(:, npar_synY))< 4) > 0; 
    [f, xi] = ksdensity(log10(fcsdat(in_syn & log10(fcsdat(:, npar_synY))< 4, npar_synY)));
    if sum(islocalmin(f))>0 & f(islocalmin(f))./sum(f) < .015; 
        minY = 10.^xi(islocalmin(f));
        minY = minY(1)
    end
    end
    gsyn_main_gate(:,2) = [minY; minY; maxY; maxY]; 
    end

    in_euk = inpolygon(fcsdatlog(:,npar_eukX),fcsdatlog(:,npar_eukY),log10(geuk_main_gate(:,1)),log10(geuk_main_gate(:,2))); 
    
    in_syn = (inpolygon(fcsdatlog(:,npar_synX),fcsdatlog(:,npar_synY),log10(gsyn_main_gate(:,1)),log10(gsyn_main_gate(:,2))));
   
   
%% PART 2 
   
    %it would be really nice if we could adjust the diagonal line in the
    %Chl PE relationship to move with the data

    subset = in_syn & log10(fcsdat(:, npar_eukX))>2 & log10(fcsdat(:, npar_eukX))<3;
    lm = fitlm(log10(fcsdat(subset, npar_eukX)), log10(fcsdat(subset, npar_eukY)));
    synGL1H2BL3Hslope = lm.Coefficients.Estimate(2) +.5;
    synGL1H2BL3Hoffset = lm.Coefficients.Estimate(1) - 2;


    if phase == 1 
        gl1_noise_thresh = 700; 
        synGL1H2BL3Hslope = 1.2; 
    elseif phase == 3
       synGL1H2BL3Hslope = 1.8; 
       synGL1H2BL3Hoffset = -1.1; %PE to CHL .4 on RB, .8?
    elseif  phase == 4
        synGL1H2BL3Hslope = 1.78; 
        synGL1H2BL3Hoffset = -1.1; %PE to CHL .4 on RB, .8?
    elseif phase ==5
       synGL1H2BL3Hslope = 1.5; 
        synGL1H2BL3Hoffset = -.1;
    elseif phase ==6
        synGL1H2BL3Hslope = 1.2; 
        synGL1H2BL3Hoffset = .5;
    elseif phase ==7 
        synGL1H2BL3Hslope = 1.2; 
        synGL1H2BL3Hoffset = .2;
    end

   if sum(in_euk)< 1000 %syneukcioncident only show up when euk concentrations are high
        synGL1H2BL3Hslope = 1.2;
        synGL1H2BL3Hoffset = -.8;
    end


  

        %% part 3

    in_nonsyn_lowPE = fcsdat(:,npar_synY) > minY & fcsdat(:,npar_synY) < maxY/2 & fcsdat(:,npar_synY)./fcsdat(:,17) < nonsynfactorB & ~(fcsdat(:,npar_synY)<1e4 & fcsdat(:,11)>1e4);
    in_nonsyn_hiPE = fcsdat(:,npar_synY) > maxY/2 & fcsdat(:,19)./fcsdat(:,18) > nonsynfactorB & fcsdat(:,19)./fcsdat(:,18) < nonsynfactorA;
        
    %assign values in class vector
    
    class(in_nonsyn_lowPE) = 3;
    class(in_nonsyn_hiPE) = 4;
    %check that euks with PE have appropriate chl signals
    class((class == 3 | class ==4) & fcsdat(:,npar_eukX)< eukminX) = 0; 
    class(in_syn) = 2; %must be done after nonsyn
    
    
    class(in_syn & (fcsdatlog(:,19)<fcsdatlog(:,15)*synGL1H2BL3Hslope+synGL1H2BL3Hoffset & fcsdat(:,15)> min(geuk_main_gate(:,1)))) = 5; %Syn&Euk coincident

    %class(find(in_nonsyn_lowPE & fcsdat(:,12)<300)) = 5; %mysterious cells with PE and SSC of Syn but high CHL
    class(in_euk) = 1; %AFTER lowPE

 

    class(class == 0 & fcsdat(:,19)>minY & fcsdat(:,19)<maxY & fcsdatlog(:,15)>fcsdatlog(:,12)*syneukBL3H2SSCHslope+syneukBL3H2SSCHoffset) = 3; %euk with PE
    class(fcsdat(:,15) < 200 & fcsdat(:,19) < 250) = 7; %noise %LAST
    

    %save gate boundaries to pass to moviemaker 
    bounds = {geuk_main_gate, gsyn_main_gate, synGL1H2BL3Hslope, synGL1H2BL3Hoffset, synGL1A2GL1Hmax, syneukBL3H2SSCHslope, syneukBL3H2SSCHoffset, nonsynfactorA, nonsynfactorB}; 
 

end
    
