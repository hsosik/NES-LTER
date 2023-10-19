%modifying (9/13/23) so that all time gates are ignored and added post-hoc.
% this will allow for uniform volume_analyzed values regardless of whether
% a time gate was used in attune workspace or not. 

% modifying so that output includes all gates, and true/false for every
% paticle and every gate. 

% This function reads the gate information out of the workspace files which
% can be created in Attune software in the form .aws
% it then applies those gates including logic to the fsc files in a folder. 
% 

%Inputs: 
% One .aws filename and location as well as the fcsdat and hdr
% variables which are the result of fca_readfcs() 
    
    % E.g: 
    % awsfilename = '\\sosiknas1\Lab_data\Attune\cruise_data\20200725_EN655\preserved2\Gated_FCS\C01N11_hbac.aws';
    % fcsfile = '\\sosiknas1\Lab_data\Attune\cruise_data\20200725_EN655\preserved2\FCS\NESLTER_EN655_Jul2020_preserved(1)_hbac_SYBR_SSC_C01N11.fcs';
    % [fcsdat, fcshdr]  = fca_readfcs([fcsfile]);
    % ApplyAWSgates(awsfilename, fcsdat, fcshdr); 


% Outputs: 
% gate_assingment is a matrix of size(gatenum, length(fcsdat)) wtih the gate value assigned to each
    % particle. ware.  
    % other outputs have details about the gate structure but are optional  for
    % downstream processing. 


% Written July 2022 by Bethany Fowler

function [gate_assignment, polygon_names, polygon_vars, polygon_vals, gate_names, gate_logic_legible] = ApplyAWSgates(awsfilename, fcsdat, fcshdr)


A = fopen(awsfilename);
C = textscan(A, '%s', 'Delimiter', '<'); 
C = C{1}; 

%Ok region names are within <Regions ... </Regions> 
% that has polygon information in it 
%then <Gates ... </Gates> has gate logic in them. 


%% First let's get the polygon info 

start_index = cellfun(@(x) contains(x,'Regions Next'), C, 'UniformOutput', 1);
end_index = cellfun(@(x) contains(x,'/Regions>'), C, 'UniformOutput', 1);

C_polygons = C(find(start_index):find(end_index));

num_polygons = sum(cellfun(@(x) contains(x,'name="&quot'), C_polygons, 'UniformOutput', 1)); 
polygon_delineations = find(cellfun(@(x) contains(x,'name="&quot'), C_polygons, 'UniformOutput', 1));
polygon_delineations = [polygon_delineations; length(C_polygons)];

polygon_names = cell(1, num_polygons); 
polygon_vars = cell(2, num_polygons); 
polygon_vals = {}; 
for p = 1:num_polygons 
    temp = C_polygons(polygon_delineations(p));
    temp = temp{:}; 

    temp = regexprep(temp, 'bacteria plus', 'bacteria_plus'); 
    temp = regexprep(temp, 'volume to include', 'volume_to_include'); 

    newstr = split(temp, ' ');

    temp = cellfun(@(x) contains(x, 'name="&quot;'), newstr, 'UniformOutput', 1);
    t = strfind(newstr{temp}, '&quot;');
    temp = newstr{temp};
    temp = temp(t(1)+6:t(2)-1); 
    polygon_names{p} = temp; 

    if isempty(newstr(cellfun(@(x) contains(x,'yparam'), newstr, 'UniformOutput', 1)))
        %special case, no y values or parameters
        
        subset = C_polygons(polygon_delineations(p):polygon_delineations(p+1));
    
        temp = newstr(cellfun(@(x) contains(x,'xparam'), newstr, 'UniformOutput', 1));
        temp = temp{:}; 
        ind = strfind(temp, '"');        
        polygon_vars{1, p} = temp(ind(1)+1:ind(2)-1);

    else

    subset = C_polygons(polygon_delineations(p):polygon_delineations(p+1));
    
    temp = newstr(cellfun(@(x) contains(x,'xparam'), newstr, 'UniformOutput', 1));
    temp = temp{:}; 
    polygon_vars{1, p} = temp(9:13);

    temp = newstr(cellfun(@(x) contains(x,'yparam'), newstr, 'UniformOutput', 1));
    temp = temp{:}; 
    polygon_vars{2, p} = temp(9:13);

    end

    xcoord_ind = find(cellfun(@(x) contains(x,'HiResCoord '), subset, 'UniformOutput', 1));
    xcoords = []; 
    ycoords = []; 
    for i = 1:length(xcoord_ind); 
        temp = subset{xcoord_ind(i)}; 

        ind = strfind(temp, '"');
        xcoords = [xcoords; str2num(temp(ind(1)+1:ind(2)-1))]; 
        ycoords = [ycoords; str2num(temp(ind(3)+1:ind(4)-1))]; 

    end
   
    polygon_vals{p} = [xcoords ycoords];

end


clear C_polygons
%% 
%%% Now we have the information about polygons. Need logic for gates. 

start_index = cellfun(@(x) contains(x,'Gates>'), C, 'UniformOutput', 1);
end_index = cellfun(@(x) contains(x,'/Gates>'), C, 'UniformOutput', 1);

C_gates = C(find(start_index):find(end_index));

num_gates = length(C_gates) - 2; 


gate_names = cell(1, num_gates); 
gate_logic_legible = cell(1, num_gates); 
for g = 1:num_gates
    newstr = C_gates{g+1}; 
    newstr = regexprep(newstr, 'bacteria plus', 'bacteria_plus'); 
    newstr = regexprep(newstr, 'volume to include', 'volume_to_include'); 
    newstr = split(newstr, ' ');

    temp = newstr{find(cellfun(@(x) contains(x,'name='), newstr, 'UniformOutput', 1))}; 
    t = strfind(temp, '&quot;');
    gate_names{g} = temp(t(1)+6:t(2)-1); 
 
    temp = newstr{find(cellfun(@(x) contains(x,'equation='), newstr, 'UniformOutput', 1))}; 
    ind = strfind(temp, '"');
    temp = temp((ind(1)+1):(ind(2)-1)); 
    temp = erase(temp, '&quot;');
    temp = regexprep(temp, '&amp;', ' and ');
    temp = regexprep(temp, '!', 'not ');

    %I think it is ok to remove parentheses, since we do not have any OR
    %gates. 
    temp = regexprep(temp, '(', '');
    temp = regexprep(temp, ')', '');
     
    gate_logic_legible{g} = temp; %make a nice legible string for the logic

end

%% go look for parents
parent_logic = cell(1, num_gates); 
for g = 1:num_gates
    newstr = C_gates{g+1}; 
    newstr = split(newstr, ' ');

    temp = newstr{find(cellfun(@(x) contains(x,'parent='), newstr, 'UniformOutput', 1))}; 
    t = strfind(temp, '&quot;');
    if ~isempty(t)
    parent_logic{g} = temp(t(1)+6:t(2)-1); 
    end
end


%% Now we want to apply these polygons and this logic to an fcs file

% Only interseted in tracking certain gates. Might want to adjust this to
% be more modular
gates_to_do = 1:length(gate_names);

gate_assignment = zeros(length(gates_to_do), length(fcsdat)); 

for g = 1:length(gates_to_do)
    parse_logic = gate_logic_legible{gates_to_do(g)}; 
    parse_logic = erase(parse_logic, ' ');
    parse_logic = split(parse_logic, 'and');

    poly_of_int = erase(parse_logic, 'not');

    logic_summary = nan((length(poly_of_int)+1), length(fcsdat)); %check each row 
    for i = 1:length(parse_logic)
        poly_name = poly_of_int(i); 
        poly_num = find(cellfun(@(x) strcmp(x, poly_name) , polygon_names, 'UniformOutput', 1)); 

        %find fcs channels corresponding to parameters of interest 
        npar_x =  strmatch(polygon_vars{1,poly_num}, {fcshdr.par.name});
        
        if ~isempty(polygon_vars{2,poly_num})
            npar_y = strmatch(polygon_vars{2,poly_num}, {fcshdr.par.name});
        
            coords = polygon_vals{poly_num};

            fcsdatlog = log(fcsdat); 
            in_poly = inpolygon(fcsdatlog(:,npar_x), fcsdatlog(:,npar_y), log(coords(:,1)), log(coords(:,2)));
            
            if contains(parse_logic{i}, 'not')
                logic_summary(i, :) = ~in_poly; 
            else
                logic_summary(i, :) = in_poly; 
            end

        else % one dimensional gate
            if npar_x == 1 %if this is a time gate
                logic_summary(i, :) = 1; %put all particles "in this polygon". 
            else %if some other gate, but this probably never happens
                lims = polygon_vals{poly_num};
                lims = lims(:,1);
                in_poly = fcsdat(:,npar_x)< lims(2) & fcsdat(:,npar_x)>lims(1);

                if contains(parse_logic{i}, 'not')
                    logic_summary(i, :) = ~in_poly;
                else
                    logic_summary(i, :) = in_poly;
                end
            end

        end
        clear in_poly

    end 

    %% now add 80% time polygon to all gates

        matdate1 = [fcshdr.date ' ' fcshdr.starttime];
        matdate1 = datenum(matdate1, 'dd-mmm-yyyy HH:MM:SS');
        matdate2 = [fcshdr.date ' ' fcshdr.stoptime];
        matdate2 = datenum(matdate2, 'dd-mmm-yyyy HH:MM:SS');
        fulltime = 1000*etime(datevec(matdate2), datevec(matdate1));

    in_poly = fcsdat(:,1)<= fulltime & fcsdat(:,1)>fulltime*.2; % first 20 percent of time is cut out

    logic_summary(end, :) = in_poly; 

    in_gate = sum(logic_summary, 1)==(length(poly_of_int)+1); %make sure it meets all logic requirements to be in gate
    gate_assignment(g, in_gate) = 1; 

end










