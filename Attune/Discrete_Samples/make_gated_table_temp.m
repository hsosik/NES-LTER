function [] = make_gated_table_temp(p,cruise)

p.awspath = strcat(p.basepath, 'FCS', filesep);
p.time_gate_fraction = 1;
load([p.outpath 'FCSfileinfo']);
gated_table = FCSfileinfo;

gated_table.awsfilename(:) = {strcat(cruise, '.aws')};

gated_table = movevars(gated_table,'awsfilename', 'After', 'fcslist');

for i = 1:height(gated_table)
    if ~ismissing(gated_table.awsfilename(i))
        %if we have an aws file proceed with gating
        [fcsdat, fcshdr]  = fca_readfcs([p.fpath gated_table.fcslist{i}]);
        [gate_assignments, polygon_names, polygon_vars, polygon_vals, gate_names, gate_logic_legible, parent_logic] = ApplyAWSgates_hierarchical(strcat(p.awspath, gated_table.awsfilename{i}), fcsdat, fcshdr, p.time_gate_fraction);
        gated_table.gate_names{i} = gate_names;
        gated_table.gate_assignments{i} = gate_assignments;
        gated_table.gate_logic{i} = gate_logic_legible;
        gated_table.polygon_names{i} = polygon_names;
        gated_table.polygon_vars{i} = polygon_vars;
        gated_table.polygon_vals{i} = polygon_vals;
    end
end
gated_table.vol_analyzed_ml = gated_table.vol_analyzed*p.time_gate_fraction/1e6; %account for time gate
no_aws_files = [];
save([p.outpath 'Gated_Table.mat'], 'gated_table', 'no_aws_files')

%%
   cut_off_pro_pop = [];

    % Here is where we decide which gates we are interested in and how we will find them

    % gates_of_interest = {'syn'; 'Syn'; 'Euks'; 'euks'; 'pro'; 'Pro'; 'bacteria'; 'Bacteria'};

    classnames = {'Euks == 1'; 'Syn == 2'; 'Bacteria == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks = 6'};

    class_all = cell(1,height(gated_table));
    syn_conc = nan(height(gated_table),1);
    euk_conc = syn_conc;
    bact_incl_pro_conc = syn_conc;
    pro_conc = syn_conc;
    lp_euk_conc = syn_conc;
    hp_euk_conc = syn_conc;

    for i = 1:height(gated_table)
        filename = gated_table.fcslist{i};
        [fcsdat, fcshdr]  = fca_readfcs([p.fpath filename]);

        a = size(gated_table.gate_assignments{i});

        if max(a) == 0
            class_all{i} = [];
            continue
        end

        %first go through and generate a single class assingment for each
        %particle
        gate_assign_i =gated_table.gate_assignments{i};
        class_i = zeros(1, max(size(gate_assign_i)));
        gate_names = gated_table.gate_names{i};


        % look for euks
        if sum(strcmp(gate_names, 'euk'))
            gate_num_1 = strcmp(gate_names, 'euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1;
        end
        if sum(strcmp(gate_names, 'Euk'))
            gate_num_1 = strcmp(gate_names, 'Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_1, :)) = 1;
        end

        if sum(strcmp(gate_names, 'LowP_Euk'))
            gate_num_5 = strcmp(gate_names, 'LowP_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5;
        end

        if sum(strcmp(gate_names, 'LowPE_Euk'))
            gate_num_5 = strcmp(gate_names, 'LowPE_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_5, :)) = 5;
        end
        if sum(strcmp(gate_names, 'HighP_Euk'))
            gate_num_6 = strcmp(gate_names, 'HighP_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6;
        end
        if sum(strcmp(gate_names, 'HighPE_Euk'))
            gate_num_6 = strcmp(gate_names, 'HighPE_Euk');
            class_i(class_i == 0' & gate_assign_i(gate_num_6, :)) = 6;
        end

        %next syn
        if sum(strcmp(gate_names, 'Syn'))
            gate_num_2 = strcmp(gate_names, 'Syn');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end
        if sum(strcmp(gate_names, 'syn'))
            gate_num_2 = strcmp(gate_names, 'syn');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end
        if sum(strcmp(gate_names, 'SYN'))
            gate_num_2 = strcmp(gate_names, 'SYN');
            class_i(class_i == 0' & gate_assign_i(gate_num_2, :)) = 2;
        end


        %next bacteria, may include prochlorococcus, but those will be
        %reassigned below
        if sum(strcmp(gate_names, 'Bacteria'))
            gate_num_b = strcmp(gate_names, 'Bacteria');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end
        if sum(strcmp(gate_names, 'bacteria'))
            gate_num_b = strcmp(gate_names, 'bacteria');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end
        if sum(strcmp(gate_names, 'BacPro'))
            gate_num_b = strcmp(gate_names, 'BacPro');
            class_i(class_i == 0' & gate_assign_i(gate_num_b, :)) = 3;
        end



        %next assign pro, allowed to overwrite bacteria label
        if sum(strcmp(gate_names, 'pro'))
            gate_num_4 = strcmp(gate_names, 'pro');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Pro'))
            gate_num_4 = strcmp(gate_names, 'Pro');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Euk_sm'))
            gate_num_4 = strcmp(gate_names, 'Euk_sm');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end
        if sum(strcmp(gate_names, 'Proc'))
            gate_num_4 = strcmp(gate_names, 'Proc');
            class_i(logical(gate_assign_i(gate_num_4, :))) = 4;
        end


        %ok now save class assignments

        class_all{i} = class_i;


        %onto calculating concentrations
        concent_i = nan(1, 6);

        sizefrac = gated_table.gate_logic{i};

        if exist('gate_num_1')
            concent_i(1) = sum(class_i == 1)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_2')
            concent_i(2) = sum(class_i == 2)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_b')
            concent_i(3) = sum(class_i == 3 | class_i ==4)./gated_table.vol_analyzed_ml(i); %here is where we combine Bacteria and Syn class assignments!
        end

        if exist('gate_num_4')
            concent_i(4) = sum(class_i == 4)./gated_table.vol_analyzed_ml(i);

            %Here is where we check if pro gate is cut off
            if ~contains(filename, 'pro', 'IgnoreCase', true) %if it's not a pro run specifically but there is a pro gate

                x_ch = strmatch('SSC-H', {fcshdr.par.name});
                y_ch = strmatch('BL3-H', {fcshdr.par.name});

                figure(2)
                [N, X] = hist(fcsdat(class_i==4, y_ch)); %check to see if cutoff on BL3
                hist(fcsdat(class_i==4, y_ch))
                if N(1) < N(2) & N(1) < N(3) & N(1) < N(4) & N(4) < N(3)
                    title('good')
                elseif N(1) < N(2) & N(2) < N(3) & N(3) < N(4)
                    title('good')
                else
                    %l = input('population cut off? y/n', 's');
                    %if l == 'y'
                    concent_i(4) = NaN;
                    title('cut off')
                    cut_off_pro_pop = [cut_off_pro_pop; string(filename)];
                end

            end

        end

        if exist('gate_num_5')
            concent_i(5) = sum(class_i == 5)./gated_table.vol_analyzed_ml(i);
        end

        if exist('gate_num_6')
            concent_i(6) = sum(class_i == 6)./gated_table.vol_analyzed_ml(i);
        end

        %disp(gated_table(i,:))
        %disp(concent_i)

        %now save concentrations
        euk_conc(i) = concent_i(1);
        syn_conc(i) = concent_i(2);
        bact_incl_pro_conc(i) = concent_i(3) ;
        pro_conc(i) = concent_i(4);
        lp_euk_conc(i) = concent_i(5);
        hp_euk_conc(i) = concent_i(6);

        clear gate_num_1 gate_num_2 gate_num_4 gate_num_b gate_assign_i gate_num_5 gate_num_6

        notes = "'Euks == 1'; 'Syn == 2'; 'Bacteria & Pro == 3'; 'Pro == 4'; 'LowPE_Euks = 5'; 'HighPE_Euks - 6'";
        notes{2} = gated_table.awsfilename{i};
        class = class_i';
        save([p.classpath regexprep(gated_table.fcslist{i}, '.fcs', '')], 'class', 'notes') %I think class files wont have volume yet if we don't have bead statistics to calibrate 
        clear class class_i
    end

    gated_table.class = class_all';
    gated_table.Euk_conc = euk_conc;
    gated_table.LowP_Euk_conc = lp_euk_conc;
    gated_table.HighP_Euk_conc = hp_euk_conc;

    gated_table.Syn_conc = syn_conc;
    gated_table.Bact_incl_pro_conc = bact_incl_pro_conc;
    gated_table.Pro_conc = pro_conc;

    save([p.outpath '\Gated_Table.mat'], 'gated_table', 'no_aws_files', 'cut_off_pro_pop');

end

 