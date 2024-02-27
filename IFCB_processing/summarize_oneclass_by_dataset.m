function [] = summarize_oneclass_by_dataset(dataset_name, classname, thresh)
% e.g., 
%    dataset_name = 'NESLTER_transect'; %'NESLTER_broadscale';
%    classname = 'Hemiaulus';
%    thresh = 0.9;
% summarize_oneclass_by_dataset(dataset_name, classname, thresh)

base_path = ['\\sosiknas1\IFCB_products\' dataset_name '\summary\'];
%base_path = ['c:\work\IFCB_products\' dataset_name '\summary\'];
if isequal(lower(dataset_name), 'mvco')
    base_path = ['\\sosiknas1\IFCB_products\' dataset_name '\summary_v4\'];
end

outname = [classname '_summary'];

%%yr = 2018:2022;
flist = dir([base_path 'summary_biovol_allHDF_min20_????.mat']);
meta_data = table;
class_summary = table;
%classC_opt = table;
%classC = table;
%classC_adhoc = table;
%classcount_opt = table;
%classcount = table;
%classcount_adhoc = table;

for ii = 1:length(flist)
    disp(flist(ii).name)
    T = load([base_path flist(ii).name], 'meta_data');
    L = load([base_path regexprep(flist(ii).name, '.mat', 'lists.mat')]);
   % primary_score_ind = find(strcmp('score_prim', T2.groupFeaList_variables));
   % second_score_ind = find(strcmp('score_sec', T2.groupFeaList_variables));
   % cellC_ind = find(strcmp('cellC', T2.groupFeaList_variables));
   % ESD_ind = find(strcmp('ESD', T2.groupFeaList_variables));
    %%
    T2 = table;
    T2.ml_analyzed = T.meta_data.ml_analyzed; %initialize the table row number
    for iii = 1:size(T.meta_data,1)
        %cc = T.strmatch(classname, L.class2use);
        %tempFea = L.classFeaList{typeind(ibin),cc};
        %temp = array2table(cat(1, tempFea{:}), 'VariableNames', L.classFeaList_variables);
        %tempFea = L.classFeaList.(classname);
        try 
            temp = array2table(L.classFeaList.(classname){iii}, 'VariableNames', L.classFeaList_variables);
            %ind = (temp.ESD>=ESDmin & temp.score_sec>=thresh);
            ind = (temp.score > thresh);
            T2.C(iii) = sum(temp.cellC(ind)); %./T.ml_analyzed(iii)/1000; %micrograms per liter;
            T2.count(iii) = length(ind);
        catch
            T2.C(iii) = NaN; %./T.ml_analyzed(iii)/1000; %micrograms per liter;
            T2.count(iii) = NaN;
            disp(['missing data for ' L.filelist{iii}])
        end
    end
    class_summary = [class_summary; T2];
    meta_data = [meta_data; T.meta_data];
end

save([base_path  outname], 'meta_data', 'class_summary', 'thresh', 'classname')

disp('Results saved:')
disp([base_path outname])

end




