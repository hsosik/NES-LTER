dataset = 'NESLTER_transect';
metaT = readtable(['https://ifcb-data.whoi.edu/api/export_metadata/' dataset]);
%%
for y = 2017:2026
    disp(y)
    f = ['\\sosiknas1.whoi.edu\IFCB_products\' dataset '\summary\summary_biovol_allHDF_min20_' num2str(y) '.mat'];
    load(f, 'meta_data')
    [ic,ia,ib] = intersect(meta_data.pid, metaT.pid);
    if ~isequal(length(ia), size(meta_data,1))
        keyboard
    end
    if ~isequal(meta_data.skip(ia), metaT.skip(ib))
        keyboard
    end
    if ~isequaln(meta_data(ia,{'skip' 'latitude' 'longitude' 'depth' 'sample_time'}), metaT(ib,{'skip' 'latitude' 'longitude' 'depth' 'sample_time'}))
        keyboard
    end
    %meta_data = removevars(meta_data, 'comment_summary'); meta_data.comment_summary(ia) = metaT.comment_summary(ib);
    if ~isequal(meta_data.ml_analyzed(ia), metaT.ml_analyzed(ib))
        disp('updating ml_analyzed')
        meta_data.ml_analyzed(ia) = metaT.ml_analyzed(ib);
        save(f, '-append', 'meta_data')
    end
end

%%
return
metaT = readtable('https://ifcb-data.whoi.edu/api/export_metadata/NESLTER_broadscale');
metaT.datetime = datetime(metaT.sample_time, 'InputFormat','uuuu-MM-dd HH:mm:ss+00:00');
%%
% load("\\sosiknas1\ifcb_products\NESLTER_broadscale\summary\summary_biovol_allHDF_min20_2017.mat",'meta_data')
ii = find(year(metaT.datetime)==2017);
[ic,ia] = setdiff(metaT.pid(ii), meta_data.pid);
metaT.n_images(ii(ia))
metaT(ii(ia),:)