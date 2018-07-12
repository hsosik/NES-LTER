datapath = '\\sosiknas1\Lab_data\MVCO\FCB\MVCO_jan2014\data\';

flist = dir([datapath 'FCB*']);
flist = {flist.name}';
t = NaN(size(flist));
z5 = t;
z6 = t;
z8 = t;
for ii = 1:length(flist)
    disp(flist{ii})
    [header, dat, hdrt datt] = fcbreadraw2(datapath, flist{ii});
    t(ii) = size(dat,1); 
    z5(ii) = length(find(dat(:,5) == 0));
    z6(ii) = length(find(dat(:,6) == 0));
    z8(ii) = length(find(dat(:,8) == 0));
end