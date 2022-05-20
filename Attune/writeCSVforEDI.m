%Attune2EDI.latitude = Attune_uw_match.gps_furuno_latitude;
%cruiseStr = 'EN668'
cruiseStr = 'EN644'
attunebase = '\\sosiknas1\Lab_data\Attune\cruise_data\';
alist = dir([attunebase '*' cruiseStr]);
attunebase = [attunebase alist.name '\bead_calibrated\'];

load([attunebase 'AttuneTable_uw_match'])

Attune2EDI = table;
Attune2EDI.cruise = repmat(cruiseStr,length(Attune_uw_match{:,1}),1);
Attune2EDI.datetime = datestr(Attune_uw_match.StartDate, 'yyyy-mm-ddTHH:MM:SS+00:00');
Attune2EDI.latitude = Attune_uw_match.gps_furuno_latitude;
Attune2EDI.longitude = Attune_uw_match.gps_furuno_longitude;
Attune2EDI.depth_m = 5*ones(size(Attune_uw_match(:,1)));
Attune2EDI.Orgpicopro_cells_per_ml = Attune_uw_match.(" Syn_count ")./Attune_uw_match.VolAnalyzed_ml;
Attune2EDI.Redpicoeuk_cells_per_ml = Attune_uw_match.(" Euk_count ")./Attune_uw_match.VolAnalyzed_ml;
Attune2EDI.volume_analyzed_ml = Attune_uw_match.VolAnalyzed_ml;
Attune2EDI.filename = Attune_uw_match.Filename;

Attune2EDI(Attune_uw_match.QC_flag==0,:) = [];

writetable(Attune2EDI, [attunebase 'Attune2EDI.csv'])
