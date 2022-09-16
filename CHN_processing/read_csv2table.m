p = 'C:\work\LTER\POC\';
C_atomic_wt = 12.0107; %g*mol^-1,  atomic weight of carbon
N_atomic_wt = 14.0067; %g*mol^-1,  atomic weight of nitrogen

CHNtable = readtable([p 'NESLTER_chn_output.csv']);
CHNtable.datetime = datetime(CHNtable.date, 'InputFormat', 'yyyy-MM-dd HH:mm:ss+00:00');
CHNtable.POC_umolperL = CHNtable.umol_C./CHNtable.CHNVol*1000; %micromolar = micromol per liter
CHNtable.PON_umolperL = CHNtable.umol_N./CHNtable.CHNVol*1000; %micromolar = micromol per liter
CHNtable.POC_ugperL = CHNtable.umol_C*C_atomic_wt./CHNtable.CHNVol*1000; %microgram per liter
CHNtable.PON_ugperL = CHNtable.umol_N*N_atomic_wt./CHNtable.CHNVol*1000; %microgram per liter
CHNtable.C_to_N_molar_ratio = CHNtable.umol_C./CHNtable.umol_N; %mol_per_mol

%POCii = find(CHNtable.depth<5);
save([p 'NESLTER_CHN_table'],'CHNtable')
