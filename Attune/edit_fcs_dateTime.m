%p = 'C:\work\Attune\20221117_AR70B\FCS\';
p = '\\sosiknas1\Lab_data\Attune\cruise_data\20221117_AR70B\FCS\';
flist = dir([p '*.fcs']);
for ii = 1:length(flist)
    %f = 'NESLTER_AR70B_17Nov2022(1)_Group_day0_Sample(1).fcs';
    f = flist(ii).name;
    deltaT = duration(24*5+1,27,0); %clock wrong by 5 days, 1 hour, 27 minutes

    fid = fopen([p f],'r+','b');
    status = fseek(fid,0,'bof');
    fcsheader_1stline = fread(fid,64,'char');
    FcsHeaderStopPos    = str2num(char(fcsheader_1stline(19:26)'));
    status = fseek(fid,0,'bof');
    fcsheader_total = fread(fid,FcsHeaderStopPos+1,'char');%read the total header

    s1 = '$BTIM/'; l1 = 7; %e.g., 21:40:00
    s2 = '$DATE/'; l2 = 10; %e.g., 11-Nov-2022
    s3 = '$ETIM/'; l3 = 7;

    startpos1 = findstr(char(fcsheader_total'),s1)+length(s1);
    btim = char(fcsheader_total(startpos1:startpos1+l1)');
    startpos2 = findstr(char(fcsheader_total'),s2)+length(s2);
    dateStr = char(fcsheader_total(startpos2:startpos2+l2)');
    startpos3 = findstr(char(fcsheader_total'),s3)+length(s3);
    etim = char(fcsheader_total(startpos3:startpos3+l3)');

    dt = datetime([dateStr ' ' btim],'Format', 'dd-MMM-yyyy HH:mm:SS');
    dt_new = dt+deltaT;
    dateStr_new = char(dt_new, 'dd-MMM-yyyy');
    btim_new = char(dt_new, 'HH:mm:SS');
    dt = datetime([dateStr ' ' etim],'Format', 'dd-MMM-yyyy HH:mm:SS');
    dt_new = dt+deltaT;
    etim_new = char(dt_new, 'HH:mm:SS');

    status = fseek(fid,startpos1-1,'bof');
    fwrite(fid,uint8(btim_new));
    status = fseek(fid,startpos2-1,'bof');
    fwrite(fid,uint8(dateStr_new));
    status = fseek(fid,startpos3-1,'bof');
    fwrite(fid,uint8(etim_new));

    fclose(fid);

end