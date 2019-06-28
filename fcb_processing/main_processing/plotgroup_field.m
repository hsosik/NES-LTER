
    disp([filetypelist(1,:) '*'])
    filelist = dir([procpath filetypelist(1,:) '*.mat']);
    for count2 = 2:size(filetypelist,1)
        disp([filetypelist(count2,:) '*'])
        filelist = [filelist; dir([procpath filetypelist(count2,:) '*.mat'])];
    end;
    date = datenum(cat(1,filelist.date));
    [temp, fileorder] = sort(date);
    clear temp date count2
    filelist = filelist(fileorder);
    
    cellCHLall = [];
    cellFLSall = [];
    cellPEall = [];
    cellSSCall = [];
    cellNUMall = [];
    cellresultsall = [];
    beadmatchall = [];
    cellCHLmodeall = [];
    cellFLSmodeall = [];
    cellPEmodeall = [];
    cellSSCmodeall = [];
    fileall = [];
    FCBnumberall = [];

    for filenum = 1:length(filelist),
        filename = filelist(filenum).name;
        disp(filename)
        eval(['load ' procpath filename])
        clear mergedwithclass
        clear ind
        if filename(1:4) == 'FCB2', 
            FCBnumberall = [FCBnumberall; repmat(2, size(cellresults,1),1)]; %FCB2
        else
            FCBnumberall = [FCBnumberall; repmat(1, size(cellresults,1),1)]; %everything else is FCB1
        end;
        fileall = [fileall; repmat(filenum,size(cellresults,1),1) (1:size(cellresults,1))'];
        cellCHLall = [cellCHLall; cellCHL(1:size(cellresults,1),:)];
        cellFLSall = [cellFLSall; cellFLS(1:size(cellresults,1),:)];
        cellPEall = [cellPEall; cellPE(1:size(cellresults,1),:)];
        cellSSCall = [cellSSCall; cellSSC(1:size(cellresults,1),:)];
        cellNUMall = [cellNUMall; cellNUM(1:size(cellresults,1),:)];
        cellresultsall = [cellresultsall; cellresults];
        beadmatchall = [beadmatchall; beadmatch];
        cellCHLmodeall = [cellCHLmodeall; cellCHLmode(1:size(cellresults,1),:)];
        cellFLSmodeall = [cellFLSmodeall; cellFLSmode(1:size(cellresults,1),:)];
        cellPEmodeall = [cellPEmodeall; cellPEmode(1:size(cellresults,1),:)];
        cellSSCmodeall = [cellSSCmodeall; cellSSCmode(1:size(cellresults,1),:)];
    
    end;
    
    
    [~,sind] = sort(cellresultsall(:,1)); cellresultsall = cellresultsall(sind,:);
    cellCHLall = cellCHLall(sind,:);
    cellFLSall = cellFLSall(sind,:);
    cellPEall = cellPEall(sind,:);
    cellSSCall = cellSSCall(sind,:);
    cellNUMall = cellNUMall(sind,:);
    beadmatchall = beadmatchall(sind,:);
    cellCHLmodeall = cellCHLmodeall(sind,:);
    cellFLSmodeall = cellFLSmodeall(sind,:);
    cellPEmodeall = cellPEmodeall(sind,:);
    cellSSCmodeall = cellSSCmodeall(sind,:);
    
    clear fileinfo* filename filenum filelist fileorder filetypelist
    clear cellCHL cellFLS cellPE cellSSC cellNUM cellresults beadmatch *mode
    clear link*
    
    %combine split hourly bins
    dv = datevec(cellresultsall(:,1));
    t = find(diff(dv(:,4))==0 & diff(cellresultsall(:,1)) < 1/24);
    temp = mean([cellresultsall(t-1,1) cellresultsall(t+2,1)],2);
    if ~isempty(find(temp - cellresultsall(t,1) > 1/24)),
        disp('time combine did not work as expected')
        keyboard
    end;
    cellresultsallA = cellresultsall;
    cellresultsall(t,1) = temp;
    cellresultsall(t,2:3) = cellresultsall(t,2:3) + cellresultsall(t+1,2:3);
    temp = cellNUMall(t,:) + cellNUMall(t+1,:);
    cellCHLall(t,:) = (cellCHLall(t,:).*cellNUMall(t,:) + cellCHLall(t+1,:).*cellNUMall(t+1,:)) ./ temp;
    cellFLSall(t,:) = (cellFLSall(t,:).*cellNUMall(t,:) + cellFLSall(t+1,:).*cellNUMall(t+1,:)) ./ temp;
    cellPEall(t,:) = (cellPEall(t,:).*cellNUMall(t,:) + cellPEall(t+1,:).*cellNUMall(t+1,:)) ./ temp;
    cellSSCall(t,:) = (cellSSCall(t,:).*cellNUMall(t,:) + cellSSCall(t+1,:).*cellNUMall(t+1,:)) ./ temp;
    cellNUMall(t,:) = temp;
    beadmatchall(t,:) = (beadmatchall(t,:) + beadmatchall(t+1,:))/2;
    %approximation for modes...
    cellCHLmodeall(t,:) = (cellCHLmodeall(t,:) + cellCHLmodeall(t+1,:))/2;
    cellFLSmodeall(t,:) = (cellFLSmodeall(t,:) + cellFLSmodeall(t+1,:))/2;
    cellPEmodeall(t,:) = (cellPEmodeall(t,:) + cellPEmodeall(t+1,:))/2;
    cellSSCmodeall(t,:) = (cellSSCmodeall(t,:) + cellSSCmodeall(t+1,:))/2;
    
    cellresultsall(t+1,:) = [];
    cellCHLall(t+1,:) = [];
    cellFLSall(t+1,:) = [];
    cellPEall(t+1,:) = [];
    cellSSCall(t+1,:) = [];
    cellNUMall(t+1,:) = [];
    beadmatchall(t+1,:) = [];
    cellCHLmodeall(t+1,:) = [];
    cellFLSmodeall(t+1,:) = [];
    cellPEmodeall(t+1,:) = [];
    cellSSCmodeall(t+1,:) = [];
    FCBnumberall(t+1,:) = [];
    
    save([procpath 'groupsum'], 'cell*', 'classnotes*', 'beadmatch*', 'FCBnumberall')
    
%     yd_fcb = [cellresultsall(:,1)+5/24];  %UTC
     yd_fcb = [cellresultsall(:,1)];  %UTC, now that abbie4 set to UTC (but with EST zone to keep time stamps from changing)
     
     k=floor(min(yd_fcb))-1;
     sunrisestart = 7.5/24; sunriseend = 7.5/24;  %for real-world data, daylength changes over a data set, which can be approximated by interpolating between these ...
     sunsetstart = sunrisestart+14/24; sunsetend = sunriseend+14/24;
     numdays = max(yd_fcb) - min(yd_fcb);
    sunrise = linspace(sunrisestart, sunriseend, ceil(numdays)+1);
    sunset = linspace(sunsetstart, sunsetend, ceil(numdays)+1);
    clear numdays
    
    ind = find(cellresultsall(:,1));
    
    figure(3)
    clf
    
    subplot(411)
%    sh = plot(yd_fcb(ind), cellNUMall(ind,1)./cellresultsall(ind,2)./cellresultsall(ind,3), '.-'); hold on
    sh = plot(yd_fcb(ind), cellNUMall(ind,1)./cellresultsall(ind,3), '.-'); hold on
    ylabel('\itSyn \rm(cells ml^{-1})')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
%    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
%        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
%        set(f1,'linestyle', 'none')
%    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellNUMall(ind,1)./cellresultsall(ind,3), '.-')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(412)
    sh = plot(yd_fcb(ind), cellNUMall(ind,2)./cellresultsall(ind,3), '.-'); hold on
    ylabel('Cryptos \rm(cells ml^{-1})')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    plot(yd_fcb(ind), cellNUMall(ind,2)./cellresultsall(ind,3), '.-')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(413)
    sh = plot(yd_fcb(ind), cellNUMall(ind,4)./cellresultsall(ind,3), '.-'); hold on
    ylabel('Euks \rm(cells ml^{-1})')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    plot(yd_fcb(ind), cellNUMall(ind,4)./cellresultsall(ind,3), '.-')
    datetick('x', 6, 'keepticks', 'keeplimits')
        
%     subplot(414)
%     sh = plot(yd_fcb(ind), cellFLSall(ind,1), '.-');  hold on
%     ylabel('\itSyn mean FLS')
%     xlabel('Yearday 2003 (local)')
%     t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
%     set(gca, 'layer', 'top')
%     for i=floor(min(yd_fcb)):floor(max(yd_fcb))
%         f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
%         set(f1,'linestyle', 'none')
%     end
%     plot(yd_fcb(ind), cellFLSall(ind,1), '.-')

    figure(4)
    clf
    
    subplot(411)
    sh = plot(yd_fcb(ind), cellSSCall(ind,1), '.-'); hold on
    ylabel('\itSyn \rmSSC mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellSSCall(ind,1), '.-')
    plot(yd_fcb(ind), cellSSCmodeall(ind,1), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(412)
    sh = plot(yd_fcb(ind), cellCHLall(ind,1), '.-'); hold on
    ylabel('\itSyn \rmCHL mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellCHLall(ind,1), '.-')
    plot(yd_fcb(ind), cellCHLmodeall(ind,1), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(413)
    sh = plot(yd_fcb(ind), cellSSCall(ind,4), '.-'); hold on
    ylabel('Euk SSC mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellSSCall(ind,4), '.-')
    plot(yd_fcb(ind), cellSSCmodeall(ind,4), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(414)
    sh = plot(yd_fcb(ind), cellCHLall(ind,4), '.-'); hold on
    ylabel('Euk CHL mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellCHLall(ind,4), '.-')
    plot(yd_fcb(ind), cellCHLmodeall(ind,4), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
   
    figure(9)
    ind = find(yd_fcb < datenum('6-28-03'));
    plot(yd_fcb(ind), cellresultsall(ind,3)./cellresultsall(ind,2)/.05, '.-')
    hold on
    ind = find(yd_fcb > datenum('6-28-03'));
    plot(yd_fcb(ind), cellresultsall(ind,3)./cellresultsall(ind,2)/.025, '.-')
    set(gca, 'ylim', [.90 1.01])
    ylabel('Syringe vol : time est. vol')
    xlabel('Yearday')
    datetick('x', 6, 'keepticks', 'keeplimits')
       
    figure(5)

    subplot(411)
    sh = plot(yd_fcb(ind), cellNUMall(ind,2)./cellresultsall(ind,3), '.-'); hold on
    ylabel('Cryptos \rm(cells ml^{-1})')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    plot(yd_fcb(ind), cellNUMall(ind,2)./cellresultsall(ind,3), '.-')
    datetick('x', 6, 'keepticks', 'keeplimits')    
    title('"Bright" cryptos')
    
    subplot(412)
    sh = plot(yd_fcb(ind), cellSSCall(ind,2), '.-'); hold on
    ylabel('Crypto SSC mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellSSCall(ind,2), '.-')
    plot(yd_fcb(ind), cellSSCmodeall(ind,2), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(413)
    sh = plot(yd_fcb(ind), cellCHLall(ind,2), '.-'); hold on
    ylabel('Crypto CHL mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellCHLall(ind,2), '.-')
    plot(yd_fcb(ind), cellCHLmodeall(ind,2), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')

    figure(6)

    subplot(411)
    sh = plot(yd_fcb(ind), cellNUMall(ind,6)./cellresultsall(ind,3), '.-'); hold on
    ylabel('Cryptos \rm(cells ml^{-1})')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    plot(yd_fcb(ind), cellNUMall(ind,6)./cellresultsall(ind,3), '.-')
    datetick('x', 6, 'keepticks', 'keeplimits')    
    title('"Dim" cryptos')
    
    subplot(412)
    sh = plot(yd_fcb(ind), cellSSCall(ind,6), '.-'); hold on
    ylabel('Crypto SSC mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellSSCall(ind,6), '.-')
    plot(yd_fcb(ind), cellSSCmodeall(ind,6), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    subplot(413)
    sh = plot(yd_fcb(ind), cellCHLall(ind,6), '.-'); hold on
    ylabel('Crypto CHL mean')
    t = axis; minX = t(1); maxX = t(2); minY = t(3); maxY = t(4);
    set(gca, 'layer', 'top')
    for i=floor(min(yd_fcb)):floor(max(yd_fcb))
        f1=fill([i+sunset(i-k)-1; i+sunrise(i-k); i+sunrise(i-k); i+sunset(i-k)-1],[minY minY maxY maxY],[.8 .8 .8]);hold on;
        set(f1,'linestyle', 'none')
    end
    %replot the data so that it is on top of the bars
    plot(yd_fcb(ind), cellCHLall(ind,6), '.-')
    plot(yd_fcb(ind), cellCHLmodeall(ind,6), '.:g')
    datetick('x', 6, 'keepticks', 'keeplimits')
    
    clear ind t sh i k fl 
clear minX minY maxX maxY sunrise* sunset* filelist fileall fileorder