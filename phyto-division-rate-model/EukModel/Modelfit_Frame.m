MR = modelresults; 
allmodelruns = cell(1,2);
allmodelruns{1,1}=modelfits;
allmodelruns{1,2}=allstarts;
allMR = allmodelruns; 


theta = modelresults(2:15); 

if exist('eukvolbins')
    volbins = eukvolbins;
elseif exist('synvolbins')
    volbins = synvolbins;
end


hr1=1;  hr2=25; xsp=0.035; ysp=0.1;
set(gcf,'Position',[119         257        1705         673])
jj = 1; 

[dirsample, simdist,Vt1,Vt2]=simdata_dirichlet_sample(Einterp,N_dist,theta,volbins,hr1,hr2, ts);
        %dirsample=dirsample./repmat(sum(dirsample,1),length(volbins),1);
        %dirsample=[nan(length(volbins),6) dirsample];
        %simdist=[nan(length(volbins),6) simdist];
        
dirsampledist = dirsample ./ sum(dirsample); 
        
        del1=(theta(4).*volbins.^theta(2))./(1+(volbins.^theta(2)));
        y1=theta(1)*ones(size(Einterp));
        ind=find(Einterp < theta(3));
        y1(ind)=(theta(1)/theta(3)) * Einterp(ind);
        
        del2=(theta(8).*volbins.^theta(6))./(1+(volbins.^theta(6)));
        y2=theta(5)*ones(size(Einterp));
        ind=find(Einterp < theta(7));
        y2(ind)=(theta(5)/theta(7)) * Einterp(ind);
 

[mu, ~, ~,] = growth_rate(Einterp,volbins,N_dist,theta,hr1,hr2, ts);


 subplot(2,7,7,'replace')
        subplot_tight(2,7,7,[ysp xsp-0.01])
        title('Parameter Values','fontsize',14)
        text(0.0,0.6,{['gmax 1: ' num2str(theta(1))];...
            ['b 1: ' num2str(theta(2))];
            ['E* 1: ' num2str(theta(3))];
            ['dmax 1: ' num2str(theta(4))];
            ['vol 1: ' num2str(theta(10))];
            ['sigma1: ' num2str(theta(12))];
            ['proportion: ' num2str(theta(9))];
            ['mu1: ' num2str(MR(jj,18))];})
        
        text(0.6,0.6,{['gmax 2: ' num2str(theta(5))];...
            ['b 2: ' num2str(theta(6))];
            ['E* 2: ' num2str(theta(7))];
            ['dmax 2: ' num2str(theta(8))];
            ['vol 2: ' num2str(theta(11))];
            ['sigma2: ' num2str(theta(13))];
            ['proportion: ' num2str(1-theta(9))];
            ['mu2: ' num2str(MR(jj,19))];})

        text(0.0, 0.1,{['s: ' num2str(100*theta(14))];
            ['mu: ' num2str(MR(jj,17))]});
            set(gca,'visible','off')
        
subplot(2,7,1,'replace')
        subplot_tight(2,7,1,[ysp xsp])
        time=0:(1/6):25;
        plot(time,Einterp,'k.')
        xlabel('Time')
        ylabel('Edata')
subplot(2,7,2,'replace')
        subplot_tight(2,7,2,[ysp xsp])
        h0=imagesc(1:25,1:length(volbins),Vhists); set(gca,'Ydir','normal')
        set(h0,'AlphaData',~isnan(Vhists));
        xlabel('Time')
        ylabel('Cell size class')
        title('Observed Data')
        title([datestr(day) ' : ' num2str(day)])
subplot(2,7,3,'replace')
        subplot_tight(2,7,3,[ysp xsp])
        h=imagesc(1:25,1:(length(volbins)),simdist);
        set(h,'AlphaData',~isnan(simdist)), set(gca,'Ydir','normal')
        xlabel('Time')
        ylabel('Cell size class')
        title('MLE Model Fit - no Dir. sample')
subplot(2,7,4,'replace')
        subplot_tight(2,7,4,[ysp xsp])
        h1=imagesc(1:25,1:length(volbins),dirsampledist);
        set(h1,'AlphaData',~isnan(dirsampledist)), set(gca,'Ydir','normal')
        xlabel('Time')
        ylabel('Cell size class')
        title('Sample from Dirichlet')
subplot(2,7,5,'replace'), hold on
        subplot_tight(2,7,5,[ysp xsp]), hold on
        [~,iy]=sort(y1);
        plot(Einterp(iy),y1(iy),'.-','color',[0 0.5 1],'markersize',8);
        [~,iy]=sort(y2);  set(gca,'box','on')
        plot(Einterp(iy),y2(iy),'.-','color',[0 0 0.8],'markersize',8);
        line([max(Einterp) max(Einterp)], ylim,'color',[0.4 0.4 0.4])
        ylabel('Fraction growing cells')
        xlabel('Edata')
subplot(2,7,6,'replace')
        subplot_tight(2,7,6,[ysp xsp]), hold on
        temp=allMR{jj};
        [~, ib]= sort(temp(:,15));
        plot(1:size(temp,1),temp(ib,16),'-','color',[0.6 0.6 0.6]), hold on
        plot(1:size(temp,1),temp(ib,16),'k.')
        xlim([0 size(temp,1)])
        YL=get(gca,'ylim');  ylim([0 YL(2)+0.1]);  set(gca,'box','on')
        line([5 5],[0 YL(2)+0.1],'color',[0.6 0.6 0.6])
        ylabel('Division rate')
        xlabel('# model runs')
subplot(2,7,8,'replace')
        subplot_tight(2,7,8,[ysp xsp])
        h1=plot(1:length(volbins),Vhists(:,hr1:hr1+4),'color',[0.5 0.5 0.5]); hold on
        h2=plot(1:length(volbins),simdist(:,1:5),'color',[0 0 0]);
        h3=plot(1:length(volbins),Vt1(:,1:5),'color',[0 0.5 1]);
        h4=plot(1:length(volbins),Vt2(:,1:5),'color',[0 0 0.8]);
        %plot(1:length(volbins),Vt1(:,1),'color',[1 0.5 0]);
        %plot(1:length(volbins),Vt2(:,1),'color',[1 0 0]);
        title('Hours 7-11')
        legend([h1(1); h2(1); h3(1); h4(1)],'Obs','Model','popn1','popn2','location','NorthWest')
        xlabel('Size class')
        ylabel('Proportion')
        

 subplot(2,7,9,'replace')
        subplot_tight(2,7,9,[ysp xsp])
        plot(1:length(volbins),Vhists(:,hr1+5:hr1+9),'color',[0.5 0.5 0.5]), hold on
        plot(1:length(volbins),simdist(:,6:10),'color',[0 0 0])
        plot(1:length(volbins),Vt1(:,6:10),'color',[0 0.5 1])
        plot(1:length(volbins),Vt2(:,6:10),'color',[0 0 0.8])
        title('Hours 12-16')
        xlabel('Size class')
        ylabel('Proportion')

subplot(2,7,10,'replace')
        subplot_tight(2,7,10,[ysp xsp])
        plot(1:length(volbins),Vhists(:,hr1+10:hr1+14),'color',[0.5 0.5 0.5]), hold on
        plot(1:length(volbins),simdist(:,11:15),'color',[0 0 0])
        plot(1:length(volbins),Vt1(:,11:15),'color',[0 0.5 1])
        plot(1:length(volbins),Vt2(:,11:15),'color',[0 0 0.8])
        title('Hours 17-21')
        xlabel('Size class')
        ylabel('Proportion')
        
 subplot(2,7,11,'replace')
        subplot_tight(2,7,11,[ysp xsp])
        plot(1:length(volbins),Vhists(:,hr1+15:hr1+18),'color',[0.5 0.5 0.5]), hold on
        plot(1:length(volbins),simdist(:,16:19),'color',[0 0 0])
        plot(1:length(volbins),Vt1(:,16:19),'color',[0 0.5 1])
        plot(1:length(volbins),Vt2(:,16:19),'color',[0 0 0.8])
        title('Hours 22-25')
        xlabel('Size class')
        ylabel('Proportion')
        

subplot(2,7,12,'replace'), hold on
        subplot_tight(2,7,12,[ysp xsp]), hold on
        plot(volbins,del1,'.-','color',[0 0.5 1],'markersize',8);
        plot(volbins,del2,'.-','color',[0 0 0.8],'markersize',8);
        ylabel('Fraction dividing cells')
        xlabel('Size bins')
        
      
subplot(2,7,14,'replace')
        subplot_tight(2,7,14,[ysp xsp]), hold on
        plot(MR(:,1),MR(:,17),'.')
        plot(MR(jj,1),MR(jj,17),'rp')
        datetick('x','mm/dd')
        ylabel('Division rate')
        set(gca,'box','on')
        

subplot(2,7,13,'replace')
        subplot_tight(2,7,13,[ysp xsp]), hold on
        temp=allMR{jj};
        [~, ib]= sort(temp(:,15));
        plot(1:size(temp,1),temp(ib,15),'-','color',[0.6 0.6 0.6]), hold on
        plot(1:size(temp,1),temp(ib,15),'.')
        line([5 5],ylim,'color',[0.6 0.6 0.6])
        ylabel('-log L')
        xlabel('# model runs')
        set(gca,'box','on')

       