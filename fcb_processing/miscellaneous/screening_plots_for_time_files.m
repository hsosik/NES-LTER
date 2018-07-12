
year2do='2009';

timelist=dir(['\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan' year2do '\data\processed\time\*.mat']);
pathname=['\\sosiknas1\Lab_data\MVCO\FCB\MVCO_Jan' year2do '\data\processed\time\'];

%%
for j=1:length(timelist)
   load([pathname timelist(j).name]) 
   tempvar=whos('FCB*'); tempvar=tempvar.name;
   eval(['temptime=' tempvar ';']) 
   
   %find all the cell or bead syringes:
   ii=find(temptime(:,6) == 3 | temptime(:,6) == 6); %flag is either cells or beads
   
   subplot(2,1,1,'replace'), hold on
   plot(temptime(:,2),temptime(:,4),'.-')
   plot(temptime(ii,2),temptime(ii,4),'.')
   datetick('x','mm/dd')
   
   ylabel('Acquisition time (s)')
   title(['file ' num2str(j) ' out of ' num2str(length(timelist)) ', ' tempvar],'Interpreter','none')
   subplot(2,1,2,'replace'), hold on
   [tt]=histc(temptime(ii,4),[1:100]);
   bar([1:100],tt);
   xlabel('Acquisition time (s)')
   
   keyboard
   eval(['clearvars ' tempvar]) %clear this variable

end
