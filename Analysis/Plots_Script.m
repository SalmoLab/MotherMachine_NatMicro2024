%% ----------------SECTION 1: Load Cycle Data  -----------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
% load the CycleData files from both replicates and combine into one struct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load CycleData_20231017_4h_R1
load CycleData_20231025_4h_R2

loop=300; %5 min loop time
%switch to spent medium occurs at frame 35->36 for R1
%switch to spent medium occurs at frame 36->37 for R2

% add +1 to all cycle times of R1, or 1x300/3600=
for ii=1:numel(cycle_list_R1)
cycle_list_R1(ii).times=cycle_list_R1(ii).times+1;
cycle_list_R1(ii).begin=cycle_list_R1(ii).begin+1;
cycle_list_R1(ii).end=cycle_list_R1(ii).end+1;
cycle_list_R1(ii).time_h=cycle_list_R1(ii).time_h+(1*loop)/3600;
cycle_list_R1(ii).timeCenter=cycle_list_R1(ii).timeCenter+(1*loop)/3600;
end

%change the range of the multipoints for R2, starting from 17 and adjust
%the cycle IDs
for ii=1:numel(cycle_list_R2)
cycle_list_R2(ii).MultipointID=cycle_list_R2(ii).MultipointID+16;
cycle_list_R2(ii).cycle_id=cycle_list_R2(ii).cycle_id+cycle_list_R1(end).cycle_id;
end


combined_cycle_list = [cycle_list_R1,cycle_list_R2];

%% ----------------SECTION 2: Sort data according to Multipoints  ----------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
% Sort the data into the categories toxin+ and toxin-
%Exclude potential artifacts, such as cycles located very close to exit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear timeCenters;
 clear timeCenterStep;

%choose subset of multipoints
condition(1).Multipoints=[1,2,3,4,5,6,7,8,17,18,19,20,21,22,23];
condition(1).title='toxin-';
condition(2).Multipoints=[9,10,11,12,13,14,15,16,24,25,26,27,28,29,30];
condition(2).title='toxin+';

timeCenterStep=0.1;
growthrateThreshold=2;
mingrowthrateThreshold=-1;
exposureStart=2.8;
TimeOfExposure=4;
recoveryStart=8; %estimated recovery start
close_to_exit_threshold=50; %distance in pixels from MM trap exit into the trap, in which the growth rate will not be calculated


for k=1:numel(condition)
condition(k).growthrates=[];
condition(k).growthrateTimeCenters=[];
condition(k).CycleLocation=[];
condition(k).CycleFinalLength=[];
condition(k).CycleEnd=[];
condition(k).CycleDuration=[];
end

%concatenate growth rates, time centers and fluorescence according to multipoint. 
for k=1:numel(condition)
for ii=1:numel(combined_cycle_list)
         if ismember(combined_cycle_list(ii).MultipointID,condition(k).Multipoints)==1

   if isempty(combined_cycle_list(ii).growthrate)==1
          condition(k).growthrates=[condition(k).growthrates 0];
          condition(k).growthrateTimeCenters=[condition(k).growthrateTimeCenters 0];
          condition(k).CycleLocation=[condition(k).CycleLocation 0];

   else
        condition(k).growthrates=[condition(k).growthrates combined_cycle_list(ii).growthrate];
        condition(k).growthrateTimeCenters=[condition(k).growthrateTimeCenters combined_cycle_list(ii).timeCenter];
        condition(k).CycleLocation=[condition(k).CycleLocation combined_cycle_list(ii).Location];
   end


   % concatenate cycle final lengths. 
    if isempty(combined_cycle_list(ii).duration)==1

      condition(k).CycleFinalLength(ii)=[condition(k).CycleFinalLength 0];
      condition(k).CycleEnd(ii)=[condition(k).CycleEnd 0];
      condition(k).CycleDuration(ii)=[condition(k).CycleDuration 0];
     else
      condition(k).CycleFinalLength=[condition(k).CycleFinalLength combined_cycle_list(ii).length(end)];
      condition(k).CycleEnd=[condition(k).CycleEnd combined_cycle_list(ii).end];
      condition(k).CycleDuration=[condition(k).CycleDuration combined_cycle_list(ii).duration];
    end 


         end
end

condition(k).growthrates=condition(k).growthrates.';
condition(k).growthrateTimeCenters=condition(k).growthrateTimeCenters.';
condition(k).CycleLocation=condition(k).CycleLocation.';
condition(k).CycleFinalLength=condition(k).CycleFinalLength.';
condition(k).CycleEnd=condition(k).CycleEnd.';
condition(k).CycleDuration=condition(k).CycleDuration.';

%filter growthrates
condition(k).growthrates(condition(k).growthrates>growthrateThreshold)=0;%filter out unprobable growthrate above a threshold by setting those values to zero
condition(k).growthrates(condition(k).growthrates<mingrowthrateThreshold)=0;%set growth rates below a threshold (especially negative values) to zero
condition(k).growthrates(condition(k).CycleLocation<close_to_exit_threshold)=0; %filter out growth rates close to MM trap exit to avoid artifacts


end


timeCenters=min(round(condition(1).growthrateTimeCenters,1)):timeCenterStep:max(round(condition(1).growthrateTimeCenters,1));
timeCenters=cat(2,timeCenters,timeCenters(end)+timeCenterStep);


 %remove empty entries
 for k=1:numel(condition)
 for ii=1:numel(condition(k).growthrates)
     if isempty(condition(k).growthrates(ii))==1
         condition(k).growthrates(ii)=0;
     end

 end
emptyEntries=find(condition(k).growthrates==0);
condition(k).growthrates(emptyEntries)=[];
condition(k).growthrateTimeCenters(emptyEntries)=[];
condition(k).CycleLocation(emptyEntries)=[];
 end

%calculate the mean growth rate and fluorescence in defined time steps
for k=1:numel(condition)

clear growthrates_in_Interval;
for ii=1:numel(timeCenters)-1
growthrates_in_Interval=condition(k).growthrates(find(timeCenters(ii)<=condition(k).growthrateTimeCenters & condition(k).growthrateTimeCenters<timeCenters(ii+1)));

condition(k).meanGrowthRate(ii)=mean(growthrates_in_Interval,'omitnan');
condition(k).stdGrowthRate(ii)=std(growthrates_in_Interval,'omitnan');

end
end

%% ----------------SECTION 3: Plot the mean growth rate  -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
clc;
close all;
exposureStart=2.8;
TimeOfExposure=4;
exposureFinish=exposureStart+TimeOfExposure;

light_blue1=[107/255 139/255 185/255]; %light blue
light_blue2=[101/255 149/255 213/255]; 
light_blue3=[140/255 184/255 230/255]; 

dark_blue1=[46/255 71/255 82/255]; %dark blue
dark_blue2=[101/255 149/255 213/255]; 

dark_gray1=[115/255 115/255 115/255]; %dark gray
light_gray1=[146/255 148/255 152/255]; %light gray
light_gray2=[205/255 205/255 205/255]; %light gray

light_green1=[183/255 214/255 194/255]; %light green

color_toxin_pos=light_gray1;
color_toxin_neg=light_blue3;

k=2; %data for toxin +
fig=figure;
fig.Position = [200 50 1000 450];

errorbar(timeCenters(1:end-1), condition(k).meanGrowthRate(1:end),condition(k).stdGrowthRate(1:end),"-o","LineWidth",0.8,"Color",color_toxin_pos,"MarkerSize",4.5,"MarkerEdgeColor","k","MarkerFaceColor",color_toxin_pos)

%add the data for toxin -
k=1;
hold on
errorbar(timeCenters(1:end-1), condition(k).meanGrowthRate(1:end),condition(k).stdGrowthRate(1:end),"-o","LineWidth",0.8,"Color",color_toxin_neg,"MarkerSize",4.5,"MarkerEdgeColor","k","MarkerFaceColor",color_toxin_neg)

lgd=legend('SN MK01 WT','SN MK01 sKO','FontSize',18,'Location','northeast','Orientation','horizontal','AutoUpdate','off','fontweight','normal');
lgd.NumColumns  = 2;

ax = gca; %get current axis
ax.FontSize=18; %font size of axis (including labels)
ax.XLim=[0 15]; %axis limits
ax.YLim =[-0.4 2.5];

%set ticks
x_ticks=0:2:14;
x_ticklabels=num2str(xticks);
set(gca,'xtick',x_ticks)

hold on
xline(exposureStart,'k-',{sprintf("")},'FontSize',10,'LineStyle','--','LineWidth',0.75,'LabelOrientation','horizontal');
xline(exposureFinish,'k-',{sprintf("")},'FontSize',10,'LineStyle','--','LineWidth',0.75,'LabelOrientation','horizontal');
hold on
patch([exposureStart exposureFinish exposureFinish exposureStart], [ax.YLim(1) ax.YLim(1), ax.YLim(2) ax.YLim(2)],[0 0 0],'FaceAlpha',0.1,'EdgeColor','none')

NE = [max(xlim) max(ylim)]-[diff(xlim) diff(ylim)]*0.05;
%text(NE(1)*0.9, NE(2),"n= "+numel(condition(k).growthrates)+"",'Interpreter', 'latex','Color','black','FontSize',12, 'VerticalAlignment','top', 'HorizontalAlignment','left');
set(gca,'fontweight','bold')
xlabel("time (h)",'FontSize',20);
ylabel("{\it S}. Tm growth rate (1/h)",'FontSize',18);

set(gca,'XMinorTick','on','YMinorTick','off')
set(gca,'TickDir','out');
set(gca,'box','off')

grid on;
   hold on
 ax = axis;

 plot(ax(2)*[1,1],ax(3:4),'k','linewidth',0.4)
 plot(ax(1:2),ax(4)*[1,1],'k','linewidth',0.4)
 
%title("growth rate vs time: toxin +",'FontSize',11);
exportgraphics(gcf,fullfile(plot_dir,"R1 and R2 - growth rate vs time.jpg"),'Resolution',900,'BackgroundColor','white')
exportgraphics(gcf,fullfile(plot_dir,"R1 and R2 - growth rate vs time.pdf"),'BackgroundColor','white')

%% ----------------SECTION 4: percentage of lysis, filamentous growth, growth arrest and regrowth  -------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%create categories of mother lineage fates and plot the corresponding
%percentages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
close all;

condition(1).toxin=0; %condition 1 is toxin negative
condition(2).toxin=1; %condition 2 is toxin psoitive

%colors:
lysisColor=[255 136 140]/255;
lysis_after_filamentationColor=[255 170 170]/255;
growthArrestColor=[75 110 170]/255;
growthArrestafter_filamentationColor=[107 139 185]/255;
regrowthColor=[160 210 170]/255;
regrowthafter_filamentationColor=[183 214 194]/255;
FilamentationColor=[200 200 200]/255;

replicateNo=2;
for k=1:numel(condition)
for replicateID=1:3 %1: R1 and R2 together. 2:R1. 3:R2
  
    clear lysis_and_regrowth_data
    clear conditionFlag
    clear filamentation
    clear lysis
    clear growth_arrest
    clear regrowth

lysis_and_regrowth_data = xlsread(cat(2,cd,'/lysis events R1 and R2.xlsx'),replicateID); %sheet 1 corresponds to all events from R1 (2023-10-17) and R2 (2023-10-25)
conditionFlag=lysis_and_regrowth_data(:,1); %1:toxin+, 0:toxin-
filamentation=lysis_and_regrowth_data(:,4);
lysis=lysis_and_regrowth_data(:,5);
growth_arrest=lysis_and_regrowth_data(:,6);
regrowth=lysis_and_regrowth_data(:,7);

condition(k).Replicate(replicateID).mother_no=numel(conditionFlag(conditionFlag==condition(k).toxin));
condition(k).Replicate(replicateID).Lysis_total=sum(lysis(conditionFlag==condition(k).toxin),'omitnan');
condition(k).Replicate(replicateID).Filamentous_total=sum(filamentation(conditionFlag==condition(k).toxin),'omitnan');
condition(k).Replicate(replicateID).growth_arrest_total=sum(growth_arrest(conditionFlag==condition(k).toxin),'omitnan');
condition(k).Replicate(replicateID).regrowth_total=sum(regrowth(conditionFlag==condition(k).toxin),'omitnan');

%categories:
condition(k).Replicate(replicateID).Lysis_after_filamenteous=sum(lysis(conditionFlag==condition(k).toxin & filamentation==1),'omitnan');
condition(k).Replicate(replicateID).Lysis_non_filamenteous=(condition(k).Replicate(replicateID).Lysis_total-condition(k).Replicate(replicateID).Lysis_after_filamenteous);
condition(k).Replicate(replicateID).growth_arrest_after_filamenteous=sum(growth_arrest(conditionFlag==condition(k).toxin & filamentation==1),'omitnan');
condition(k).Replicate(replicateID).growth_arrest_non_filamenteous=(condition(k).Replicate(replicateID).growth_arrest_total-condition(k).Replicate(replicateID).growth_arrest_after_filamenteous);
condition(k).Replicate(replicateID).regrowth_after_filamenteous=sum(regrowth(conditionFlag==condition(k).toxin & filamentation==1),'omitnan');
condition(k).Replicate(replicateID).regrowth_non_filamenteous=(condition(k).Replicate(replicateID).regrowth_total-condition(k).Replicate(replicateID).regrowth_after_filamenteous);
condition(k).Replicate(replicateID).JustFilamentation=(condition(k).Replicate(replicateID).Filamentous_total-condition(k).Replicate(replicateID).Lysis_after_filamenteous-condition(k).Replicate(replicateID).growth_arrest_after_filamenteous-condition(k).Replicate(replicateID).regrowth_after_filamenteous);

end
end

for k=1:numel(condition)

condition(k).RegrowthMatrix(:,1)=[condition(k).Replicate.Lysis_non_filamenteous].';
condition(k).RegrowthMatrix(:,2)=[condition(k).Replicate.Lysis_after_filamenteous].';
condition(k).RegrowthMatrix(:,3)=[condition(k).Replicate.growth_arrest_non_filamenteous].';
condition(k).RegrowthMatrix(:,4)=[condition(k).Replicate.growth_arrest_after_filamenteous].';
condition(k).RegrowthMatrix(:,5)=[condition(k).Replicate.regrowth_non_filamenteous].';
condition(k).RegrowthMatrix(:,6)=[condition(k).Replicate.regrowth_after_filamenteous].';
condition(k).RegrowthMatrix(:,7)=[condition(k).Replicate.JustFilamentation].';


%reduced matrix with only 4 catergories:
%lysis
%growth arrest
%growth
%filamentous escape
condition(k).RegrowthMatrix_reduced(:,1)=[condition(k).Replicate.Lysis_after_filamenteous].'+[condition(k).Replicate.Lysis_non_filamenteous].';
condition(k).RegrowthMatrix_reduced(:,2)=[condition(k).Replicate.growth_arrest_after_filamenteous].'+[condition(k).Replicate.growth_arrest_non_filamenteous].';
condition(k).RegrowthMatrix_reduced(:,3)=[condition(k).Replicate.regrowth_after_filamenteous].'+[condition(k).Replicate.regrowth_non_filamenteous].';
condition(k).RegrowthMatrix_reduced(:,4)=[condition(k).Replicate.JustFilamentation].';
end


for k=1:numel(condition)
for replicateID=1:3
condition(k).RegrowthMatrix(replicateID,:)=condition(k).RegrowthMatrix(replicateID,:)./condition(k).Replicate(replicateID).mother_no;  
condition(k).RegrowthMatrix_reduced(replicateID,:)=condition(k).RegrowthMatrix_reduced(replicateID,:)./condition(k).Replicate(replicateID).mother_no;
end
end

%calculate errors between all replicates
for k=1:numel(condition)
condition(k).RegrowthMatrix_reduced_Errors=std(condition(k).RegrowthMatrix_reduced(2:replicateNo+1,:));
end

%-------------------Bar Plots-------------------------------
%-----------------------------------------------------------
clc
close all

light_blue3=[140/255 184/255 230/255]; 
light_gray1=[146/255 148/255 152/255]; %light gray


%light blue - grey
colors=[light_gray1;light_blue3];

    clear err
    err=[condition(1).RegrowthMatrix_reduced_Errors;condition(2).RegrowthMatrix_reduced_Errors].';

%for replicateID=1:3
for replicateID=1 %replicateID 1 corresponds to R1 and R2 combined

    clear b
    clear temp
    clear y
    temp=replicateID-1;
fig=figure;
fig.Position = [500 150 1000 400];
y=[condition(1).RegrowthMatrix_reduced(replicateID,:);condition(2).RegrowthMatrix_reduced(replicateID,:)].';
%if the repllicateID-1 is zero, i.e. all replicates together, then add errorbars to the bar plot
if replicateID==1
y_replicates=[condition(1).RegrowthMatrix_reduced(replicateID+1:end,:);condition(2).RegrowthMatrix_reduced(replicateID+1:end,:)].';

b= bar(y,'LineWidth',0.5,'FaceColor','flat','FaceAlpha',1,'EdgeColor','flat','LineWidth',1);

hold on
for ii = 1:size(y,2)
    % get x positions per group
    xpos = b(ii).XData + b(ii).XOffset;
    % draw errorbar
    errorbar(xpos, y(:,ii), err(:,ii), 'LineStyle', 'none','Color', 'k', 'LineWidth', 1);
    if ii==2
    scatter(xpos+rand*0.06*(-1),y_replicates(:,3),60,'^','MarkerFaceColor','white','MarkerEdgeColor','k','LineWidth',1)
    scatter(xpos+rand*0.06,y_replicates(:,4),60,'^','MarkerFaceColor','white','MarkerEdgeColor','k','LineWidth',1)
    else
    scatter(xpos+rand*0.06*(-1),y_replicates(:,1),60,'v','MarkerFaceColor','white','MarkerEdgeColor','k','LineWidth',1)    
    scatter(xpos+rand*0.06,y_replicates(:,2),60,'v','MarkerFaceColor','white','MarkerEdgeColor','k','LineWidth',1)    
    end
end
bar(y,'LineWidth',0.5,'FaceColor','flat','FaceAlpha',0,'EdgeColor','k','LineWidth',0.8);
hold off

else
b=bar(y,'LineWidth',0.5,'FaceColor','flat','EdgeColor','flat','LineWidth',1);    
hold on
bar(y,'LineWidth',0.5,'FaceColor','flat','FaceAlpha',0,'EdgeColor','k','LineWidth',0.8);
end

 b(1).CData(1,:) = colors(2,:);
 b(1).CData(2,:) = colors(2,:);
 b(1).CData(3,:) = colors(2,:);
 b(1).CData(4,:) = colors(2,:);
 b(2).CData(1,:) = colors(1,:);
 b(2).CData(2,:) = colors(1,:);
 b(2).CData(3,:) = colors(1,:);
 b(2).CData(4,:) = colors(1,:);

 ax = gca; %get current axis
 ax.FontSize=18; %font size of axis (including labels)
 ax.YLim =[0 1.1];
 set(gca,'xticklabel',{'lysis','growth arrest','regrowth','filamentation'});
 set(gca,"fontweight","bold");
 ylabel("relative amount",'FontSize',20);

lgd=legend('SN MK01 sKO','SN MK01 WT','FontSize',18,'Location','northeast','AutoUpdate','off',"fontweight","normal");
lgd.NumColumns  = 1;

line([1.5, 1.5], ylim,'LineStyle','--','Color', 'k', 'LineWidth', 1.5);
line([2.5, 2.5], ylim,'LineStyle','--','Color', 'k', 'LineWidth', 1.5);
line([3.5, 3.5], ylim,'LineStyle','--','Color', 'k', 'LineWidth', 1.5);

grid on
ax.XAxis.TickLength = [0 0];
ax.XGrid = 'off';
exportgraphics(gcf,fullfile(plot_dir,"R"+temp+" - Percentage of lysis regrowth filamentation 4h-no errorbar.jpg"),'Resolution',900,'BackgroundColor','white');
exportgraphics(gcf,fullfile(plot_dir,"R"+temp+" - Percentage of lysis regrowth filamentation 4h-no errorbar.pdf"),'BackgroundColor','white');

end