%% ----------------SECTION 1: DEFINITION OF DIRECTORIES AND COLORS ---------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
close all;
LastCycleID=0;

clc;
workingdir=uigetdir();%select the main folder
cd(workingdir);
addpath(cat(2,workingdir,'\Analysis'));
colormap_directory=cat(2,workingdir,'\Analysis\Colormaps');
addpath(colormap_directory);

analysis_dir=cat(2,workingdir,'\Analysis');
mkdir([analysis_dir,'\PLOTS']); %directory for graphs and plots
plot_dir=cat(2,analysis_dir,'\PLOTS');
kymo_plot_dir=cat(2,workingdir,'\Analysis'); %directory for kymograph-jpgs with overlays of cycles

mu=char(181);

%Determine color maps
%Read all colormaps from the corresponding folder:
Files_in_Colors_Folder=dir(fullfile(""+colormap_directory));
for ii = 1:size(Files_in_Colors_Folder,1)
    A = strfind(Files_in_Colors_Folder(ii).name, 'mat');
      if ~isempty(A)
tmp=load([Files_in_Colors_Folder(ii).name(1:end-4),'.mat']); %in the name, remove the file ending (4 characters: .mat)
color(ii).cmap=tmp.(Files_in_Colors_Folder(ii).name(1:end-4));
color(ii).name=Files_in_Colors_Folder(ii).name(1:end-4);
      end
end

%remove empty entries
temp=[];
for ii=1:numel(color)
    if isempty(color(ii).cmap)==true
        temp=[temp ii];
    end
end
color(temp) = [];

%add colormaps from Matlab:
color(end+1).cmap=colormap(hsv(256));
color(end).name='hsv';
color(end+1).cmap=colormap(jet(256));
color(end).name='jet';
color(end+1).cmap=colormap(autumn(256));
color(end).name='autumn';
close();

  batlow1=[5 31 90]/255;
 batlow2=[55 107 88]/255;
  batlow3=[179 142 47]/255;
  batlow4=[253 171 154]/255;
  batlow5=[251 202 242]/255;
  cyan=[0.06 0.545 0.97];
  yellow=[0.97 0.7 0.06];
  lightblue=[0 0.4470 0.7410];
  bordeaux=[0.6350 0.0780 0.1840];

for m=1:1000000
colorValue(m)=floor(rand*256)+1;
end


%% ----------------SECTION 2: DELETION OF DATA FROM PREVIOUS PROCESSING ----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
% Run only in case if the detection in one of the Kymographs belonging to
% the last multipoint processed in SECTION 3 was faulty. Correct the
% corresponding binary masks and start again in section 3 with the ID of
% the the affected multipoint

%if no data needs to be deleted, continue with SECTION 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delete cycles of last multipoint:
last_multipoint=15; % choose the last multipoint to be deleted
Multipoint(last_multipoint)=[];
temp1=[];
for ii=1:numel(cycle_list)
    if cycle_list(ii).MultipointID==last_multipoint
        temp1=[temp1 ii];
    end
end
cycle_list(min(temp1):LastCycleID) = [];
LastCycleID=min(temp1)-1;
%% ----------------SECTION 3: AUTOMATED TRACKING MODULE --------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The script opens every kymograph of the selected multipoint (with the corresponding binary mask) 
% and detects all Objects (i.e., cells) and links them into cell division cycles(cycle_list struct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
DivisionCriterium=0.88; %if the next cell is e.g. 15% smaller than the previous one, count it as a division event
snakeLengthCriterium=15; %if the cell is filamentous, a higher Division Criterium of at least 0.95 needs to be applied
nextMultipoint=1;
while nextMultipoint==1

% -----------------SELECT THE MULTIPOINT AND IMAGE PARAMETERS--------------
%--------------------------------------------------------------------------
path = uigetdir;
prompt = {'Enter multipoint ID:','trap width [px]:','time step [s]:','Phase Channel:','pixel size [µm]:','size threshold [px]:'};
dlgtitle = 'input: multipoint ID';
dims = [1 40];
definput = {'','30','300','1','0.072','100'};

getfile = inputdlg(prompt,dlgtitle,dims,definput);
multipointID = getfile(1);
trapWidth = getfile(2);
loop=getfile(3);
PhaseContrastChannel=getfile(4);
pixelWidth=getfile(5);
sizeThreshold=getfile(6);

multipointID=cell2mat(multipointID);
multipointID=str2double(multipointID);
trapWidth=cell2mat(trapWidth);
trapWidth=str2double(trapWidth);
loop=cell2mat(loop);
loop=str2double(loop);
PhaseContrastChannel=cell2mat(PhaseContrastChannel);
PhaseContrastChannel=str2double(PhaseContrastChannel);
if PhaseContrastChannel==1
FluorescenceChannel=2;
elseif PhaseContrastChannel==2
FluorescenceChannel=1; 
else
error('Phase contrast channel can be 1 or 2.')
end

pixelWidth=cell2mat(pixelWidth);
pixelWidth=str2double(pixelWidth);
sizeThreshold=cell2mat(sizeThreshold);
sizeThreshold=str2double(sizeThreshold);

%get the tif files and the corresponding masks
Files_in_Multipoint_Folder=dir(fullfile(""+path));

AllFiles = {};
AllMasks = {};
for ii = 1:size(Files_in_Multipoint_Folder,1)
    A = strfind(Files_in_Multipoint_Folder(ii).name, 'tif');
    B = strfind(Files_in_Multipoint_Folder(ii).name, '_Simple Segmentation');
      if ~isempty(A)
        AllFiles = [AllFiles;  Files_in_Multipoint_Folder(ii).name];
      end

    if ~isempty(B)
        AllMasks = [AllMasks;  Files_in_Multipoint_Folder(ii).name];
    end
end


AllTifs = AllFiles(~contains(AllFiles,'_Simple Segmentation')); % All Kymo tifs without the segmentation masks
trapIDs=regexp(AllTifs,'\d*','Match');
filenumber=size(AllMasks,1);

%-----------------READ THE KYMOGRAPHS--------------------------------------
%--------------------------------------------------------------------------
for file_index  = 1:size(trapIDs,1)

trapID=trapIDs{file_index}(end);
Masks=strfind(AllMasks,cat(2,'channel_',cell2mat(trapID),'_'));%find Mask which has the correct trap ID in the name
for ii = 1:size(AllMasks,1)
if ~isempty(cell2mat(Masks(ii)))
    Mask=AllMasks(ii);
end
end
trapID=str2double(cell2mat(trapID));
Mask=cell2mat(Mask);

imageFile=cell2mat(AllTifs(file_index));
imageFile=cat(2,path,'\',imageFile);
MaskFile=cat(2,path,'\',Mask);

PhaseChannel= imread(imageFile,PhaseContrastChannel);
Fluorescence_Ch1=imread(imageFile,FluorescenceChannel);
BW = imread(MaskFile); %read the Binary Mask (data type uint8)

[filepath, name, ext] = fileparts(imageFile);
[height, width, numberOfColorChannels] = size(PhaseChannel);


timepointNo=size(PhaseChannel,2)/trapWidth;
timematrix=zeros(4,timepointNo);
timematrix(1,:)=1:timepointNo; %frames
timematrix(2,:)=timematrix(1,:).*loop; %time in sec
timematrix(3,:)=(timematrix(1,:)-1).*trapWidth; %left boundary of kymo position
timematrix(4,:)=timematrix(1,:).*trapWidth; %right boundary of kymo position


%add zeros to the left and bottom edges in order to exclude the pixels only at
%the top. 
BW=cat(1,BW,zeros(1,size(BW,2)));
BW=cat(2,zeros(size(BW,1),1),BW);
BW=cat(2,BW,zeros(size(BW,1),1));

BW =imclearborder(BW); % exclude objects touching the top border. this command does not change the data type, but changes the BW values to 0 and 1

BW = imerode(BW, strel('disk',1));%erode
BW=imgaussfilt(BW,1.5);%smoothen the image
BW =imdilate(BW, strel('rectangle',[3,3]));%dilate
BW = imfill( BW ,'holes');%fill holes in the binary image

BW = bwareaopen(BW,sizeThreshold); %remove all objects with a smaller area than "sizeThreshold" pixels. This command changes data type to logical
%BW = bwareafilt(BinaryMask,sizeRange); % remove all objects out of sizeRange. changes data type to logical

%remove edges to match the size of the original image
BW= imcrop(BW,[1 0 width-1 height]);


%-----------------DETECT ALL OBJECTS IN THE KYMOGRAPH----------------------
%--------------------------------------------------------------------------
Obj = regionprops(BW,PhaseChannel,{'Centroid','PixelValues','Area','BoundingBox','ConvexHull','ConvexImage','FilledImage','MajorAxisLength'});
ObjFI_1=regionprops(BW,Fluorescence_Ch1,{'Centroid','PixelValues','Area','BoundingBox','ConvexHull','ConvexImage','FilledImage','MajorAxisLength'});

allcentroids = cat(1,Obj.Centroid);
centroidXvalues=allcentroids(:,1);
centroidYvalues=allcentroids(:,2);
allareas= cat(1,Obj.Area);
allboxes=cat(1,Obj.BoundingBox);
allmajAxislengths=cat(1,Obj.MajorAxisLength);
allnumObj=numel(Obj);

% %---------------Adjust brightness of displayed phase contrast image--------
 max_scale=0.7;
 AdjustedImage_ph = imadjust(PhaseChannel,[0 max_scale],[0 1]);
% %--------------------------------------------------------------------------

%----------normalized coordinates of x axis: highlight every 10. Kymo-position
 labelpositions=timematrix(4,:);
 labelpositions=labelpositions(10:10:end)-trapWidth/2;
 labeltexts=string(10:10:numel(labelpositions)*10);


%---------find first and last non-empty Kymo position:
AllPositions=(1:timepointNo)';
xl_Kymo=(1:trapWidth:width)';
xr_Kymo=(trapWidth:trapWidth:width)';
KymoPosBWmean=zeros(timepointNo,1);
for ii=1:timepointNo
KymoPosBWmean(ii)=mean(BW(1:height,xl_Kymo(ii):xr_Kymo(ii)),'all');
end
emptyKymoPositions=find(KymoPosBWmean==0);
emptyKymoPositions=[0;emptyKymoPositions;timepointNo+1];


InitialTime=emptyKymoPositions(diff(emptyKymoPositions)==max(diff(emptyKymoPositions)))+1;
Start=InitialTime;
LastTime=emptyKymoPositions(find(diff(emptyKymoPositions)==max(diff(emptyKymoPositions)))+1)-1;
Finish=LastTime;


%--------------------------------------------------------------------------------------
%--------- AUTOMATED LINEAGE ASSIGNMENT (bottom to top) - with cycle list--------------
%--------------------------------------------------------------------------------------

close all; 
clc;
lineageID=1;
Multipoint(multipointID).Trap(trapID).Lineage = [];

%create positions array of bottom objects
clear centroidsAtTime;
clear BottomPositionsAtTime;

figure;
fig1=gcf;
set(gcf, 'position', [500 2*height width height]);
imshow(AdjustedImage_ph)
ax=gca;
set(gca,'LooseInset',get(gca,'TightInset'));
set(gca, 'visible', 'on')
set(gca,'ytick',[])
set(gca,'FontSize',9)
xticks(labelpositions);
xticklabels(labeltexts)
ax.TickLength = [0.005, 1];
ax.TickDir = 'out';
ax.FontSize=9;
title("multipoint "+multipointID+" - trap "+trapID);
pause(0.25);

for m = Start:Finish
centroidsAtTime(m).YValues=sort(centroidYvalues(centroidXvalues<=timematrix(4,m) & centroidXvalues>timematrix(3,m)));
centroidsAtTime(m).XValues=centroidXvalues(centroidXvalues<=timematrix(4,m) & centroidXvalues>timematrix(3,m));
CellPositionInTrap=numel(centroidsAtTime(m).YValues);
BottomPositionsAtTime(m).YValues=centroidsAtTime(m).YValues(CellPositionInTrap);

end 
  
centroidsAtTimeI=centroidsAtTime;
BottomYValues=zeros(size(centroidsAtTimeI,2),1); %find all Y values of the mother lineage (to flag mother cells)
for ii=1:size(centroidsAtTimeI,2)
if ~isempty(centroidsAtTimeI(ii).YValues)
BottomYValues(ii)=max(centroidsAtTimeI(ii).YValues);
end
end



while isempty(centroidsAtTime(Start).YValues)==0


        while isempty(centroidsAtTime(LastTime).YValues)==1
       LastTime=LastTime-1;  
        end

    for m=1:10000
colorValue1(m)=floor(rand*256)+1;
    end

for m = InitialTime:LastTime

    if isempty(centroidsAtTime(m).YValues)==0
     
    CellPositionInTrap=numel(centroidsAtTime(m).YValues);
    BottomPositionsAtTime(m).YValues=centroidsAtTime(m).YValues(CellPositionInTrap);
    else
    LastTime=m-1;
    break
    end

end

%find the position (ID) of occupied centroids
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs=1; %Initialize DivisionLocs
colorVal=floor(rand*256)+1;
for n=InitialTime:LastTime
XsOfcentroidsAtTime_n=centroidXvalues(centroidXvalues<=timematrix(4,n) & centroidXvalues>timematrix(3,n));
IDsMatchingCentroidYpos=find(centroidYvalues==max(centroidsAtTime(n).YValues)); %IDs matching the object's centroid Y position

for ii=1:numel(IDsMatchingCentroidYpos)
if ismember(Obj(IDsMatchingCentroidYpos(ii)).Centroid(1),XsOfcentroidsAtTime_n)==1
Object.ID(n)=IDsMatchingCentroidYpos(ii);
end
end

%cell stats
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).ObjectID(n)=Object.ID(n);
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n)=Obj(Object.ID(n)).MajorAxisLength*pixelWidth;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellArea(n)=Obj(Object.ID(n)).Area*(pixelWidth)*(pixelWidth);
ObjFI_1(Object.ID(n)).Mean = mean(double(ObjFI_1(Object.ID(n)).PixelValues));
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Intensity_F1(n)=ObjFI_1(Object.ID(n)).Mean;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).XCentroids(n)=Obj(Object.ID(n)).Centroid(1);
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).YCentroids(n)=Obj(Object.ID(n)).Centroid(2);

%-----------------determine the division locations-------------
  if n>1
    
  if Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n)<DivisionCriterium*Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n-1) || (Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n)<0.975*Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n-1) && Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(n-1)>=snakeLengthCriterium)

DivisionLocation=n-1;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs=cat(2,Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs,DivisionLocation);   
  end
  end
  
end

%create centroids-Array without the assigned lineage
for m = InitialTime:LastTime
centroidsAtTime(m).YValues(numel(centroidsAtTime(m).YValues))=[];

end 

Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs(1)=[]; %remove the initialization value
  

DivisionVector=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs;
DivisionVector=cat(2,0,DivisionVector);

%cell stats (2)
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Timespan=InitialTime:LastTime;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Time=(Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Timespan-1)*loop;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).ParentID=0;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).LineageID=lineageID;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).BirthLoc=InitialTime;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID).FinalLoc=LastTime;

%-----------------------HERE BEGINS THE DIVISION OF LINEAGES INTO CYCLES-----------------------------
%----------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------
DivisionLocs=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs;
if isempty(Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs)==1
  cycles=1;  
else
cycles=1:numel(Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs)+1;
end
cycleNo=numel(cycles);

for nn=1:cycleNo
cycle_list(LastCycleID+nn).cycle_id=LastCycleID+nn;
end
LastCycleID=cycle_list(end).cycle_id;


if isempty(Multipoint(multipointID).Trap(trapID).Lineage(lineageID).DivisionLocs)==1
cycle_list(LastCycleID).MultipointID=multipointID;
cycle_list(LastCycleID).trapID=trapID;   
cycle_list(LastCycleID).lineageID=lineageID;
cycle_list(LastCycleID).times= Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Timespan; 
cycle_list(LastCycleID).begin=InitialTime;
cycle_list(LastCycleID).end=LastTime;
cycle_list(LastCycleID).duration=numel(InitialTime:LastTime);
cycle_list(LastCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(InitialTime:LastTime); 
cycle_list(LastCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellArea(InitialTime:LastTime); 
cycle_list(LastCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Intensity_F1(InitialTime:LastTime);
cycle_list(LastCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).XCentroids(InitialTime:LastTime);
cycle_list(LastCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).YCentroids(InitialTime:LastTime);
if min(ismember(cycle_list(LastCycleID).YCentroids,BottomYValues))==1
cycle_list(LastCycleID).motherLineageFlag=1; 
else
cycle_list(LastCycleID).motherLineageFlag=0;     
end
cycle_list(LastCycleID).parentID=[];
cycle_list(LastCycleID).daughterIDs=[];

hold on;
 for n=cycle_list(LastCycleID).begin:cycle_list(LastCycleID).end
 pause(0.001);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue(1+lineageID),:),'LineWidth',1);
 end
 
else

InitialCycleID=LastCycleID-(cycleNo-1);
cycle_list(InitialCycleID).MultipointID=multipointID;
cycle_list(InitialCycleID).trapID=trapID;
cycle_list(InitialCycleID).lineageID=lineageID;
cycle_list(InitialCycleID).times=InitialTime:DivisionLocs(1);
cycle_list(InitialCycleID).begin=InitialTime;
cycle_list(InitialCycleID).end=DivisionLocs(1);
cycle_list(InitialCycleID).duration=numel(cycle_list(InitialCycleID).times);
cycle_list(InitialCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(InitialTime:DivisionLocs(1)); 
cycle_list(InitialCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellArea(InitialTime:DivisionLocs(1)); 
cycle_list(InitialCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Intensity_F1(InitialTime:DivisionLocs(1));
cycle_list(InitialCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).XCentroids(InitialTime:DivisionLocs(1));
cycle_list(InitialCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).YCentroids(InitialTime:DivisionLocs(1));
if min(ismember(cycle_list(InitialCycleID).YCentroids,BottomYValues))==1
cycle_list(InitialCycleID).motherLineageFlag=1; 
else
cycle_list(InitialCycleID).motherLineageFlag=0;     
end
cycle_list(InitialCycleID).parentID=[];
cycle_list(InitialCycleID).daughterIDs=InitialCycleID+1;

cycle_list(LastCycleID).MultipointID=multipointID;
cycle_list(LastCycleID).trapID=trapID;
cycle_list(LastCycleID).lineageID=lineageID;
cycle_list(LastCycleID).times=DivisionLocs(end)+1:LastTime;
cycle_list(LastCycleID).begin=DivisionLocs(end)+1;
cycle_list(LastCycleID).end=LastTime;
cycle_list(LastCycleID).duration=numel(cycle_list(LastCycleID).times);
cycle_list(LastCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(DivisionLocs(end)+1:LastTime); 
cycle_list(LastCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellArea(DivisionLocs(end)+1:LastTime); 
cycle_list(LastCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Intensity_F1(DivisionLocs(end)+1:LastTime);
cycle_list(LastCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).XCentroids(DivisionLocs(end)+1:LastTime);
cycle_list(LastCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).YCentroids(DivisionLocs(end)+1:LastTime);
if min(ismember(cycle_list(LastCycleID).YCentroids,BottomYValues))==1
cycle_list(LastCycleID).motherLineageFlag=1; 
else
cycle_list(LastCycleID).motherLineageFlag=0;     
end
cycle_list(LastCycleID).parentID=LastCycleID-1;

hold on;
 for n=cycle_list(InitialCycleID).begin:cycle_list(InitialCycleID).end
 pause(0.001);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue(colorValue(1+lineageID)),:),'LineWidth',1);
 end
 
for mm=1:cycleNo-2
cycle_list(InitialCycleID+mm).MultipointID=multipointID;
cycle_list(InitialCycleID+mm).trapID=trapID;
cycle_list(InitialCycleID+mm).lineageID=lineageID;
cycle_list(InitialCycleID+mm).times= DivisionLocs(mm)+1:DivisionLocs(mm+1);   
cycle_list(InitialCycleID+mm).begin= DivisionLocs(mm)+1;
cycle_list(InitialCycleID+mm).end= DivisionLocs(mm+1);
cycle_list(InitialCycleID+mm).duration=numel(cycle_list(InitialCycleID+mm).times);
cycle_list(InitialCycleID+mm).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellLength(DivisionLocs(mm)+1:DivisionLocs(mm+1)); 
cycle_list(InitialCycleID+mm).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).cellArea(DivisionLocs(mm)+1:DivisionLocs(mm+1)); 
cycle_list(InitialCycleID+mm).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).Intensity_F1(DivisionLocs(mm)+1:DivisionLocs(mm+1));
cycle_list(InitialCycleID+mm).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).XCentroids(DivisionLocs(mm)+1:DivisionLocs(mm+1));
cycle_list(InitialCycleID+mm).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID).YCentroids(DivisionLocs(mm)+1:DivisionLocs(mm+1));
if min(ismember(cycle_list(InitialCycleID+mm).YCentroids,BottomYValues))==1
cycle_list(InitialCycleID+mm).motherLineageFlag=1; 
else
cycle_list(InitialCycleID+mm).motherLineageFlag=0;     
end
cycle_list(InitialCycleID+mm).parentID=cycle_list(InitialCycleID+mm-1).cycle_id;
cycle_list(InitialCycleID+mm).daughterIDs=cycle_list(InitialCycleID+mm+1).cycle_id;


 for n=cycle_list(InitialCycleID+mm).begin:cycle_list(InitialCycleID+mm).end
 pause(0.002);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue1(mm),:),'LineWidth',1);
 end
 
end

  for n=cycle_list(LastCycleID).begin:cycle_list(LastCycleID).end
 pause(0.002);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue1(cycleNo),:),'LineWidth',1);
 end

end
%----------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------



%-----------------------------------------------
%-----------------------------------------------
%-------------- further lineages:---------------
%-----------------------------------------------
%-----------------------------------------------

while isempty(DivisionVector)==0

%break condition:
ind=[];
for ii=1:numel(centroidsAtTime)
    if isempty(centroidsAtTime(ii).YValues)
        ind=[ind, ii];
    end
end

   if DivisionVector==0
       lineageID=lineageID+1;
        break
   elseif numel(ind)==numel(centroidsAtTime)
       break
   elseif numel(DivisionVector)==2 && DivisionVector(1)==InitialTime-Start && isempty(centroidsAtTime(DivisionVector(2)+1).YValues)==1
        lineageID=lineageID+1;
        break 
   end
   
InitialTimeDaugther=DivisionVector(end)+1;

noDaughter= isempty(centroidsAtTime(InitialTimeDaugther).YValues);
while noDaughter==1
InitialTimeDaugther=DivisionVector(end-1)+1;
DivisionVector(end)=[];

if isempty(centroidsAtTime(InitialTimeDaugther).YValues)==0
InitialTimeDaugther=DivisionVector(end)+1;
noDaughter=0;    
end
end

%LastTimeDaughter corresponds to the last non-empty row of
%centroidsAtTime-struct
clear emptyTrapIndex;
for m = InitialTime:LastTime 
    if isempty(centroidsAtTime(m).YValues)==1
    emptyTrapIndex(m)=1;
    else
    emptyTrapIndex(m)=0;    
    end
end


zeroCounter=[];
     for nn=InitialTimeDaugther:numel(emptyTrapIndex)
         
         if emptyTrapIndex(nn)==0
           zeroCounter=cat(2,zeroCounter,nn);  
           
         end
     end
     
     
     if isempty(find(diff(zeroCounter)~=1))==1
       LastTimeDaugther= max(zeroCounter); 
     else
 LastTimeDaugther=zeroCounter(min(find(diff(zeroCounter)~=1)));   
     end

     if numel(zeroCounter)==1 && zeroCounter>InitialTimeDaugther+1
      LastTimeDaugther=InitialTimeDaugther;
     end
     
       if InitialTimeDaugther==numel(emptyTrapIndex)
         LastTimeDaugther=numel(emptyTrapIndex);
       end


Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).BirthLoc=InitialTimeDaugther;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).FinalLoc=LastTimeDaugther;

%find the bottom positions of the daughter lineage:
  for m = InitialTimeDaugther:LastTimeDaugther 
  BottomPositionsAtTime(m).YValues=max(centroidsAtTime(m).YValues);
  end
  
  %assign the bottom position to the daughter lineage:
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs=1; %Initialize DivisionLocs
colorVal=floor(rand*256)+1;

for n=InitialTimeDaugther:LastTimeDaugther
XsOfcentroidsAtTime_n=centroidXvalues(centroidXvalues<=timematrix(4,n) & centroidXvalues>timematrix(3,n));
IDsMatchingCentroidYpos=find(centroidYvalues==max(centroidsAtTime(n).YValues)); %IDs matching the object's centroid Y position

for ii=1:numel(IDsMatchingCentroidYpos)
if ismember(Obj(IDsMatchingCentroidYpos(ii)).Centroid(1),XsOfcentroidsAtTime_n)==1
Object.ID(n)=IDsMatchingCentroidYpos(ii);
end
end

%cell stats
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).ObjectID(n)=Object.ID(n);
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n)=Obj(Object.ID(n)).MajorAxisLength*pixelWidth;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellArea(n)=Obj(Object.ID(n)).Area*(pixelWidth)*(pixelWidth);
ObjFI_1(Object.ID(n)).Mean = mean(double(ObjFI_1(Object.ID(n)).PixelValues));
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Intensity_F1(n)=ObjFI_1(Object.ID(n)).Mean;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).XCentroids(n)=Obj(Object.ID(n)).Centroid(1);
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).YCentroids(n)=Obj(Object.ID(n)).Centroid(2);

%-----------------determine the division locations-------------
  if n>1
 
  if Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n)<DivisionCriterium*Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n-1) || (Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n)<0.975*Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n-1) && Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(n-1)>=snakeLengthCriterium)

DivisionLocation=n-1;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs=cat(2,Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs,DivisionLocation);   
  end
 
  end
  
end

DivisionVector(end)=[];
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs(1)=[]; %remove the initialization value
DivisionVector=cat(2,DivisionVector,Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs);

%-----------------Identify ParentID-------------
sisterLocation=find(centroidsAtTimeI(InitialTimeDaugther).YValues==Obj(Object.ID(InitialTimeDaugther)).Centroid(2))+1;

for n=1:numel(Obj)
     
   sisterYpos = find(Obj(n).Centroid(:,2)==centroidsAtTimeI(InitialTimeDaugther).YValues(sisterLocation));
    
     if isempty(sisterYpos)==1
     sisterYpos=nan;   
    else
     sisterYpos=n;
     end
 
     sisterID(n)=sisterYpos;

end
 
 sisterID(isnan(sisterID))=[];

      for ii=1:numel(sisterID)
      if ismember(Obj(sisterID(ii)).Centroid(1),centroidsAtTimeI(InitialTimeDaugther).XValues)==1
      sister_ID=sisterID(ii);
      end
     end
     
for n=1:numel(Multipoint(multipointID).Trap(trapID).Lineage)
lineagefinder=find(Multipoint(multipointID).Trap(trapID).Lineage(n).ObjectID==sister_ID);

if isempty(lineagefinder)==1
    lineagefinder=nan;
      else
    lineagefinder=n;
end    
LineageOfSister(n)=lineagefinder;
end
LineageOfSister(isnan(LineageOfSister))=[];

Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).ParentID=LineageOfSister;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).LineageID=lineageID+1;
%-----------------------------------

%create centroids-Array without the assigned lineage

for m = InitialTimeDaugther:LastTimeDaugther
    if isempty(centroidsAtTime(m).YValues)==0 
centroidsAtTime(m).YValues(numel(centroidsAtTime(m).YValues))=[];    
    end
end 


%cell stats (2)
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Timespan=InitialTimeDaugther:LastTimeDaugther;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Time=(Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Timespan-1)*loop;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).BirthLoc=InitialTimeDaugther;
Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).FinalLoc=LastTimeDaugther;



%-----------------------HERE BEGINS THE DIVISION OF LINEAGES INTO CYCLES-----------------------------
%----------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------
DivisionLocs=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs;
if isempty(Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs)==1
  cycles=1;  
else
cycles=1:numel(Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs)+1;
end
cycleNo=numel(cycles);

for nn=1:cycleNo
cycle_list(LastCycleID+nn).cycle_id=LastCycleID+nn;
end
LastCycleID=cycle_list(end).cycle_id;


if isempty(Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).DivisionLocs)==1
cycle_list(LastCycleID).MultipointID=multipointID;
cycle_list(LastCycleID).trapID=trapID;    
cycle_list(LastCycleID).lineageID=lineageID+1;
cycle_list(LastCycleID).times= Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Timespan; 
cycle_list(LastCycleID).begin=InitialTimeDaugther;
cycle_list(LastCycleID).end=LastTimeDaugther;
cycle_list(LastCycleID).duration=numel(InitialTimeDaugther:LastTimeDaugther);
cycle_list(LastCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(InitialTimeDaugther:LastTimeDaugther); 
cycle_list(LastCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellArea(InitialTimeDaugther:LastTimeDaugther); 
cycle_list(LastCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Intensity_F1(InitialTimeDaugther:LastTimeDaugther);
cycle_list(LastCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).XCentroids(InitialTimeDaugther:LastTimeDaugther);
cycle_list(LastCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).YCentroids(InitialTimeDaugther:LastTimeDaugther);
%-----------------Identify ParentID-------------

centroidPos=find(ismember(cycle_list(LastCycleID).XCentroids,centroidsAtTimeI(cycle_list(LastCycleID).begin).XValues));
sisterLoc=find(centroidsAtTimeI(cycle_list(LastCycleID).begin).YValues==cycle_list(LastCycleID).YCentroids(centroidPos))+1;

for nn=1:numel(Obj)
     
   sisterYposition = find(Obj(nn).Centroid(:,2)==centroidsAtTimeI(cycle_list(LastCycleID).begin).YValues(sisterLoc));
    
     if isempty(sisterYposition)==1
     sisterYposition=nan;   
    else
     sisterYposition=nn;
     end
 
     sisterID_(nn)=sisterYposition;

end

 sisterID_(isnan(sisterID_))=[];
 
      for kk=1:numel(sisterID_)
      if ismember(Obj(sisterID_(kk)).Centroid(1),centroidsAtTimeI(cycle_list(LastCycleID).begin).XValues)==1
      sister_ID_=sisterID_(kk);
      end
      end
     
     
for zz=1:numel(cycle_list)
cyclefinder=find(cycle_list(zz).YCentroids==Obj(sister_ID_).Centroid(2) & cycle_list(zz).XCentroids==Obj(sister_ID_).Centroid(1));

if isempty(cyclefinder)==1
    cyclefinder=nan;
      else
    cyclefinder=zz;
end   

CycleOfSister(zz)=cyclefinder;
end
CycleOfSister(isnan(CycleOfSister))=[];
%if several cycles have identical X and Y - centroid coordinates, it
%happens that CycleOfSister gives more than one cycle ID. In this case,
%take simply the last ID, since only the last one can be the sister:
CycleOfSister=CycleOfSister(end);

cycle_list(LastCycleID).parentID=CycleOfSister-1;
%------------------------Add daughterID to the parent cycle------------------------------------
cycle_list(cycle_list(LastCycleID).parentID).daughterIDs=[cycle_list(cycle_list(LastCycleID).parentID).daughterIDs,cycle_list(LastCycleID).cycle_id];
%------------------------Identify ParentID (end)---------------------------

hold on;
 for n=cycle_list(LastCycleID).begin:cycle_list(LastCycleID).end
 pause(0.002);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue(LastCycleID),:),'LineWidth',1);
 end
 
else

InitialCycleID=LastCycleID-(cycleNo-1);
cycle_list(InitialCycleID).MultipointID=multipointID;
cycle_list(InitialCycleID).trapID=trapID;
cycle_list(InitialCycleID).lineageID=lineageID+1;
cycle_list(InitialCycleID).times=InitialTimeDaugther:DivisionLocs(1);
cycle_list(InitialCycleID).begin=InitialTimeDaugther;
cycle_list(InitialCycleID).end=DivisionLocs(1);
cycle_list(InitialCycleID).duration=numel(cycle_list(InitialCycleID).times);
cycle_list(InitialCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(InitialTimeDaugther:DivisionLocs(1)); 
cycle_list(InitialCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellArea(InitialTimeDaugther:DivisionLocs(1)); 
cycle_list(InitialCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Intensity_F1(InitialTimeDaugther:DivisionLocs(1));
cycle_list(InitialCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).XCentroids(InitialTimeDaugther:DivisionLocs(1));
cycle_list(InitialCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).YCentroids(InitialTimeDaugther:DivisionLocs(1));
%-----------------Identify ParentID-------------
sisterLoc=find(centroidsAtTimeI(cycle_list(InitialCycleID).begin).YValues==cycle_list(InitialCycleID).YCentroids)+1;

for nn=1:numel(Obj)
     
   sisterYposition = find(Obj(nn).Centroid(:,2)==centroidsAtTimeI(cycle_list(InitialCycleID).begin).YValues(sisterLoc));
    
     if isempty(sisterYposition)==1
     sisterYposition=nan;   
    else
     sisterYposition=nn;
     end
 
     sisterID_(nn)=sisterYposition;

end
 
 sisterID_(isnan(sisterID_))=[];
 
      for kk=1:numel(sisterID_)
      if ismember(Obj(sisterID_(kk)).Centroid(1),centroidsAtTimeI(cycle_list(InitialCycleID).begin).XValues)==1
      sister_ID_=sisterID_(kk);
      end
      end
     
     
for zz=1:numel(cycle_list)
cyclefinder=find(cycle_list(zz).YCentroids==Obj(sister_ID_).Centroid(2) & cycle_list(zz).XCentroids==Obj(sister_ID_).Centroid(1));

if isempty(cyclefinder)==1
    cyclefinder=nan;
      else
    cyclefinder=zz;
end    
CycleOfSister(zz)=cyclefinder;
end
CycleOfSister(isnan(CycleOfSister))=[];
%if several cycles have identical X and Y - centroid coordinates, it
%happens that CycleOfSister gives more than one cycle ID. In this case,
%take simply the last ID, since only the last one can be the sister:
CycleOfSister=CycleOfSister(end);

if cycle_list(InitialCycleID).begin~=InitialTime
cycle_list(InitialCycleID).parentID=CycleOfSister-1;
%------------------------Add daughterID to the parent cycle------------------------------------
cycle_list(cycle_list(InitialCycleID).parentID).daughterIDs=[cycle_list(cycle_list(InitialCycleID).parentID).daughterIDs,cycle_list(InitialCycleID).cycle_id];
end
%------------------------Identify ParentID (end)-----------------------------------------------

%cycle_list(InitialCycleID).daughterIDs=InitialCycleID+1;
cycle_list(LastCycleID).MultipointID=multipointID;
cycle_list(LastCycleID).trapID=trapID;
cycle_list(LastCycleID).lineageID=lineageID+1;
cycle_list(LastCycleID).times=DivisionLocs(end)+1:LastTimeDaugther;
cycle_list(LastCycleID).begin=DivisionLocs(end)+1;
cycle_list(LastCycleID).end=LastTimeDaugther;
cycle_list(LastCycleID).duration=numel(cycle_list(LastCycleID).times);
cycle_list(LastCycleID).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(DivisionLocs(end)+1:LastTimeDaugther); 
cycle_list(LastCycleID).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellArea(DivisionLocs(end)+1:LastTimeDaugther); 
cycle_list(LastCycleID).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Intensity_F1(DivisionLocs(end)+1:LastTimeDaugther);
cycle_list(LastCycleID).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).XCentroids(DivisionLocs(end)+1:LastTimeDaugther);
cycle_list(LastCycleID).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).YCentroids(DivisionLocs(end)+1:LastTimeDaugther);
cycle_list(LastCycleID).parentID=LastCycleID-1;
cycle_list(cycle_list(LastCycleID).parentID).daughterIDs=[cycle_list(cycle_list(LastCycleID).parentID).daughterIDs,cycle_list(LastCycleID).cycle_id];

hold on;
 for n=cycle_list(InitialCycleID).begin:cycle_list(InitialCycleID).end
 pause(0.001);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue(1+lineageID),:),'LineWidth',1);
 end
 
for mm=1:cycleNo-2
cycle_list(InitialCycleID+mm).MultipointID=multipointID;
cycle_list(InitialCycleID+mm).trapID=trapID;
cycle_list(InitialCycleID+mm).lineageID=lineageID+1;
cycle_list(InitialCycleID+mm).times= DivisionLocs(mm)+1:DivisionLocs(mm+1);   
cycle_list(InitialCycleID+mm).begin= DivisionLocs(mm)+1;
cycle_list(InitialCycleID+mm).end= DivisionLocs(mm+1);
cycle_list(InitialCycleID+mm).duration=numel(cycle_list(InitialCycleID+mm).times);
cycle_list(InitialCycleID+mm).length= Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellLength(DivisionLocs(mm)+1:DivisionLocs(mm+1)); 
cycle_list(InitialCycleID+mm).areas=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).cellArea(DivisionLocs(mm)+1:DivisionLocs(mm+1)); 
cycle_list(InitialCycleID+mm).Intensity_F1=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).Intensity_F1(DivisionLocs(mm)+1:DivisionLocs(mm+1));
cycle_list(InitialCycleID+mm).XCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).XCentroids(DivisionLocs(mm)+1:DivisionLocs(mm+1));
cycle_list(InitialCycleID+mm).YCentroids=Multipoint(multipointID).Trap(trapID).Lineage(lineageID+1).YCentroids(DivisionLocs(mm)+1:DivisionLocs(mm+1));
cycle_list(InitialCycleID+mm).parentID=cycle_list(InitialCycleID+mm-1).cycle_id;
cycle_list(InitialCycleID+mm).daughterIDs=cycle_list(InitialCycleID+mm+1).cycle_id;
%------------------------Add daughterID to the parent cycle------------------------------------
cycle_list(cycle_list(InitialCycleID+mm).parentID).daughterIDs=[cycle_list(cycle_list(InitialCycleID+mm).parentID).daughterIDs,cycle_list(InitialCycleID+mm).cycle_id];


 for n=cycle_list(InitialCycleID+mm).begin:cycle_list(InitialCycleID+mm).end
 pause(0.001);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue1(mm),:),'LineWidth',1);
 end
 
end

  for n=cycle_list(LastCycleID).begin:cycle_list(LastCycleID).end
 pause(0.001);
 plot(Obj(Object.ID(n)).ConvexHull(:,1),Obj(Object.ID(n)).ConvexHull(:,2),'Color',color(32).cmap(colorValue1(cycleNo),:),'LineWidth',1);
 end

end
%----------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------
 
lineageID=lineageID+1;


end

% find the new LastTime
    emptyPositionsCounter=[];
     for nn=InitialTime:numel(emptyTrapIndex)
         
         if emptyTrapIndex(nn)==0
           emptyPositionsCounter=cat(2,emptyPositionsCounter,nn);  
           
         end
     end
     
        if isempty(find(diff(emptyPositionsCounter)~=1))==1

       LastTime= max(emptyPositionsCounter); 
     
     else
 LastTime=emptyPositionsCounter(min(find(diff(emptyPositionsCounter)~=1)));   
        end


%break condition

ind=[];
for ii=1:numel(centroidsAtTime)
    if isempty(centroidsAtTime(ii).YValues)
        ind=[ind, ii];
    end
end

 
if numel(ind)==numel(centroidsAtTime)
   break 
end



end

%remove double daughter assignments
for ii=1:numel(cycle_list)
 daughters= cycle_list(ii).daughterIDs;
 UniqueDaughters=unique(daughters);    
 cycle_list(ii).daughterIDs=UniqueDaughters;
end

%remove empty lineages or cycles:
emptyLineage=0;
for ii=1:numel(Multipoint(multipointID).Trap(trapID).Lineage)
 if  isempty(Multipoint(multipointID).Trap(trapID).Lineage(ii).ObjectID) ==1
   emptyLineage=cat(2,emptyLineage,ii);  
 end
end
emptyLineage(1)=[];

emptyCycle=0;
for ii=1:numel(cycle_list)
 if  cycle_list(ii).duration ==0
   emptyCycle=cat(2,emptyCycle,ii);  
 end
end

emptyCycle(1)=[];

Multipoint(multipointID).Trap(trapID).Lineage(emptyLineage)=[];
cycle_list(emptyCycle)=[];
LastCycleID=cycle_list(end).cycle_id;


 for ii=1:numel(cycle_list)
    if isempty(cycle_list(ii).motherLineageFlag)==1
        cycle_list(ii).motherLineageFlag=0;
    end
 end 

 for ii=1:numel(cycle_list)
   cycle_list(ii).loglength=log(cycle_list(ii).length);
   cycle_list(ii).time_h =(cycle_list(ii).times-1)*loop/3600; %time in h, taken from the time positions
   cycle_list(ii).doublingtime=(cycle_list(ii).end-cycle_list(ii).begin+1)*loop/3600;
   if cycle_list(ii).duration>1
   cycle_list(ii).doublingtimePos=round(mean([cycle_list(ii).end cycle_list(ii).begin]-1)*loop/3600,1);
   else
   cycle_list(ii).doublingtimePos= round((cycle_list(ii).begin-1)*loop/3600,1);
   end
 end

%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
 if isempty(cycle_list(end).cycle_id)==1
  LastCycleID=cycle_list(end-1).cycle_id;
  cycle_list(end) = [];
 end
%%%%%%%%%%%%%%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
exportgraphics(gcf,fullfile(kymo_plot_dir,"Lineages_XY"+multipointID+"_Trap"+trapID+".jpg"),'Resolution',300,'BackgroundColor','white')
close;
%%%%%%%%%%%%%%%%%$$$$$$$$$$4
save("CycleData.mat","cycle_list","Multipoint","LastCycleID","loop","trapWidth","pixelWidth","timepointNo","timematrix","height","width");
%%%%%%%%%%%%%%%%%$$$$$$$$$$4
end
%------------------------------------------------------

%------------------------------------------------------
ContinueMultipointSelectionDlg = questdlg("continue selecting the next multipoint? Last Multipoint ID: "+multipointID+"");
switch ContinueMultipointSelectionDlg
    case 'Yes'
   nextMultipoint=1;

    case 'No'
nextMultipoint=0;
end

end

%%%%%%%%%%%%%%%%%$$$$$$$$$$4
save("CycleData.mat","cycle_list","Multipoint","LastCycleID","loop","trapWidth","pixelWidth","timepointNo","timematrix","height","width");
%%%%%%%%%%%%%%%%%$$$$$$$$$$4
%% CYCLE DATA - TOTAL
clc;
%-----------------growth rate interval calculation method----------------------------------------
%Determine the growth rate in intervals of a size of max. Int_N positions

minimalCellNo=3; %minimal no of cells to calculate a growth rate
Int_Pos=4; %preferrable amount of cells in one Interval to calculate the growth rate
Interval_PositionCounter=1:1000;

%minimalCellNo=3; %minimal no of cells to calculate a growth rate
%Int_Pos=5; % e.g. 5 positions correspond to 4*loop=1200s with loop=300s. I.e., the growth rate is calculated over a time of 1200s
% between 3 and 5   (+2) positions: one Interval (6-7 would be mean that first one is 5 and the second interval would be too short)
% between 8 and 10  (+2) positions: two intervals(5 and the second one up to 7)
% between 13 and 15 (+2) positions: three intervals (5,5,3..7)
% between 18 and 20 (+2) positions: four intervals (5,5,5,3..7)
%...
%Interval_LowerBoundaries=Int_Pos*(Interval_PositionCounter-1)+3;
%Interval_UpperBoundaries=Int_Pos*Interval_PositionCounter+2;

%minimalCellNo=8; %minimal no of cells to calculate a growth rate
%Int_Pos=12; %preferrable amount of cells in one Interval to calculate the growth rate
%between 8 and 12  (+7) positions: one Interval
%between 20 and 24 (+7) positions: two Intervals (12 and the second one up to 19)
%between 32 and 36 (+7) positions: three Intervals (12,12,8..19)


Interval_LowerBoundaries=Int_Pos*(Interval_PositionCounter-1)+minimalCellNo;
Interval_UpperBoundaries=Int_Pos*Interval_PositionCounter+(minimalCellNo-1);

%--------------------------------------------------------------------------------
h = waitbar(0,'1','Name','calculating growth rates',...
    'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
setappdata(h,'canceling',0);

for ii=1:numel(cycle_list)
    if getappdata(h,'canceling')
        break
    end
   % waitbar((ii-1)/numel(cycle_list),h,sprintf('cycle %d',ii))
 
   if numel(cycle_list(ii).time_h)>=minimalCellNo
   %number of intervals within cycle
   NoOfIntervals=max(find(Interval_LowerBoundaries<=numel(cycle_list(ii).times)));
   
   if NoOfIntervals==1
    f=fit(cycle_list(ii).time_h.',cycle_list(ii).loglength.','poly1');
    cycle_list(ii).growthrate=f.p1;
    cycle_list(ii).timeCenter = mean(cycle_list(ii).time_h);
    cycle_list(ii).cycleMeanFluorescence1=mean(cycle_list(ii).Intensity_F1);
    cycle_list(ii).cycleMaxFluorescence1=max(cycle_list(ii).Intensity_F1);
    cycle_list(ii).cycleDeltaFluorescence1=cycle_list(ii).Intensity_F1(end)-cycle_list(ii).Intensity_F1(1);
    cycle_list(ii).Location=mean(cycle_list(ii).YCentroids);

   elseif NoOfIntervals>1

for m=1:NoOfIntervals-1
clear x;
clear y; 
timeIntvlStart=cycle_list(ii).times(Int_Pos*(m-1)+1);
timeIntvlEnd=cycle_list(ii).times(Int_Pos*m);
y=cycle_list(ii).loglength(Int_Pos*(m-1)+1:Int_Pos*m);
x=((1:numel(y))-1)*loop/3600;

f=fit(x.',y.','poly1');
cycle_list(ii).growthrate(m)=f.p1;
Intensities1_in_Interval=cycle_list(ii).Intensity_F1(find(cycle_list(ii).times==timeIntvlStart):find(cycle_list(ii).times==timeIntvlEnd));

cycle_list(ii).timeCenter(m)=mean(([timeIntvlStart timeIntvlEnd]-1)*loop/3600);
cycle_list(ii).cycleMeanFluorescence1(m)=mean(Intensities1_in_Interval);
cycle_list(ii).cycleMaxFluorescence1(m)=max(Intensities1_in_Interval);
cycle_list(ii).cycleDeltaFluorescence1(m)=Intensities1_in_Interval(end)-Intensities1_in_Interval(1);
cycle_list(ii).Location(m)=mean(cycle_list(ii).YCentroids(find(cycle_list(ii).times==timeIntvlStart):find(cycle_list(ii).times==timeIntvlEnd)));
end
%consider the last interval:
clear x;
clear y; 
timeIntvlStart=cycle_list(ii).times(Int_Pos*(NoOfIntervals-1)+1);
timeIntvlEnd=cycle_list(ii).times(end);
Intensities1_in_Interval=cycle_list(ii).Intensity_F1(find(cycle_list(ii).times==timeIntvlStart):find(cycle_list(ii).times==timeIntvlEnd));
y=cycle_list(ii).loglength(Int_Pos*(NoOfIntervals-1)+1:numel(cycle_list(ii).loglength));
x=((1:numel(y))-1)*loop/3600;

f=fit(x.',y.','poly1');
cycle_list(ii).growthrate(NoOfIntervals)=f.p1;
cycle_list(ii).timeCenter(NoOfIntervals)=mean(([timeIntvlStart timeIntvlEnd]-1)*loop/3600);
cycle_list(ii).cycleMeanFluorescence1(NoOfIntervals)=mean(Intensities1_in_Interval);
cycle_list(ii).cycleMaxFluorescence1(NoOfIntervals)=max(Intensities1_in_Interval);
cycle_list(ii).cycleDeltaFluorescence1(NoOfIntervals)=Intensities1_in_Interval(end)-Intensities1_in_Interval(1);
cycle_list(ii).Location(NoOfIntervals)=mean(cycle_list(ii).YCentroids(find(cycle_list(ii).times==timeIntvlStart):find(cycle_list(ii).times==timeIntvlEnd)));

   end

   end

  waitbar(ii/numel(cycle_list),h,sprintf('cycle %d of %d',ii,numel(cycle_list)))

end
delete(h)

save("CycleData.mat","cycle_list","Multipoint","LastCycleID","loop","trapWidth","pixelWidth","timepointNo","timematrix","height","width","minimalCellNo","Int_Pos","NoOfIntervals");
