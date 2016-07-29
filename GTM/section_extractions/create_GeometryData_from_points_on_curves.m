% script to read wing section points from Draco extracted data and create
% a GeometryData object.

close all;
clear all;

%%%%%%%%%%%%%%%%%%%%%%%%% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scaleFactor = 1.;
tol = .1;
spanDirection = 2;             % 
chordDirection = 1;                  % 
thicknessDirection = 3;          % 

addpath('/Users/ktrinh/ESAC/SVN/esac/meshing/src/');
%%%%%%%%%%%%%%%%%% inputs: m file from Draco output %%%%%%%%%%%%%%
GTM_wing;
nodes = GTM_wing_data*scaleFactor;
% fine leading edge data

const_stations = const_stations*scaleFactor;

[r_st c] = size(const_stations);
constCoords = nodes(:,spanDirection);

numWingSection = 0;
sections = {};
figure;
grid on;
hold on;
axis equal;
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'.k');

completeSections = 10;     % from post-processing
for i=1:completeSections
    skipList = [];
    if find(skipList==i)>0
        continue;
    end    
    curConstCoord = const_stations(i);
    fprintf('section %d, at %f\n',i,curConstCoord);
    ids = constCoords < (curConstCoord+tol);
    curNodes = nodes(ids,:);
    curConstCoords = curNodes(:,spanDirection);
    ids = find(curConstCoords > (curConstCoord-tol));
    curNodes = curNodes(ids,:);    
    % WingSectionData instantiation
    curWingSection = WingSectionData();
    meanThickPercentage = 90.;
    changeMeanThickList = [14, 16, 18, 20];
    if find(changeMeanThickList==i)>0
        meanThickPercentage = 90.;
    end  
    curWingSection.populateSectionPoints(curNodes,spanDirection,chordDirection,thicknessDirection,meanThickPercentage);
    numWingSection = numWingSection + 1;
    meanU = mean([curWingSection.leadingEdge(spanDirection);...
        curWingSection.trailingEdge(spanDirection)]);
    curWingSection.curves(1).u = meanU;
    curWingSection.curves(2).u = meanU;
    curWingSection.u = meanU;     
    curWingSection.wingDataUpdated = 0;
    sections{i} = curWingSection;
end

numSection = length(sections);


GTM_wing_geom = GeometryData(sections);
save('GTM_wing_geom', 'GTM_wing_geom');
g = GTM_wing_geom;
numSections = length(g.sections);
for i=1:numSections
    plt = 0;
%     checkList = [2;17;20;39];
%     if find(checkList==i)>0
%         plt = 1;
%     end
    g.sections{i}.checkSection(i,plt);
end
g.plot();



