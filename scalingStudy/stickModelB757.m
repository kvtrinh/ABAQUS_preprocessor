% script to generate stick model for lattice & shell B757 wing

close all;
clear all;

load('C:\Users\ktrinh\MADCAT\scaling_study\gtm_stiff_wing_0_percent_fuel_data.mat');
Lambda = 22.896;

% convert to inch
x = xd*cosd(Lambda)*12.;
y = xd*sind(Lambda)*12.;
z = xd*atand(Lambda)*12.;
numNodes = length(x);
numElems = numNodes-1;

figure;
hold on;
axis equal;
grid on;
plot3(x,y,z,'*r');
plot3(x,y,z,'b');

outfile = 'VoxelShell_GTM_stick.inp';
fid=fopen(outfile,'w');
fprintf(fid,'*Heading\n');
fprintf(fid,'Stick model for GTM wing with voxel substructure\n');
fprintf(fid,'*Part, name=GTM\n');
fprintf(fid,'**\n');
fprintf(fid,'*Node\n');
for i=1:numNodes
    fprintf(fid,'%d, %f, %f, %f\n',i,x(i), y(i), z(i));
    
end
fprintf(fid,'*ELEMENT, TYPE=B31, ELSET=Beams\n');
for i=1:numElems
    fprintf(fid,'%d, %d, %d\n',i,i,i+1);
end
fprintf(fid,'*BEAM SECTION, ELSET=Beams, SECTION=CIRC, MATERIAL=MAT1\n');
fprintf(fid,'0.5\n');
fprintf(fid,'*End Part\n');
fprintf(fid,'*Assembly, name=Assembly\n');
fprintf(fid,'*Instance, name=beam-1, part=GTM\n');
fprintf(fid,'*End Instance\n');
fprintf(fid,'*Nset, nset=fixed_nodes, instance=beam-1\n');
fprintf(fid,'1\n');
fprintf(fid,'*Nset, nset=moved_nodes, instance=beam-1\n\n');
fprintf(fid,'101\n');
fprintf(fid,'**\n');
fprintf(fid,'**\n');
fprintf(fid,'**\n');
fprintf(fid,'**\n');
fprintf(fid,'**\n');
fprintf(fid,'*End Assembly\n');
fprintf(fid,'*MATERIAL, NAME=MAT1\n');
fprintf(fid,'*ELASTIC\n');
fprintf(fid,'10.6e6,0.33\n');
fprintf(fid,'*BOUNDARY\n');
fprintf(fid,'fixed_nodes, 1, 6, 0.\n');
fprintf(fid,'*STEP, NAME=STEP-1, PERTURBATION\n');
fprintf(fid,'*STATIC\n');
fprintf(fid,'*CLOAD\n');
fprintf(fid,'moved_nodes,2,100\n');
fprintf(fid,'*NODE PRINT\n');
fprintf(fid,'U\n');
fprintf(fid,'RF\n');
fprintf(fid,'*END STEP\n');
fclose(fid);