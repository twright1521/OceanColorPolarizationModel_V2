function [outFile] = Ray_XP_File_Converter_A1(mie_file,CHL,RayXP_file,Sol_Dir)

% mie_file = 'Simulation_1\RayXP_CHL_1_WS_0_V_1\Mie Files\Sol Files\RayXP_CHL_1_WS_0_V_1_440.033';
% CHL = 1;
% RayXP_file = 'PTH_mmVF.440';
% Sol_Dir = 'C:\Users\Trevor Wright\Documents\Research\RayXP\RAYdelivery_1\Ray\SolLib';

% Get RayXP file and create a data array
fid = fopen(strcat(Sol_Dir,filesep,RayXP_file));
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = c{1};
RayXP_data = [];
for z = 1:length(lines)
    if ~isempty(lines{z})
        parsenumbers = textscan(lines{z},'%f','delimiter',',');
        numbers = parsenumbers{1};
        if ~isempty(numbers)
            RayXP_data = [RayXP_data; numbers.'];
        end
    end
    clearvars parsenumbers numbers
end

clearvars lines c fid

% Get RayXP file header info
fid = fopen(strcat(Sol_Dir,filesep,RayXP_file));
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = c{1};
Header_lines = lines(1:18,1);

clearvars lines c fid

% Get Mie file and create a data array
fid = fopen(mie_file);
c = textscan(fid,'%s','delimiter','\n');
fclose(fid);
lines = c{1};
Mie_data = [];
for z = 1:length(lines)
    if ~isempty(lines{z})
        parsenumbers = textscan(lines{z},'%f','delimiter',',');
        numbers = parsenumbers{1};
        if ~isempty(numbers)
            Mie_data = [Mie_data; numbers.'];
        end
    end
    clearvars parsenumbers numbers
end

clearvars lines c fid

% Calculate new data set
New_RayXP_data = RayXP_data(:,1);
for x = 2:size(RayXP_data,2)
    New_RayXP_data(:,x) = Mie_data(:,2).*(RayXP_data(:,x)./RayXP_data(:,2));
end

% Output file
k = textscan(RayXP_file, '%s', 'delimiter', '.');
rayxp_name = k{1};
c = textscan(mie_file, '%s', 'delimiter', '\');
mie_name = c{1};

clearvars c k

if length(mie_name) == 1
    new_file_path = [];
else
    for x = 1:length(mie_name)-1
        if x == 1
            new_file_path = mie_name{1};
        else
            new_file_path = strcat(new_file_path,filesep,mie_name{x});
        end
    end
end

new_file_name = strcat('Scat_Mat_Normalized_',rayxp_name{1},'_CHL_',num2str(CHL),'.',rayxp_name{end});

% Writing New File
outFile = strcat(new_file_path,filesep,new_file_name);
FID  = fopen(outFile, 'wt');

% Header info from rayXP_file
for x = 1:length(Header_lines)
    fprintf(FID, '%s\n', Header_lines{x});
end

% Print Data array

for x = 1:size(New_RayXP_data,1)
    for y = 1:size(New_RayXP_data,2)
        if y == size(New_RayXP_data,2)
            fprintf(FID, '%16.9e \n', New_RayXP_data(x,y));
        else
            fprintf(FID, '%16.9e, ', New_RayXP_data(x,y));
        end
    end
end

% Print ancillary information at the bottom
fprintf(FID, 'TEXT=\n');
fprintf(FID, 'Converted %s using %s to Generate this new scattering matrix\n',RayXP_file,mie_name{end});
fprintf(FID, 'END_TEXT\n');

fclose(FID);

copyfile(outFile,Sol_Dir)
end