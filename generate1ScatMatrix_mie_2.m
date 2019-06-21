function [ outputFilePath ] = generate1ScatMatrix_mie_2(outputFileName,file_path,varargin)
%generate1ScatMatrix_mie Generates a single scattering Matrix in the form
%that Jacek's GAP program expects.
%   Input Description needed.  Takes the same inputs as the 
%   writeMie3InputFile function.

%Define folders to store files
logFileDir = strcat(file_path,filesep,'Mie Output');
outputDir = strcat(file_path,filesep,'Expansion Coeffs');
inputDir = strcat(file_path,filesep,'Mie Input');

%If these folders don't exist, create them.
if ~exist(logFileDir, 'dir') 
    mkdir(logFileDir);  
end

if ~exist(outputDir,  'dir')
    mkdir(outputDir);   
end

if ~exist(inputDir,   'dir')
    mkdir(inputDir);    
end

%Remove pre-existing mie_sc files 
if exist('mie_sc','file')
    delete('mie_sc');
end

%Write the input file
inputFile = writeMie3InputFile(varargin{:});

%Move the input file to the inputfile directory
inputFilePath = strcat(inputDir,filesep, outputFileName);
movefile(inputFile, inputFilePath);

%Perform the Mie calculation
outputFilePath = strcat(logFileDir, filesep, outputFileName);
cmd = strcat('"',pwd,filesep,'mie3.exe" < ','"',inputFilePath,'"',' > ','"',outputFilePath,'"');
system(cmd);

end

