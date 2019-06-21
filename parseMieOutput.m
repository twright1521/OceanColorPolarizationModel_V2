function [  F11, F21, F33, F34, ScattAngle, ...
            alpha1, alpha2, alpha3, alpha4, beta1, beta2, ...
            Wavelen_nm, Rem, Imm, SSA, Qsca, Qext, G, Vavg] = ...
                                      parseMieOutput( mieFileName )
%parseMieOutput Parses the results of the Meerhoff Mie code

%Check that the mieFilename exists
if(exist(mieFileName, 'file') == 0)
    error('Filename: %s doesnt exist!\n', mieFileName);
end

%Read the file into a cell array
FID = fopen(mieFileName, 'rt');

inputFile = cell(5000,1);
linenum = 1;
while(~feof(FID))
    
    inputFile{linenum, 1} = fgetl(FID);
    
    linenum = linenum + 1;
end

%Close the input file
fclose(FID);

%Remove unused lines from the struct
inputFile(cellfun('isempty', inputFile), :) = [];

%Determine whether this is only a file with expansion coefficients, or
%contains both
if(strcmpi(inputFile{1}, '1OUTPUT OF THE MEERHOFF MIE PROGRAM VERSION 3.0'))
    fileHasDirectF = true;
elseif(strcmpi(inputFile{1}, ' EXPANSION COEFFICIENTS SCATTERING MATRIX'))
    fileHasDirectF = false;
else
    error('I dont understand this file header!\n');
end

if(fileHasDirectF)
    %Search for first line which contains scat.angle
    %These angles are a higher resolution than the ones at the bottom.
    headerlinematch = regexp(inputFile, 'scat.angle', 'match');
    headerLineNum = find(~cellfun('isempty', headerlinematch), 1);

    %Define start and end lines for the Scat Matrices
    dataStartLine = headerLineNum+2;

    %Search for the line which contains alpha1
    coeffheaderlinematch = regexp(inputFile, 'alpha1', 'match');
    coeffheaderLineNum = find(~cellfun('isempty', coeffheaderlinematch));
    dataEndLine = coeffheaderLineNum - 1;  %Last line doesn't count.

    %Convert the string data to numbers
    strinput = strcat(inputFile(dataStartLine:dataEndLine, :), ';');
    ScatMatrix = str2num(cat(1,strinput{:}));

    %Break out the individual components
%     ScattAngle  = ScatMatrix(:,2);
%     F11         = ScatMatrix(:,3);
%     F21         = ScatMatrix(:,4);
%     F33         = ScatMatrix(:,5);
%     F34         = ScatMatrix(:,6);
else
%     [ScattAngle, F11, F21, F33, F34] = deal(NaN);
%     warning('No Direct Scattering Matrix found!\n');
end

    %Grab the expansion coefficients
    %Search for the line which contains alpha1
    coeffheaderlinematch = regexp(inputFile, 'alpha1', 'match');
    coeffheaderLineNum = find(~cellfun('isempty', coeffheaderlinematch));

    %Define start and end lines for the legendre coefficients
    dataStartLine = coeffheaderLineNum+1;
    
    %Find the lines with scat.angle again.  If this is an output file that
    %also has direct Fs, there will be 2, otherwise there will be none
    headerlinematch = regexp(inputFile, 'scat.angle', 'match');
    headerLineNum = find(~cellfun('isempty', headerlinematch));
    if(length(headerLineNum) == 2)
        dataEndLine = headerLineNum(2)-1;
    else
        %This file contains only coefficients, so you can go to the end.
        dataEndLine = size(inputFile, 1) - 1;  %Last line doesn't count.
    end

    %Convert the string data to numbers
    strinput = strcat(inputFile(dataStartLine:dataEndLine, :), ';');
    
    %Replace any equal signs with nothing
    %(only applied to files containing direct Scat Matrices)
    strinput = strrep(strinput,'l=', '');
    
    %Turn it into a character array
    strinput = cat(1,strinput{:});
    

    Coefficients = str2num(strinput);

    %Break out the individual components
    alpha1 = Coefficients(:,2);
    alpha2 = Coefficients(:,3);
    alpha3 = Coefficients(:,4);
    alpha4 = Coefficients(:,5);
    beta1  = Coefficients(:,6);
    beta2  = Coefficients(:,7);
       
    %Precompute the Wigner-d functions
    ScattAngle = [ 0:0.05:1, 1:0.1:5, 5:0.25:10, 10:0.35:179, 179:0.1:180];
    ScattAngle = unique(ScattAngle)';  %Remove duplicates because I'm lazy
    
    s    = size(alpha1, 1)-1;
    d00  = WignerD(0, 0, 0:s, ScattAngle);
    d02  = WignerD(0, 2, 2:s, ScattAngle);
    d22  = WignerD(2, 2, 2:s, ScattAngle);
    d2m2 = WignerD(2,-2, 2:s, ScattAngle);
    
    %Convert the coefficients back into scattering matrix elements
    %(Better than reading the values directly!)
    OnesAng     = ones(1,size(d00,2));
    F11         = sum((alpha1*OnesAng).*d00, 1)';
    F11plusF33  = sum(((alpha2(3:end)+alpha3(3:end))*OnesAng).*d22, 1)';
    F11minusF33 = sum(((alpha2(3:end)-alpha3(3:end))*OnesAng).*d2m2, 1)';
    F21         = -1.*sum((beta1(3:end)*OnesAng).*d02, 1)';
    F34         = -1.*sum((beta2(3:end)*OnesAng).*d02, 1)';    
    F33         =  0.5.*(F11plusF33 - F11minusF33);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parse the extra data

%Cut the input to the first 60 lines (to increase speed)
headers = inputFile(1:60);
% headers = inputFile(1:30);


%%%%%%%%%%%%%%%%%% Extract the wavelength %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for "wavelength"
expr = 'wavelength';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));

%If not found, search for 'lambda'
if(~any(matchedRow))
    % lambda=   0.6660000 Re(m) =   1.3880000 Im(m) =   0.0000000
    expr = 'lambda';
    matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
    line = headers{matchedRow};
    numOnly = regexprep(line, '[^\d.\s]', '');
    doublevals = str2num(numOnly);
    
    %Assign the wavelength, Rem, and Imm
    Wavelen_nm = doublevals(1)*1000;
    
    %If there is no Rem or Imm, this is a mixture file, so put NaN for Rem
    %and Imm
    if(length(doublevals) >= 3)
        Rem = doublevals(2);
        Imm = doublevals(3);
    else
        Rem = NaN;
        Imm = NaN;
    end
    
    HaveRemImmAlready = true;
else
    % wavelength                       0.6660000   0.0000000   0.0000000   0.0000000 
    line = headers{matchedRow};
    
    %Extract the numbers only.
    numOnly = regexprep(line, '[^\d.\s]', '');
    doublevals = str2num(numOnly);
    
    %Assign the wavelength
    Wavelen_nm = doublevals(1)*1000;
    
    HaveRemImmAlready = false;
end

%%%%%%%%%%%%%%%%%% Extract the Rem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If we don't have the Rem and Imm already, look for them
if(~HaveRemImmAlready)
    %Search for Rem
    expr = 'Re\(m\)';
    matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
    line = headers{matchedRow};
    
    %Extract the numbers only.
    numOnly = regexprep(line, '[^\d.\s]', '');
    doublevals = str2num(numOnly);
    
    %Assign the Rem
    Rem = doublevals(1);
    
    %Search for Imm
    expr = 'abs\(Im\(m\)\)';
    matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
    line = headers{matchedRow};
    
    %Extract the numbers only.
    numOnly = regexprep(line, '[^\d.\s]', '');
    doublevals = str2num(numOnly);
    
    %Assign the Rem
    Imm = doublevals(1);    
end

%%%%%%%%%%%%%%%%%% Extract the SSA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Search for the term 'single scattering albedo'
expr = 'single scattering albedo';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
line = headers{matchedRow};

%Extract the numbers only.
numOnly = regexprep(line, '[^\d.\sE+-]', '');
doublevals = str2num(numOnly);

%Assign the SSA
SSA = doublevals(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract Qsca %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for the term 'Qsca'
expr = 'Qsca';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
line = headers{matchedRow};

%Extract the numbers only.
numOnly = regexprep(line, '[^\d.\sE+-]', '');
doublevals = str2num(numOnly);

%Assign the Qsca
Qsca = doublevals(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract Qext %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for the term 'Qext'
expr = 'Qext';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
line = headers{matchedRow};

%Extract the numbers only.
numOnly = regexprep(line, '[^\d.\sE+-]', '');
doublevals = str2num(numOnly);

%Assign the Qext
Qext = doublevals(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract G %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for the term 'geometrical cross section'
expr = 'geometrical cross section';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
line = headers{matchedRow};

%Extract the numbers only.
numOnly = regexprep(line, '[^\d.\sE+-]', '');
doublevals = str2num(numOnly);

%Assign the G
G = doublevals(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Extract Vavg %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Search for the term 'average volume'
expr = 'average volume';
matchedRow = ~cellfun('isempty', regexp(headers, expr, 'match'));
line = headers{matchedRow};

%Extract the numbers only.
numOnly = regexprep(line, '[^\d.\sE+-]', '');
doublevals = str2num(numOnly);

%Assign the Vavg
Vavg = doublevals(1);

end

