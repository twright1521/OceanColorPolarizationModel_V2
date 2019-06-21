function [Results] =  parseRayXP( inputFileName )


IFILE = fopen(inputFileName);
while (~feof(IFILE)) 
    %Read the line
    line = fgetl(IFILE);
    
    %If you reach a line of asterisks, the parameters section is done
    if(strcmp(line, '**********************************'))
       break; 
    end
end

lineNum = 1;
resultsRaw = cell(1100000, 1);
while (~feof(IFILE)) 
    %Read the line
    line = fgetl(IFILE);
    
    %Save the parameters to a cell array
    resultsRaw{lineNum} = line;
    lineNum = lineNum + 1;
    
end
fclose(IFILE);

%Delete empty lines from the results cell
resultsRaw(find(cellfun('isempty',resultsRaw))) = [];


%% Parse the data section
%Save the detector position and full thickness
%expr = '\s*detector position = (?<detPos>\d*\.\d*)\s*\[\s*full thickness = (?<fullThick>\d*\.\d*)\]';
% expr = '\s*detector position = (?<detPos>(\d*\.\d*)|(\d*))\s*\[\s*full thickness = (?<fullThick>\d*\.\d*)\]';
expr = '\s*detector position = (?<detPos>(\d*\.\d*)|(\d*))\s*\[\s*full thickness = (?<fullThick>(\d*\.\d*)|(\d*))\]';
str = resultsRaw{~cellfun('isempty', regexp(resultsRaw, 'detector position'))};
names = regexp(str, expr, 'names');
Results.DetectorPosition = str2double(names.detPos);
Results.FullThickness = str2double(names.fullThick);

%Split the rest into sections starting with sun zenith angle = 
expr = '\w*Sun zenith angle cosine = ';
[locations] = regexp(resultsRaw, expr);
[sectionStartLines] = find(~cellfun('isempty', locations));

for sect = 1:size(sectionStartLines,1)
    %Read the solar zenith cosine angle
    offsetFromSectionStart = 0;
    line = resultsRaw{sectionStartLines(sect) + offsetFromSectionStart};
    expr = '\w*Sun zenith angle cosine\s*=\s*(?<sZenRad>\d*\.\d*)\s*\[\s*angle\,\s*degree\s*=\s*(?<sZenDeg>\d*\.\d*)\]';
    names = regexp(line, expr, 'names');
    solarZenithRad = str2double(names.sZenRad);
    solarZenithDeg = str2double(names.sZenDeg);
    
    %Read the solar azimuth angle
    offsetFromSectionStart = 1;
    line = resultsRaw{sectionStartLines(sect) + offsetFromSectionStart};
    expr = '\w*azimuth, radian\s*=\s*(?<sAziRad>\d*\.\d*)\s*\[\s*azimuth\,\s*degree\s*=\s*(?<sAziDeg>\d*\.\d*)\]';
    names = regexp(line, expr, 'names');
    solarAzimuthRad = str2double(names.sAziRad);
    solarAzimuthDeg = str2double(names.sAziDeg);
    
    %Read the column headers
    offsetFromSectionStart = 2;
    line = resultsRaw{sectionStartLines(sect) + offsetFromSectionStart};
    headers = regexp(line, '\s*', 'split');
    headers(cellfun('isempty', headers)) = [];   %Clear whitespace (first element returns whitespace)
    headers = regexprep(headers, '/', '_over_');  %/ in headers is not a valid struct field name
    
    
    %Calculate the remaining rows 
    offsetFromSectionStart = 4;
    if(sect < size(sectionStartLines,1))
        rows = sectionStartLines(sect+1)-sectionStartLines(sect) - offsetFromSectionStart;
    else
        rows = size(resultsRaw,1) - sectionStartLines(sect) - offsetFromSectionStart + 1;
    end
    
    %Read the rest of the data
    try
    r = 0:rows-1;
    line = cat(1,resultsRaw{sectionStartLines(sect) + r + offsetFromSectionStart});
    [data, success] = str2num(line);
    
    %Sometimes, when simulating just below the surface some geometries fail
    %for some reason.  The numbers include '#IND' in them.  In this case
    %the str2num operation fails.  Detect this problem and circumvent it.
    if(~success)
        
        %Get the bad rows
        badrowlistlogical = ~cellfun('isempty', regexp(cellstr(line), '#IND'));
        badrowlist = find(badrowlistlogical);
        
        warning('Found bad data in RayXP Simulation: SolZen %u, Azimuth %u, Row %s.\n', ...
            solarZenithDeg, solarAzimuthDeg, sprintf('%u ', badrowlist));
        
        %Initialize data array
        data = nan(size(line,1), 20);
        
        %Convert all good rows
        for eachline = 1:size(line,1)
            %Skip bad rows
            if(any(eachline == badrowlist))
               continue; 
            end            
        
            %Split the row into invdividual columns          
            splitline = regexp(line(eachline,:), '\s{1,2}', 'split');
            splitline(cellfun('isempty', splitline)) = [];
            
            % Copy the line over
            data(eachline,:) = str2double(splitline);
        end
            
         
        %For each bad row, interpolate the data
        for badline = 1:length(badrowlist)
            %Split the row into invdividual columns          
            splitline = regexp(line(badrowlist(badline),:), '\s{1,2}', 'split');
            splitline(cellfun('isempty', splitline)) = [];
             
            for bc = 1:20
                %Get the column value
                thisValue = splitline{bc};
                    
                %Is this value bad?
                badvalue = ~isempty(regexp(thisValue, '#IND'));
                    
                if(badvalue)
                    %Interpolate the bad values
                    angles = data(~badrowlistlogical, 1);
                    OKvalues = data(~badrowlistlogical, bc);
                    interpAngle = data(badrowlist(badline), 1);                    
                    data(badrowlist(badline), bc) = ...
                                interp1(angles, OKvalues, interpAngle, ...
                                        'spline', 'extrap');
                else
                    %If its fine, just add it
                    data(badrowlist(badline), bc) = str2double(thisValue);
                end %End if bad values
            end %End loop over columns
        end %End loop over bad rows            
    end %End if not sucessful conversion
            
    catch er
       error('Something went wrong!  Possible you didnt delete the old .out file?');
       rethrow(er);
    end
    
    %Add the data to the results struct
    for i = 1:size(headers,2)
%         Results.(sprintf('SolZen_%03u', solarZenithDeg)).(['RelAzimuth_', sprintf('%03u', solarAzimuthDeg)]).(headers{i}) = data(:,i);
        Results.(sprintf('SolZen_%03.0f', solarZenithDeg)).(['RelAzimuth_', sprintf('%03.0f', solarAzimuthDeg)]).(headers{i}) = data(:,i);
    end


end
end
%clearvars -except Results
