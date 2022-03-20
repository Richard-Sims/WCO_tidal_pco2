function [Std_nom_values, Std_runtimes] = getCO2calSeasonalStudy(fullpathLicorcal,j,fList)
%%%%%%%%%%%%%%%%%%%   CO2 calibration import and linear fitting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if fList(j,1:10)==('2016-04-27');
    delimiter = '\t';
    formatSpec = '%s%s%[^\n\r]';
    fileID = fopen(fullpathLicorcal,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
    fclose(fileID);
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = dataArray{col};
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));
    rawData = dataArray{1};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, 1) = numbers{1};
                raw{row, 1} = numbers{1};
            end
        catch me
        end
    end
    rawNumericColumns = raw(:, 1);
    rawCellColumns = raw(:, 2);
    R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
    rawNumericColumns(R) = {NaN}; % Replace non-numeric cells
    Std_nom_values = cell2mat(rawNumericColumns(:, 1));
    Std_runtimes = rawCellColumns(:, 1);
    clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns R;
else
    delimiter = '\t';startRow = 2;formatSpec = '%*s%s%[^\n\r]'; fileID = fopen(fullpathLicorcal,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines' ,startRow-1, 'ReturnOnError', false);
    fclose(fileID);rawNumericColumns = dataArray{:, 1};
    Std_runtimes=cell2mat(rawNumericColumns(:, 1));Std_runtimes=cellstr(Std_runtimes);
    Std_nom_values = [250,450,250,450]';
end
        