% mqlLogParser function
% Matlab function for reading in and parsing MaxQuant.Live log files from targeting runs. Takes
% a log file filename, and splits the output into a table covering peptides
% identified, as well as mz, rt and intensity deviations.
% 
% Example function call:
%   log = mqlLogParser('filename.txt');
%
% Note: this function has only been used to analyse output log files
% generated with the maxquant.live targeting app so may not work or provide
% useful output for the other maxquant.live apps.
%
% Ed Emmott, May 2019, Northeastern U.

% Todo: script can take a couple of minutes to run. regexp conversion of 
% the 'continued processing' to replace the slow loop with faster data 
% splitting is planned.

function [mqLiveLog] = mqlLogParser(filename);

%% Import files
    delimiter = '\t';
    if nargin<=2
        startRow = 1;
        endRow = inf;
    end

    formatSpec = '%s%s%s%s%s%[^\n\r]';

    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for block=2:length(startRow)
        frewind(fileID);
        dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        for col=1:length(dataArray)
            dataArray{col} = [dataArray{col};dataArrayBlock{col}];
        end
    end
    fclose(fileID);

    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    rawData = dataArray{4};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
        
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, 4) = numbers{1};
                raw{row, 4} = numbers{1};
            end
        catch
            raw{row, 4} = rawData{row};
        end
    end

    % Convert the contents of columns with dates to MATLAB datetimes using the
    % specified date format.
    try
        dates{1} = datetime(dataArray{1}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
    catch
        try
            % Handle dates surrounded by quotes
            dataArray{1} = cellfun(@(x) x(2:end-1), dataArray{1}, 'UniformOutput', false);
            dates{1} = datetime(dataArray{1}, 'Format', 'yyyy-MM-dd', 'InputFormat', 'yyyy-MM-dd');
        catch
            dates{1} = repmat(datetime([NaN NaN NaN]), size(dataArray{1}));
        end
    end

    anyBlankDates = dataArray{1} == '';
    anyInvalidDates = isnan(dates{1}.Hour) - anyBlankDates;
    dates = dates(:,1);

    % Split data into numeric and string columns.
    rawNumericColumns = raw(:, 4);
    rawStringColumns = string(raw(:, [2,3,5]));

    % Exclude rows with non-numeric cells
    I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),rawNumericColumns),2);  %Find rows with non-numeric cells
    K = I | anyInvalidDates | anyBlankDates;
    dates = cellfun(@(x) x(~K,:), dates, 'UniformOutput', false);
    rawNumericColumns(K,:) = [];
    rawStringColumns(K,:) = [];

    % Make sure any text containing <undefined> is properly converted to an <undefined> categorical
    idx = (rawStringColumns(:, 2) == "<undefined>");
    rawStringColumns(idx, 2) = "";

    Untitled = table;
    Untitled.Date = dates{:, 1};
    Untitled.Time = rawStringColumns(:, 1);
    Untitled.ScanProtocol = categorical(rawStringColumns(:, 2));
    Untitled.RT = cell2mat(rawNumericColumns(:, 1));
    Untitled.Scan = rawStringColumns(:, 3);

    mqLiveLog     = Untitled;

    %% mqLiveLog Processing
    % Filter on 'Scan protocol'
    mqLiveLog.ScanProtocol = categorical(mqLiveLog.ScanProtocol);
    mqLiveLog = mqLiveLog(mqLiveLog.ScanProtocol=='ScanProtocol',:);
 
 
 % Remove all 'NEW SCAN ARRIVED' rows
    for ii = 1:numel(mqLiveLog(:,1))
        isNewScan(ii,1) =  strcmp(mqLiveLog.Scan{ii}(1:16), 'NEW SCAN ARRIVED');
    end
 
    mqLiveLog = mqLiveLog(~isNewScan,:);

%% Continued processing
    clear mqLiveLogStruct
    

    mqLiveLogStruct = struct();
    for ii = 1:numel(mqLiveLog(:,1))
        mqLiveLogStruct.under.(['logline_',num2str(ii)]) = strfind(mqLiveLog.Scan{ii} , '_');
        mqLiveLogStruct.space.(['logline_',num2str(ii)]) = strfind(mqLiveLog.Scan{ii} , ' ');
        mqLiveLogStruct.colon.(['logline_',num2str(ii)]) = strfind(mqLiveLog.Scan{ii} , ':');
    end

numrows = numel(mqLiveLog(:,1));

f = waitbar(0,'Starting log parsing');

    for ii = 1:numrows
         waitbar(ii/numrows , f , ['Parsing Log File Rows (',num2str(numrows),')']);
        mqLiveLog.id(ii)          = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(1):mqLiveLogStruct.under.(['logline_',num2str(ii)])(2)-1));
        mqLiveLog.Sequence{ii}    = mqLiveLog.Scan{ii}(mqLiveLogStruct.under.(['logline_',num2str(ii)])(2):mqLiveLogStruct.under.(['logline_',num2str(ii)])(3));
        mqLiveLog.Mass(ii)        = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(3)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(4)-1));
        mqLiveLog.z(ii)           = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(5)+1));
        mqLiveLog.mz(ii)          = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(7)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(8)-1));
        mqLiveLog.mzDeviation(ii) = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(9)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(10)-1));
        mqLiveLog.mzCorrection(ii)= str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(11)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(12)-1));
        mqLiveLog.mzWindow(ii)    = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(13)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(14)-1));
        %rt % Ignoring as have the earlier rt.
        mqLiveLog.rtDeviation(ii) = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(17)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(18)-1));
        mqLiveLog.rtCorrection(ii)= str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(19)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(20)-1));
        mqLiveLog.rtWindow(ii)    = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(21)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(22)-1));
        mqLiveLog.intensity(ii)   = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(23)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(24)-1));
        mqLiveLog.intRatio(ii)    = str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(25)+1:mqLiveLogStruct.space.(['logline_',num2str(ii)])(26)-1));
        mqLiveLog.intCorrection(ii)= str2num(mqLiveLog.Scan{ii}(mqLiveLogStruct.space.(['logline_',num2str(ii)])(27)+1:end));
    
    end
close(f)
end