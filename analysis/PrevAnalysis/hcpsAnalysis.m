% just plotting the data?

clear all; close all;

HCPS_DIR = '/Users/eatai/Google Drive/hcps';
DATA_DIR = [HCPS_DIR filesep 'data'];

hcps = hcpsTools();

%dataDir = uigetdir(DATA_DIR)
dataDir = [DATA_DIR filesep 'sam-tue-jul-25'];

rmErrorTrials = @(fileList) fileList(arrayfun(@(s) isempty(strfind(s.name, '_rst')), fileList));

RT_fileList = dir([dataDir filesep '*r-sos_d-zero*.csv']);
RT_fileList = rmErrorTrials(RT_fileList);

DR_fileList = dir([dataDir filesep '*r-zero_d-sos*.csv']);
DR_fileList = rmErrorTrials(DR_fileList);

% Make trial struct for Ref Tracking
for k = 1:length(RT_fileList)
    trial(k) = hcps.csv2struct(fullfile(dataDir, RT_fileList(k).name));
end

figure(1)
h = hcps.plotstruct(trial)


