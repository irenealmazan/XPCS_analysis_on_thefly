% information for scan 100 file 2018_0813_1

% temperatures in C
TCV = [800]; % Equilibrium Thermocouple (new heater)

% power in Watts
PWV = [15];

% Filenames
specfilenameM = ['2018_0813_1'];

% scans number (see lab book)
SCNstrM = ['100']; 

% Center pixels in between the CTRS
%XCENV = [230;236;240;240];% del
%YCENV = [120;118;120;120]; % nu

XCENV = [93];% del
YCENV = [138]; % nu



% Width of the larger ROI
XWIDV = [25];
% Pilatus detector on sevchex arm (X is del and Y is nu)
YWIDV = [45];


% Offset of two of the ROIS around the Positions of CTRs
ymax = [8];

% Min and max on time range for delta-time average
tminv = [1000];
tmaxv = [4000];
%tminv = [720; 2500;1110;200];
%tmaxv = [1400;3500;2000;600];
