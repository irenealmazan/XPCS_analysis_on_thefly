% info to analyze scans 51 - 64 of 2018_0813_1 file


 % temperatures in C
TCV = [850;850;850;850;850]; % Equilibrium Thermocouple (new heater)

% power in Watts
PWV = [0; 10 ; 15 ; 20 ; 25];

% Filenames
specfilenameM = ['2018_0813_1';'2018_0813_1';'2018_0813_1';'2018_0813_1';'2018_0813_1'];

% scans number (see lab book)
SCNstrM = ['051';'054';'057';'061';'064']; % Equilibrium



XCENV = [93;93;93;93;93];% del
YCENV = [138;138;138;138;138]; % nu


% Width of the larger ROI
XWIDV = [25;25;25;25;25];
% Pilatus detector on sevchex arm (X is del and Y is nu)
YWIDV = [45;45;45;45;45];


% Offset of two of the ROIS around the Positions of CTRs
ymax = [8;8;8;8;8];

% Min and max on time range for delta-time average
tminv = [100; 100; 100;100;100];
tmaxv = [4000;4000;4000;4000;4000];
