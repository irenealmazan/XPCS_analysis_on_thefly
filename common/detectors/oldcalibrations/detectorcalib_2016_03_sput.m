function DetCalib = detectorcalibrationfile

%%% Use the ROI/epics or imageJ convention for X and Y and indexing of pixels
%%% Note - indexing starts at ZERO for imageJ or the ROI/epics.
%%%  size(X, Y) and indexing for length N, ---  [0 : (N-1)]
%%% For the matlab programs that use these remember for matrices/images
%%%			(matlab) X = columns Y=rows, and indexing starts at ONE
%%% Matlab, size(Y(rows),X(columns)), and length N ---- [1 : N]
DetCalib.docu = strvcat(...
'2016-02_Pilatus file working on structure CT 2016-2 version compatible with rotated detector and non rectangular');
%% I am using XY below to remind myself it isn't in column/row 
DetCalib.pixXYSIZE	= [487 195];   % [X Y] from imageJ 
%	Xcen  Ycen
DetCalib.pixXYCEN		= [122 92];
% relative scale	  [Xscale Yscale  XYangle]
DetCalib.pixXYSCALE	= [1 1 90];

DetCalib.ADU2photon 	= -1; 
%% ADU2photon - Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) 
%%		(note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons

%%  In terms of the angles (vertical) and (horizontal) 
DetCalib.HVbeam		= [0   0];  % [Hor  Vert] angle where beam would be in defined center ROI (this allows scattered beam to be used, also takes account that direct beam is at delta(Vertical angle) = -rho on sevchx
DetCalib.Hname 		= 'Nu';
DetCalib.HRotType	= 'Rot';
% Vcal 2015_0208_3 #5 nu calibration scan about beam zero 
% These are from old one need to find newere checks, but assume close enough
%					  [angle	Xpix	Ypix]
DetCalib.Hcal		= [ -0.50	122		 92-47
						+0.50	122		 92+47];


DetCalib.Vname 		= 'Delta';
DetCalib.VRotType	= 'Rot';
% Vcal 2015_0208_3 #4 delta calibration scan about beam zero 
%					  [at_angle	Xpix	Ypix]
DetCalib.Vcal		= [ -0.50	122+47		92
						+0.50	122-47		92];
						
					
DetCalib.ADU2photon 	= -1; 

%% ADU2photon
%%	Calibrated efficiency of 'counts' in pixel to photons; 
%%		-1 (negative number) denotes uncertain, uncalibrated
%% 		+1  (1 count in pixel is 1 photon) (note: Pilatus are not 1 to 1 efficient at high energies.
%% 		0.5 means 1 count in pixel per 2 photons

%#CR ROI  2: MinX= 116  SizeX=   5 MinY=  78  SizeY=  12  
%   Xscan is 116 117 118 119 120   Y 78 79 80 81 82 83 84 85 86 87 88 89 
% #CR ROI  2: MinX= 121  SizeX=   3 MinY=  86  SizeY=  11
end


