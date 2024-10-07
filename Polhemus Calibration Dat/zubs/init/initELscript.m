% Initilize eyetracker
global wndw
%addpath('C:\toolbox\Psychtoolbox\PsychBasic\MatlabWindowsFilesR2007a\'); % Add path to Eyelink.mexw64

if Eyelink('Initialize') ~= 0
	fprintf('Problem initializing eyelink\n');
	return;
end

el=EyelinkInitDefaults(wndw);

dummymode=0;
if ~EyelinkInit(dummymode)
	fprintf('Eyelink Init aborted.\n');
	cleanup;  % cleanup function
	return; end

%Eyelink('command','screen_pixel_coords = %ld %ld %ld %ld', 0, 0, Sexp.resx-1, Sexp.resy-1);
Eyelink('command','screen_pixel_coords = %d, %d, %d, %d', 0, 0, Sexp.resx-1, Sexp.resy-1);
%Eyelink('message','DISPLAY_COORDS %ld %ld %ld %ld', 0, 0, Sexp.resx-1, Sexp.resy-1);
Eyelink('message','DISPLAY_COORDS %d, %d, %d, %d', 0, 0, Sexp.resx-1, Sexp.resy-1);
%%%
%%%
Eyelink('command',sprintf('calibration_area_proportion = %f %f',.35,.35))
Eyelink('command',sprintf('validation_area_proportion = %f %f',.3,.3))
%%%
%%%
Eyelink('command','calibration_type = HV13');
Eyelink('command', 'file_sample_data  = LEFT,RIGHT,GAZE,HREF,AREA,GAZERES,STATUS,INPUT');
Eyelink('command', 'link_sample_data  = LEFT,RIGHT,GAZE,AREA');

% 
% mglEyelinkCMDPrintF('screen_pixel_coords = 0, 0, %d, %d',...
%             mglGetParam('screenWidth'), mglGetParam('screenHeight'));
% 
% mglEyelinkCMDPrintF(sprintf('calibration_area_proportion = %f %f',eyelinkParams.calibrationAreaX,eyelinkParams.calibrationAreaY));
% mglEyelinkCMDPrintF(sprintf('validation_area_proportion = %f %f',eyelinkParams.calibrationAreaX,eyelinkParams.calibrationAreaY));


% make sure we're still connected.
if Eyelink('IsConnected')~=1 && dummymode == 0
	fprintf('not connected, clean up\n');
	Eyelink( 'Shutdown');
	Screen('Close');
	cleanup
	return; end

el.backgroundcolour=bgcolor;
el.calibrationtargetcolour=BLU;
el.calibrationtargetsize=1.25;
el.calibrationtargetwidth=.75;

% parameters are in frequency, volume, and duration
% set the second value in each line to 0 to turn off the sound
el.cal_target_beep=[600 0.5 0.05];
el.drift_correction_target_beep=[600 0.5 0.05];
el.calibration_failed_beep=[400 0.5 0.25];
el.calibration_success_beep=[800 0.5 0.25];
el.drift_correction_failed_beep=[1 0 0];% [400 0.5 0.25];
el.drift_correction_success_beep=[1 0 0]; %[800 0.5 0.25];

EyelinkUpdateDefaults(el); %apply the changes from above

count=0; t0=GetSecs;
% while count<1111,
% 	count=count+(Eyelink('NewFloatSampleAvailable')>0); end
Sexp.ELnominalrate=ceil(count/(GetSecs-t0));