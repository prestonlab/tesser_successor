function data = tesser_ori(header)

Screen('Preference', 'SkipSyncTests', 1);

par = header.parameters;
data = struct('subjNum', header.subjNum);
data.header = header;
data.startTime = fix(clock);
addpath(header.path.stim);

rng(sum(100*clock)); %initialize the random number generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Screen setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCREEN SETUP
allscreens = Screen('screens'); % How many screens are connected to the computer?
screennum = max(allscreens); % Choose the last one (not the main screen)
[window, screenrect] = Screen('OpenWindow', screennum, par.bgcol);
screenLeft = screenrect(RectLeft); %pixel value for left of screen
screenRight = screenrect(RectRight); %pixel value for right of screen
screenTop = screenrect(RectTop); %pixel value for top of screen
screenBottom = screenrect(RectBottom); %pixel value for bottom of screen
screenWidth = screenRight - screenLeft; %pixel value for screen width
screenHeight = screenBottom - screenTop; %pixel value for screen height
ycenter = (screenBottom - screenTop)/2; %ycenter pixel value
xcenter = (screenRight - screenLeft)/2; %xcenter pixel value
xoffset_text = xcenter - 170; %for old new judgment text
yoffset_text = ycenter + 220;

Screen(window, 'TextFont', 'Arial');
HideCursor;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START SCREEN
Screen(window, 'TextSize', 50);
DrawFormattedText(window,'Observe the orientation of these objects carefully.\nPress any key to begin!',...
    'center','center',par.txtcol);
Screen(window, 'Flip');
Screen(window, 'TextSize', 50);
WaitSecs(1);

%wait for button press to start the experiment
clear keyCode;
clear keyIsDown;
KbWait(-1);

%start time for stimulus and image folders
imageFolder = 'stimuli';
objs = par.pics;
reps = 2; %will see two repeitions of this list

for r = 1:reps
    %Run experimental trials
    for t = 1:length(objs) %remember this is the number of items in list
        
        %% Stimuli details and coordinates
        %Load image
        objname = objs{t}; %pick out that number of item from the list of shuff object names
        
        img = imread(fullfile(imageFolder,objname));
        imageDisplay = Screen('MakeTexture', window, img); %display said image
        
        %% Actually presenting the stimuli
        Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
        Screen('DrawTexture', window, imageDisplay); %putting the image on the screen with specified dimensions
        Screen('Flip', window); % Start of trial
        WaitSecs(3);
        
        %ITI b/w the stimuli
        Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
        Screen(window, 'TextSize', par.txtsize_fix)
        DrawFormattedText(window,'+','center','center',par.txtcol);
        Screen('Flip', window); % Start of trial
        WaitSecs(.5);
        
    end
end

DrawFormattedText(window,'Great job!','center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(2);
sca;

