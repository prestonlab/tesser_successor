function data = tesser_struct_one(header, run)

Screen('Preference', 'SkipSyncTests', 1);

par = header.parameters;
data = struct('subjNum', header.subjNum);
data.type = 'StructLearn_part1';
data.run = run; 
data.header = header;
data.startTime = fix(clock);
addpath(header.path.stim);

rng(sum(100*clock)); %initialize the random number generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keyboard setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET RESPONSE OPTIONS
yeskey = KbName('1!'); %button press if correct orientation .
nopkey = KbName('2@'); %button press if correct orientation

%% RESPONSE BOX CHECK
scanner_respbox;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up results file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up the output file
outfname = sprintf('%s/%s_StructLearn_Part1_Run_%d',header.path.subjinfo,header.subjinfo, run);
if exist([outfname '.txt'],'file') == 2;
    error('The data file for this task already exists.  Please check your input parameters.');
end
fid=fopen([outfname '.txt'], 'a'); % open the file
fprintf(fid, 'SubjNum\t run\t trial\t objnum\t file\t orientnam\t orientnum\t resp\t respnum\t acc\t rt\n'); %Create header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START SCREEN
Screen(window, 'FillRect', par.bgcol, screenrect);
Screen(window, 'TextSize', par.txtsize_instruct);
DrawFormattedText(window,'Press 1 if the object is oriented correctly.\nPress 2 if the object''s orientation is rotated.\nPress any key to begin!',...
    'center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(1);

%wait for button press to start the experiment
clear keyCode; %clearing any prior key presses
clear keyIsDown;
KbWait(-1);

%start time for stimulus onset times
imageFolder = 'stimuli';
objs = par.picshuf;

%if ori+nback task
seqsplit =  par.Seq_Part2a_div;
thiseq = seqsplit{run};

rundata = [];

%Run experimental trials
for t = 1:size(thiseq, 1) %remember this is the number of items in each run
    acctotal = 0;
    
    %% Stimuli details and coordinates
    %Load image
    objnum = thiseq(t, 1); %pick out the t image number from the sequence
    objpos = thiseq(t, 2); %pick out the t image's orientation
    objname = objs{objnum}; %pick out that number of item from the list of shuff object names
    
    if objpos == 0
        img = imread(fullfile(imageFolder,objname)); %read out that image
        orient = 'cor';
        orientnum = 1;
    elseif objpos == 1
        img = imread(fullfile(imageFolder,objname)); %read out that image
        img = imrotate(img, 90); %rotate the image by 90 degreees
        orient = 'rot';
        orientnum = 0;
    end
    
    imageDisplay = Screen('MakeTexture', window, img); %display said image
    
    %Stimulus display time
    maxpictime = par.structTimeout;
    
    %% Actually presenting the stimuli
    pictime = GetSecs; %when the item should begin
    Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
    Screen('DrawTexture', window, imageDisplay); %putting the image on the screen with specified dimensions
    Screen(window, 'Flip'); % Start of trial
    
    %Get keypress response
    RT = 0;
    resp = [];
    
    %response
    respstart = GetSecs; %start of response time
    while GetSecs < pictime + maxpictime
        [keyIsDown, ~, keyCode] = KbCheck(-1); %check for key response
        if keyIsDown
            if keyCode(yeskey)
                resp = 'c';
                respnum = 1;
                respstop = GetSecs;
                RT = respstop - respstart;
            elseif keyCode(nopkey)
                resp = 'n';
                respnum = 2;
                respstop = GetSecs;
                RT = respstop - respstart;
            end %end of keycode check
        end %end of keycode check
    end  %end of kbcheck
    
    if isempty(resp)
        resp = NaN;
        respnum = 0;
        acc = 0;
        RT = NaN;
    end
    
    %CALCULATE ACCURACY HERE%
    if respnum == 1
        if objpos == 0
            acc = 1;
            acctotal = acctotal + 1;
        else
            acc = 0;
        end
    elseif respnum == 2
        if objpos == 1
            acc = 1;
            acctotal = acctotal + 1;
        else
            acc = 0;
        end
    end
    
    %% Save trial into matrix
    %Col1: Subj#
    %Col2: Run#
    %Col3: Trial#
    %Col4: Object#
    %Col5: Rot appearance (1/0)
    %Col6: RotResp(1)/NotRotResp(0)
    %Col7: Rot Accuracy (1/0)
    %Col8: Rot RT
    thistrl = [header.subjNum run t objnum orientnum, respnum, acc, RT];
    rundata = [rundata; thistrl]; %#ok<AGROW>
    
    %% Save results to file
    %Col1: Subj#
    %Col2: Run#
    %Col3: Trial#
    %Col4: Object#
    %Col5: Object File
    %Col6: Rot appearance (rot/not)
    %Col7: Rot appearance num (1/0)
    %Col8: Rot Resp (c/n)
    %Col9: RotResp(1)/NotRotResp(0)
    %Col10: Parse Accuracy (1/0)
    %Col11: Parse RT
    fprintf(fid, '%d\t %d\t %d\t %d\t %s\t %s\t %d\t %s\t %d\t %d\t %f\n',...
        header.subjNum, run, t, objnum, objname, orient, orientnum, resp, respnum, acc, RT);
    
    %% save your matrix stuff
    data.rundata = rundata;
    
    %% output mat file for matrix stuff
    save(sprintf('%s/%s_StructLearn_Part1_Run_%d',header.path.subjinfo,header.subjinfo, run),'data');
    %end of trials
    
end


%prop correct for this block of trials
totalcorr = nansum(data.rundata(:, 7));
propcorr = totalcorr/length(thiseq);
data.propcorr = propcorr;

%adding feedback at end of run
Screen(window, 'FillRect', par.bgcol)
Screen(window, 'TextSize', par.txtsize_instruct);
disptext = sprintf('End of this run!\n\n\nAccuracy was %.2f %.', propcorr); %text to show on the screen
DrawFormattedText(window,disptext,'center','center',par.txtcol); %show the time left
Screen(window, 'Flip');
%final second
WaitSecs(2);

close all
sca;