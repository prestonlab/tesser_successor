function data = tesser_quest(header)

Screen('Preference', 'SkipSyncTests', 1);

%% ready, set, go
par = header.parameters;
data = struct('subjNum', header.subjNum);
data.type = 'InductGen'; 
data.header = header;
data.startTime = fix(clock);
addpath(header.path.stim);

rng(sum(100*clock)); %initialize the random number generator

%% SET UP RESULTS FILE
%Set up the output file
outfname = sprintf('%s/%s_InductGen',header.path.subjinfo,header.subjinfo);
if exist([outfname '.txt'],'file') == 2;
    error('The data file for this task already exists.  Please check your input parameters.');
end
fid=fopen([outfname '.txt'], 'a'); % open the file
fprintf(fid, 'SubjNum\t TrialNum\t QuestType\t Environment\t CueNum\t CueObject\t Opt1Num\t Option1\t Opt2Num\t Option2\t Resp\t Acc\t RT\n'); %Create header


%% SCREEN SETUP
allscreens = Screen('screens'); % How many screens are connected to the computer?
screennum = max(allscreens); % Choose the last one (not the main screen)
[window, screenrect] = Screen('OpenWindow',screennum,par.bgcol);
screenLeft = screenrect(RectLeft); %pixel value for left of screen
screenRight = screenrect(RectRight); %pixel value for right of screen
screenTop = screenrect(RectTop); %pixel value for top of screen
screenBottom = screenrect(RectBottom); %pixel value for bottom of screen
screenWidth = screenRight - screenLeft; %pixel value for screen width
screenHeight = screenBottom - screenTop; %pixel value for screen height
ycenter = (screenBottom - screenTop)/2; %ycenter pixel value
xcenter = (screenRight - screenLeft)/2; %xcenter pixel value
Screen(window, 'TextSize', 45);
Screen(window, 'TextFont', 'Arial');

%% IMAGE SIZING
objheight = par.objheight*(3/4); %height in pixels
objwidth = par.objwidth*(3/4); %width in pixels

%cue
Cx1 = screenWidth/2 - objwidth/2; %left side of image
Cy1 = screenHeight/4 - objheight/2; %top of image
Cx2 = screenWidth/2 + objwidth/2; %right side of image
Cy2 = screenHeight/4 + objheight/2; %bottom of image
cuerect = [Cx1 Cy1 Cx2 Cy2]; %rectangular coordinates for the cue location

%object 1
O1x1 = screenWidth/4 - objwidth/2;
O1y1 = screenHeight*(3/4) - objheight/2;
O1x2 = screenWidth/4 + objwidth/2;
O1y2 = screenHeight*(3/4) + objheight/2;
obj1rect = [O1x1 O1y1 O1x2 O1y2];

%object 2
O2x1 = screenWidth*(3/4) - objwidth/2;
O2y1 = screenHeight*(3/4) - objheight/2;
O2x2 = screenWidth*(3/4) + objwidth/2;
O2y2 = screenHeight*(3/4) + objheight/2;
obj2rect = [O2x1 O2y1 O2x2 O2y2];

%% SET RESPONSE OPTIONS
leftkey1 = KbName('1'); %button to select the left object
leftkey2 = KbName('1!');
rightkey1 = KbName('2'); %button to select the right object
rightkey2 = KbName('2@');

%% START SCREEN
Screen(window, 'FillRect', par.bgcol, screenrect);
Screen(window, 'TextSize', par.txtsize_instruct);
DrawFormattedText(window,'Press 1 for the left object.\nPress 2 for the right object.\nPress any key to begin!',...
    'center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(1);


%wait for button press to start the experiment
clear keyCode;
clear keyIsDown;
KbWait(-1);

%where the stimuli are located
imageFolder = 'stimuli';
objs = par.picshuf;

%all of the items in a community
ocean_comm = par.comm1; %ocean community
ocean_comm = strread(num2str(ocean_comm), '%s'); %ocean community into strings

desert_comm = par.comm2; %desert community
desert_comm = strread(num2str(desert_comm), '%s'); %desert community into strings

forest_comm = par.comm3; %forest community
forest_comm = strread(num2str(forest_comm), '%s'); %forest community into strings

rundata = [];

for i = 1:length(par.induct)
    %timing (max time is 8s for each question; 0.5s for fixation after)
    
    %onsettime = par.onset.induc(1, i);
    pictime = par.inductTimeout;
    fixtime = par.inductTimefix;
    
    %the type of induction question
    questtype = par.induct{i, 1};
    questtypenum = par.induct{i, 9}; 
    
    %the cue
    cuenum =  par.induct{i, 2}; %pick out the t image number from the sequence
    cuename = objs{cuenum}; %pick out that number of item from the list of shuff object names
    
    %option 1 %all of the correct choices of now
    opt1num = par.induct{i, 3};
    opt1name = objs{opt1num}; %#ok<NASGU>
    opt1pos = par.induct{i, 5}; %this is how the correct answer gets switched around later
    
    %option2
    opt2num = par.induct{i, 4};
    opt2name = objs{opt2num}; %#ok<NASGU>
    
    %figure out how they will be cued with each objec thing here:
    %Cue:
    cuestim = imread(fullfile(imageFolder,cuename));
    
    %Multiple Choices:
    %choice 1 (bottom left)
    obj1stim = imread(fullfile(imageFolder,opt1name));
    
    %choice 2 (bottom right)
    obj2stim = imread(fullfile(imageFolder,opt2name));
    
    %show the questions....
    %remember at this point, all the correct answers are in option1 and
    %gets switched around below
    if ismember(cellfun(@num2str, par.induct(i, 2), 'UniformOutput', false), ocean_comm) == 1
        env = 'ocean';
        envnum = 1;
        cuetext = 'This object above is found in the OCEAN.';
        if opt1pos == 1
            obj1 = obj1stim;
            obj1num = opt1num;
            obj1name = opt1name;
            obj2 = obj2stim;
            obj2num = opt2num;
            obj2name = opt2name;
        elseif opt1pos == 2
            obj1 = obj2stim;
            obj1num = opt2num;
            obj1name = opt2name;
            obj2 = obj1stim;
            obj2num = opt1num;
            obj2name = opt1name;
        end
    elseif ismember(cellfun(@num2str, par.induct(i, 2), 'UniformOutput', false), desert_comm) == 1
        env = 'desert';
        envnum = 2;
        cuetext = 'This object above is found in the DESERT.';
        if opt1pos == 1
            obj1 = obj1stim;
            obj1num = opt1num;
            obj1name = opt1name;
            obj2 = obj2stim;
            obj2num = opt2num;
            obj2name = opt2name;
        elseif opt1pos == 2
            obj1 = obj2stim;
            obj1num = opt2num;
            obj1name = opt2name;
            obj2 = obj1stim;
            obj2num = opt1num;
            obj2name = opt1name;
        end
    elseif ismember(cellfun(@num2str, par.induct(i, 2), 'UniformOutput', false), forest_comm) == 1
        env = 'forest';
        envnum = 3; 
        cuetext = 'This object above is found in the FOREST.';
        if opt1pos == 1
            obj1 = obj1stim;
            obj1num = opt1num;
            obj1name = opt1name;
            obj2 = obj2stim;
            obj2num = opt2num;
            obj2name = opt2name;
        elseif opt1pos == 2
            obj1 = obj2stim;
            obj1num = opt2num;
            obj1name = opt2name;
            obj2 = obj1stim;
            obj2num = opt1num;
            obj2name = opt1name;
        end
    end
   
    
    %Time for Cue, Statement, Option 1, and Option 2
    Screen(window, 'FillRect', par.bgcol, screenrect)
    Screen(window, 'TextSize', par.txtsize_instruct)
    DrawFormattedText(window, cuetext, 'center', 'center', par.txtcol);
    Screen('PutImage', window, cuestim, cuerect);
    Screen('PutImage', window, obj1, obj1rect);
    Screen('PutImage', window, obj2, obj2rect);
    Screen(window, 'Flip')
    
    %...and get the response
    resp = [];
    RT = 0;
    
    %start time for stimulus onset times
    respstart = GetSecs;
    endTime = respstart + pictime; 
    
    respwait = 1; 
    %start of response time
    while GetSecs < endTime 
        [keyIsDown, ~, keyCode] = KbCheck(-1); %check for key response
        if keyIsDown
            if keyCode(leftkey1)|| keyCode(leftkey2)
                respwait = 0; %#ok<NASGU>
                resp = 1;
                respstop = GetSecs;
            elseif keyCode(rightkey1)|| keyCode(rightkey2)
                respwait = 0; %#ok<NASGU>
                resp = 2;
                respstop = GetSecs;
            end      
        end %end of kb check
        if respwait == 0
            break;
        end        
    end   %end of while
    
    %CALCULATE ACCURACY DOWN HERE    
    if resp == 1
       if opt1pos == 1
          acc = 1;
          RT = respstop - respstart;
       elseif opt1pos == 2
          acc = 0;
          RT = respstop - respstart;
       end
    elseif resp == 2
       if opt1pos == 2
          acc = 1;
          RT = respstop - respstart;
       elseif opt1pos == 1
          acc = 0;
          RT = respstop - respstart;
       end 
    else
        resp = nan;
        RT = nan;
        acc = 0;
    end
    
    %% Save trial to matrix
    %Col1: SubjNum\t
    %Col2: TrialNum\t
    %Col3: QuestType\t
    %Col4: Environment\t
    %Col5: CueObject\t
    %Col6: Option1\t
    %Col7: Option2\t
    %Col8: Resp\t
    %Col9: Acc\t
    %Col10: RT\n
    thistrl = [header.subjNum, i, questtypenum, envnum, cuenum, obj1num, obj2num, resp, acc, RT];
    rundata = [rundata; thistrl]; %#ok<AGROW>
    
    %% Save results to file
    %Col1: SubjNum\t
    %Col2: TrialNum\t
    %Col3: QuestType\t
    %Col4: Environment\t
    %Col5: CueObject#\t
    %Col6: CueObjectName\t
    %Col7: Option1#\t
    %Col8: Option1Name\t
    %Col9: Option2#\t
    %Col10: Option2Name\t
    %Col11: Resp\t
    %Col12: Acc\t
    %Col13: RT\n
    fprintf(fid, '%d\t %d\t %s\t %s\t %d\t %s\t %d\t %s\t %d\t %s\t %d\t %d\t %f\n',...
        header.subjNum, i, questtype, env, cuenum, cuename, obj1num, obj1name, obj2num, obj2name, resp, acc, RT);
    
    %% save your matrix stuff
    data.rundata = rundata;
    
    %% output mat file for matrix stuff
    save(sprintf('%s/%s_InductGen',header.path.subjinfo,header.subjinfo),'data');
    
    %ITI b/w the questions
    Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
    Screen(window, 'TextSize', par.txtsize_fix);
    DrawFormattedText(window,'+','center','center',par.txtcol);
    Screen(window, 'Flip');
    WaitSecs(fixtime);
    
end

Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
Screen(window, 'TextSize', par.txtsize_instruct);
DrawFormattedText(window,'End of this phase!','center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(2);
sca;

close all
sca;