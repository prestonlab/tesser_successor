function data = tesser_parse(header)

Screen('Preference', 'SkipSyncTests', 1);

par = header.parameters;
data = struct('subjNum', header.subjNum);
data.type = 'StructParse';
data.header = header;
data.startTime = fix(clock);
addpath(header.path.stim);

rng(sum(100*clock)); %initialize the random number generator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Keyboard setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET RESPONSE OPTIONS
yeskey = KbName('space'); %button if to parse
contkey = KbName('7&');

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
outfname = sprintf('%s/%s_StructParse',header.path.subjinfo,header.subjinfo);
if exist([outfname '.txt'],'file') == 2;
    error('The data file for this task already exists.  Please check your input parameters.');
end
fid=fopen([outfname '.txt'], 'a'); % open the file
fprintf(fid, 'SubjNum\t run\t trial\t objnum\t file\t objseq\t objtype\t respnum\t resp\t acc\t rt\t parsetype\n'); %Create header

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run experiment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START SCREEN
Screen(window, 'FillRect', par.bgcol);
Screen(window, 'TextSize', par.txtsize_instruct);
DrawFormattedText(window,'Press the spacebar at times in the sequence\nwhich feel like natural breaking points.\n\nIf you''re unsure, just try to use your intuition.\nPress the spacebar to begin!',...
    'center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(1);


%wait for button press to start the experiment
clear keyCode;
clear keyIsDown;
KbWait(-1);

%start time for stimulus onset times
imageFolder = 'stimuli';

objs = par.picshuf;
seq = par.sequence_parse;

%dividing the sequence into segments
seqsize = length(par.sequence_parse);
numbreaks = 3;

seqbreak = seqsize/numbreaks;
seqarray = seqbreak*(1:numbreaks);
seqarray = [0 seqarray];
seqdiv = [];

for i = 1:numbreaks
    this = seqarray(i)+1;
    seqeach = seq(:, this:seqarray(i+1))';
    seqdiv = horzcat(seqdiv, seqeach);
end

rundata = [];

%% Start sequence
for s = 1:numbreaks %1:3
    
    %Run experimental trials
    for t = 1:size(seqdiv, 1) %remember this is the number of items in list when divided by the breaks
        
        %% Stimuli details and coordinates
        
        %Load image
        objnum = seqdiv(t, (s*3)-2); %pick out the t image number from the sequence
        objseq = seqdiv(t, (s*3)-1); %pick out if the item part of random sequence [1] or Hamiltonian Forward [2] or Hamiltonian Backward [3]
        objtype = seqdiv(t, (s*3));  %pick out if the item boundary [1] or nonboundary [0]
        
        objname = objs{objnum}; %pick out that number of item from the list of shuff object names
        
        img = imread(fullfile(imageFolder,objname)); %read out that image
        imageDisplay = Screen('MakeTexture', window, img); %display said image
        
        %Stimulus display time
        maxpictime = par.parseTimeout;
        
        %% Actually presenting the stimuli
        pictime = GetSecs; %when the item should begin
        Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
        Screen('DrawTexture', window, imageDisplay); %putting the image on the screen with specified dimensions
        Screen('Flip', window); % Start of trial
        
        %Get keypress response
        RT = 0;
        resp = [];
                
        %response
        respstart = GetSecs; %start of response time
        while GetSecs < pictime + maxpictime
            [keyIsDown, ~, keyCode] = KbCheck(-1); %check for key response
            if keyIsDown
                if keyCode(yeskey)
                    resp = 'PARSED';
                    respnum = 1;
                    respstop = GetSecs;
                    RT = respstop - respstart;
                end
            end %end of keycode check
        end  %end of kbcheck
        
        if isempty(resp)
            resp = 'NONE';
            respnum = 0;
            RT = NaN;
        end
        
        %CALCULATE ACCURACY DOWN HERE
        if respnum == 1
            if objseq == 1 && objtype == 1
                acc = 1;
                parsetype = 'Bound_Ran';
                parsetypenum = 11;
            elseif objseq == 1 && objtype == 0
                acc = 0;
                parsetype = 'Non_Ran';
                parsetypenum = 21;
            elseif objseq == 2 && objtype == 1
                acc = 1;
                parsetype = 'Bound_ForH';
                parsetypenum = 12;
            elseif objseq == 2 && objtype == 0
                acc = 0;
                parsetype = 'Non_ForH';
                parsetypenum = 22;
            elseif objseq == 3 && objtype == 1
                acc = 1;
                parsetype = 'Bound_BckH';
                parsetypenum = 13;
            elseif objseq == 3 && objtype == 0
                acc = 0;
                parsetype = 'Non_BckH';
                parsetypenum = 23;
            end
        elseif respnum == 0
            acc = NaN;
            parsetype = 'NONE';
            parsetypenum = NaN;
        end
        
        %% Save trial into matrix
        %Col1: Subj#
        %Col2: Run#
        %Col3: Trial#
        %Col4: Object#
        %Col5: Object sequence number: 1 (random), 2 (hamfor), 3 (hamback)
        %Col6: Object type: 1 (boundary), 0 (nonboundary)
        %Col7: Response: 1 (parsed) or 0 (no parse)
        %Col8: Parse accuracy (1 or 0)
        %Col9: Parse RT
        %Col10: Parse Type Descriptor
        thistrl = [header.subjNum, s, t, objnum, objseq, objtype, respnum, acc, RT, parsetypenum];
        rundata = [rundata; thistrl]; %#ok<AGROW>
        
        %% Save results to text file
        %Col1: Subj#
        %Col2: Run#
        %Col3: Trial#
        %Col4: Object#
        %Col5: Object File
        %Col6: Object sequence number: 1 (random), 2 (hamfor), 3 (hamback)
        %Col7: Object type: 1 (boundary), 0 (nonboundary)
        %Col8: Response: PARSED or NaN
        %Col9: Parse accuracy (1 or 0)
        %Col10: Parse RT
        %Col11: Parse Type Descriptor
        fprintf(fid, '%d\t %d\t %d\t %d\t %s\t %d\t %d\t %d\t %s\t %d\t %f\t %s\n',...
            header.subjNum, s, t, objnum, objname, objseq, objtype, respnum, resp, acc, RT, parsetype);
        
        %% save your matrix stuff
        data.rundata = rundata;
        
        %% output mat file for matrix stuff
        save(sprintf('%s/%s_StructParse',header.path.subjinfo,header.subjinfo),'data');
    end
    %end of trials
    
    %adding a break between trials
    %break between runs
    if s < numbreaks %if between runs, break
        clear keyCode; %clearing any prior key presses
        clear keyIsDown;
        
        Screen(window, 'FillRect', par.bgcol)
        Screen(window, 'TextSize', par.txtsize_instruct);
        disptext = sprintf('BREAK TIME!\n\nPress 7 when you are ready to continue.'); %text to show on the screen
        DrawFormattedText(window,disptext,'center','center',par.txtcol); %show the time left
        Screen(window, 'Flip');
        
        r = 1;
        while r == 1
            [keyIsDown, ~, keyCode] = KbCheck(-1); %check for key response to stop break
            if keyIsDown
                if keyCode(contkey) %if press spacebar
                    break
                end
            end
        end
        
        %final second
        WaitSecs(1);
        
    else  %if it's the final run, go to end screen
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% End the experiment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Screen(window, 'FillRect', par.bgcol); %creating a blank screen for the images
        Screen(window, 'TextSize', par.txtsize_instruct);
        DrawFormattedText(window,'End of this phase!','center','center',par.txtcol);
        Screen(window, 'Flip');
        WaitSecs(2);
        sca;
        
    end
    
end

close all
sca;