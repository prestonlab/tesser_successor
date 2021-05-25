function data = tesser_group(header)

Screen('Preference', 'SkipSyncTests', 1);

par = header.parameters;
data = struct('subjNum', header.subjNum);
data.type = 'Grouping';
data.header = header;
data.startTime = fix(clock);
addpath(header.path.stim);

rng(sum(100*clock)); %initialize the random number generator

%% get together all of the pictures
imageFolder = 'stimuli';
objs = par.picshuf;
objnum = num2cell([1:size(objs)]'); %#ok<NBRAK>
objs = [objs objnum];

%% put the pictures on an invisible grid first

%object and grid parameters
griddim_len = 19; %length of the grid
griddim_wid = 11; %width of the grid (will show up as this-1)
griddim = griddim_len * (griddim_wid-1);

% PSYCHTOOLBOX STUFF
par.txtcol = [0 0 0]; %black text
par.bgcol = [255 255 255]; %white background
par.progcol = [255 0  0]; %color of the progress bar
par.linecol = [200 200 200]; %black grid lines
par.txtsize = 75; %text size for almost everything
par.txtsize_test = 50; %text size for the test instructions
par.linesize = 1 ; %width of the grid lines in the rating task
par.progsize = 50; %height of the progress bar in pixels


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


Screen(window, 'TextFont', 'Arial');


%% CALCULATE MAXIMUM GRID SIZE
if screenWidth > screenHeight
    maxgrid = screenHeight;
    %maxgrid = screenWidth;
else
    maxgrid = screenWidth;
    %maxgrid = screenHeight;
end

%use maxgrid to figure out how big the grid spaces should be according to the display
gridsize = floor(maxgrid/griddim_wid-6); %divide the remaining space aka size of each grid
%this is taking the height(width) of the space, and dividing it by the number of width spaces designated
%since this space will be limited by the smallest dimension

%% IMAGE SIZING
%taking the divided number from above, and assigning the square sizes for
%the object dimensions too so that they fit within the squares
objheight = gridsize; %object height for study and test in pixels
objwidth = gridsize; %object width for study and test in pixels

%% MAKE THE GRID
linesout_len = griddim_len; %number of lines we need to go out in any direction for length
linesout_wid = griddim_wid; %number of lines we need to go out in any direction for width

%grid spacing
grid_len_dim = [gridsize/2 ((gridsize/2)+gridsize):gridsize:(((linesout_len/2)*gridsize)+gridsize/2)]; %positive values for length
grid_wid_dim = [gridsize/2 ((gridsize/2)+gridsize):gridsize:(((linesout_wid/2)*gridsize)+gridsize/2)]; %positive values for width
grid_len = [fliplr(grid_len_dim)*-1 grid_len_dim]; %add in the negative values for length
grid_wid = [fliplr(grid_wid_dim)*-1 grid_wid_dim]; %add in the negative values for width

%x direction (length)
xgrid = grid_len + xcenter;
xmin = min(xgrid);
xmax = max(xgrid);

%y direction (width)
ygrid = grid_wid  + ycenter;
ymin = min(ygrid);
ymax = max(ygrid);

%key press
contkey = KbName('7&');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up results file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Set up the output file
outfname = sprintf('%s/%s_Grouping',header.path.subjinfo,header.subjinfo);
if exist([outfname '.txt'],'file') == 2;
    error('The data file for this task already exists.  Please check your input parameters.');
end
fid=fopen([outfname '.txt'], 'a'); % open the file
fprintf(fid, 'SubjNum\t clicknum\t file\t objnum\t rt\n'); %Create header


%% START SCREEN
Screen(window, 'FillRect', par.bgcol);
Screen(window, 'TextSize', par.txtsize_instruct);
DrawFormattedText(window,'Press any key to start grouping!','center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(1);

%wait for button press to start the experiment
clear keyCode;
clear keyIsDown;
KbWait(-1);

%initial fixation
Screen(window, 'FillRect', par.bgcol, screenrect);
Screen(window, 'TextSize', par.txtsize_fix);
DrawFormattedText(window,'+','center','center',par.txtcol);
Screen(window, 'Flip');
WaitSecs(1);

%Create an empty 11x19 matrix
emptmat = zeros(11,19);

%figure out where to place the grid object
%the starting locations for all 21 objects
%[length, width]
gridobjall = [10 9;...
              18 4;...
              16 8;...
              1 6;...
              14 10;...
              16 3;...
              10 6;...
              4 9;...
              6 2;...
              19 6;...
              12 8;...
              2 2;...
              14 6;...
              18 10;...
              10 3;...
              6 6;...
              8 10;...
              4 4;...
              8 4;...
              2 8;...
              12 2;];

for i = 1:length(objs)
    
    gridobjpos_xlen = gridobjall(i, 1); %change back to not shuf if want the OG positions
    gridobjpos_ywid = gridobjall(i, 2); %change back to not shuf if want the OG poistions
    
    %positioning for the grid objects
    G1x1 = xgrid(gridobjpos_xlen); %the pixel vlaues of the x length lines
    xcor = gridobjpos_xlen; %the grid values of the y width lines
    G1y1 = ygrid(gridobjpos_ywid);
    ycor = gridobjpos_ywid;
    G1x2 = xgrid(gridobjpos_xlen+1);
    G1y2 = ygrid(gridobjpos_ywid+1);
    Grect = [G1x1 G1y1 G1x2 G1y2]; %making each square with the above and below x and y pixel lines
    
    %random grid/place assigment
    gridstim = objs{i};
    objid = objs{i, 2};
    emptmat(ycor, xcor) = i;
    startmat = emptmat;
    
    %load in stim
    [gridstimvalues, ~, ~] = imread(fullfile(imageFolder,gridstim));
    
    %show everything and get the response
    Screen('PutImage', window, gridstimvalues, Grect); %show the grid object within the grid
end

%draw the grid
for gridlines = 1:length(xgrid)
    Screen('DrawLine',window,par.linecol, xgrid(gridlines), ymin, xgrid(gridlines), ymax, par.linesize); %vertical lines
end

%draw the grid
for gridlines = 1:length(ygrid)
    Screen('DrawLine',window,par.linecol, xmin, ygrid(gridlines), xmax, ygrid(gridlines), par.linesize); %horizontal lines
end

%show the display
Screen(window, 'Flip');

%% Save out the graph as it looks initially
rundata = emptmat; %the first matrix with original positions
data.initgraph = rundata;
save(sprintf('%s/%s_GroupingInitial',header.path.subjinfo,header.subjinfo),'data');

%%
%Figure out what the object types are that are within these positions
%Community membership, primal or boundary node
%all of the items in a community
ocean_comm = par.comm1; %ocean community

desert_comm = par.comm2; %desert community

forest_comm = par.comm3; %forest community
%may be better to do this afterward

%%
%Now that all 21 items are in their orignal positions,
%figure out how to move them around, first to other unoccupied squares
%then figure out how to move them to the occupied squares

%response
xposmax = [];
xposmin = [];
yposmax = [];
yposmin = [];
respwait = 1;
respstart = GetSecs;

currobj = [];
currobjloc = [];
thismat = [];
thistrltotal = [];
respstart = GetSecs;
rtaccum = [];
clicknumtotal = 0;
respwait = 1;

while respwait == 1
    %wait for responses from the subject
    %[keyIsDown, ~, keyCode] = KbCheck; %check for key response
    [x, y, buttons] = GetMouse(window);  %check for mouse response
    
    [keyIsDown, ~, keyCode] = KbCheck(-1); %check for key response to stop break
    if keyIsDown
        if keyCode(contkey) %if press 7
            respwait = 0;
            %% Save out the graph as it looks finally
            data.finalgraph = startmat;
            save(sprintf('%s/%s_GroupingFinal',header.path.subjinfo,header.subjinfo),'data');
            
            %Then display the end screen
            Screen(window, 'FillRect', par.bgcol);
            Screen(window, 'TextSize', par.txtsize_instruct);
            DrawFormattedText(window,'End of this phase!','center','center',par.txtcol);
            Screen(window, 'Flip');
            WaitSecs(2);
            sca;
        end
    end
    
    %if they clicked the left mouse button, update the screen display
    if buttons(1)
        clicknum = 1;
        clicknumtotal = clicknumtotal + clicknum; %the number of 'trials'
        
        %figure out which grid the subject clicked
        xposmax = find(x <= xgrid,1,'first');
        xposmin = find(x >= xgrid,1,'last');
        yposmax = find(y <= ygrid,1,'first');
        yposmin = find(y >= ygrid,1,'last');
        %here this corresponds to the min and max of the top/bottom of the
        %square that makes up the actual grid position
        %the min x and y positions are needed just to be sure that the
        %participant isn't clicking outside of the grid
        
        %make sure they are actually clicking inside the grid
        if ~isempty(xposmax) && ~isempty(xposmin) && ~isempty(yposmax) && ~isempty(yposmin)
            %now that they have clicked inside the grid
            %we are now making sure that the place they clicked
            %contains an object
            
            
            if startmat(yposmin, xposmin) ~= 0 && isempty(currobj)
                currobj = startmat(yposmin, xposmin); %this is the object identity
                currobjloc = [yposmin, xposmin];  %this is the coordinate of where the object is
            elseif startmat(yposmin, xposmin) == 0 && ~isempty(currobj)
                startmat(yposmin, xposmin) = currobj;
                currobj = [];
                startmat(currobjloc(1), currobjloc(2)) = 0;
                currobjloc = [];
            end
            
            %% Save out the moves made
            respstop = GetSecs;            
            rt = respstop - respstart;  %#ok<NASGU>
                  
            
            if isempty(currobj) == 1
                objnum = 0;
                objname = '';
            else
                objnum = currobj;
                objname = objs{objnum};
            end
            
            %sequential mat file of all the clicks that were made
            thistrl = [header.subjNum, clicknumtotal, objnum, rt];
            thistrltotal = [thistrltotal; thistrl]; %#ok<AGROW>
            data.clicksgraph = thistrltotal;
            save(sprintf('%s/%s_GroupingClicks',header.path.subjinfo,header.subjinfo),'data');
            
            %sequential txt file of all the clicks that were made
            fprintf(fid, '%d\t %d\t %s\t %d\t %f\n',...
                header.subjNum, clicknumtotal, objname, objnum, rt);
            
            %here we need to redraw all of the objects again
            %basically redraw it referring to the startmat
            
            %show the objects
            for xlen = 1:length(startmat(1, : ))
                for ywid = 1:length(startmat(:, 1))
                    if startmat(ywid, xlen) ~= 0
                        %this is the object ID
                        objid = startmat(ywid, xlen); %current object number
                        gridstimput = objs{objid};
                        [placestimvalues, ~, ~] = imread(fullfile(imageFolder,gridstimput));
                        Screen('PutImage', window, placestimvalues, [xgrid(xlen) ygrid(ywid) xgrid(xlen+1) ygrid(ywid+1)])  %now the place object in the grid
                    end
                end
            end
            
            %now draw the grid again
            
            %draw the grid horz
            for gridlines = 1:length(xgrid)
                Screen('DrawLine',window,par.linecol, xgrid(gridlines), ymin, xgrid(gridlines), ymax, par.linesize); %vertical lines
            end
            
            %draw the grid vert
            for gridlines = 1:length(ygrid)
                Screen('DrawLine',window,par.linecol, xmin, ygrid(gridlines), xmax, ygrid(gridlines), par.linesize); %horizontal lines
            end
            
            %% Save out the .mat for graph of moves made
            thismat = cat(3, thismat, startmat);
            data.placegraph = thismat;
            save(sprintf('%s/%s_GroupingPlaces',header.path.subjinfo,header.subjinfo),'data');
            
            
            %show the display
            Screen(window, 'Flip');
            
        end %end grid object/place object overlap check
        
    end %end inside grid check
    
end %end of mouse button check
end
