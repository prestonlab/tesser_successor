% Athula Pudhiyidath, June 2017
% header file for Tesser Scanning
% creates a subject data folder, mat file, and workspace structure called 'header' that is read into sub functions

% ##input##
% none immediately, but will prompt:
% subject number
% subject initials
% subject age
% subject gender
%
% ##output##
% data directory, TesserScan/data/
% header file, data/tesser_xxx_yy_header.mat and also tesser_xxx_yy_header.txt
% workspace structure, header, that contains all of the experiment parameters for subfunctions
%
% ##current design (as of 06/13/17)##
%---PART 1:---
%---orientation learning---
% participants view all the objects in the canonical orientation
%
%
%---PART 2:--- [A and B]
%---community structure learning---
%--[A]: 
% the community structure is learned through only random walks
% participants will press the spacebar whenever they see an item upside down
% 1575 trials 
% 1.5s for each trial, no ITI
% 5 runs to divide this phase (7.875min each)
%
%--[B]:
% the community structure is continued to be learned
% this time through true-random and random walks
% orientation task performed again
% 882 trials
% 1s for each trial
% Random jitter of 1s, 3s, 5s 
% 6 runs to divide this phase (9.8min each)
%
%
%---PART 3:---
%---inductive generalization questions---
% questions about the community structure
% 42 total questions with 2 AFC 
% 
%
%---PART 4:---
%---parsing task---
% participants will now press the spacebar at changing intervals
% 756 trials (aka 36 blocks: 18 random blocks, 18 hamfor and hambak)
% 1.5s for each trial 
% 3  breaks to divide this phase (6.3min each)
%
%
%%---PART 5:---
%---grouping task---
% participants will now group the stimuli in accordance to groups they see
% 21 items on the screen at once randomly dispersed
% self-paced

%% clear any preexisting stuff in the workspace and command window
clear all;
clc;

%% SETUP
header = struct('exp', 'tesserScan','version','May 2017','parameters',struct); %experiment info
header.path.exp = pwd; %path for experiment files
par = struct;
header.path.data = [header.path.exp '/data/']; %path for data folder
header.path.stim = [header.path.exp '/stimuli/']; %path for stimuli
header.setuptime = fix(clock); %time this was run

%% SUBJECT INFO
subjNum = input('Subject Number:  '); %subject number
subjName = input('Subject Initials:  ', 's'); %subject initials
subjAge= input('Subject Age:  '); %subject initials
subjGen= input('F(1) or M(0)?:  '); %subject initials

header.subjinfo = ['tesserScan_', num2str(subjNum), '_', subjName]; %add subject info to header
header.path.subjinfo = [header.path.data sprintf('%s', header.subjinfo)];
header.subjNum = subjNum;
header.subjAge = subjAge;
header.subjGen = subjGen; 
rng(subjNum); %seed the random number generator by subject number

%check whether a directory with that name exists already
if exist(header.path.subjinfo)~=0; %#ok<EXIST>
    disp('!!! FILE for this subject exists already! Aborting...');
    clear all; return;  %abort experiment, rather than overwrite existing data
end

%if we're good, make the data folder
mkdir(sprintf('%s/%s',header.path.data, header.subjinfo));

%% THINGS YOU MIGHT WANT TO EDIT
% PSYCHTOOLBOX STUFF
par.txtcol = [0 0 0]; %black text
par.bgcol = [255 255 255]; %white background
par.progcol = [255 0  0]; %color of the progress bar
par.linecol = [0 0 0]; %black grid lines
par.txtsize_instruct = 50; %text size for almost everything
par.txtsize_fix = 75; 
par.linesize = 2; %width of the grid lines in the rating task
par.progsize = 50; %height of the progress bar in pixels

%% OBJECTS
%now I will import list of novel object names
par.pics = {'object_1.jpg'; 'object_2.jpg'; 'object_3.jpg'; 'object_4.jpg'; 'object_6.jpg'; 'object_7.jpg'; 'object_9.jpg';...
            'object_10.jpg'; 'object_12.jpg'; 'object_13.jpg'; 'object_16.jpg'; 'object_17.jpg'; 'object_18.jpg'; 'object_23.jpg';...
            'object_24.jpg'; 'object_25.jpg'; 'object_28.jpg'; 'object_30.jpg'; 'object_34.jpg'; 'object_35.jpg'; 'object_36.jpg';};
%this will be used when participants are making the similarity ratings

par.picshuf = par.pics(randperm(size(par.pics, 1)), :);
%this will be used when participants are looking at community structure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------PART 0: Image Setup---------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMAGE SIZING
par.objheight = 400; %object height for study and test in pixels
par.objwidth = 400; %object width for study and test in pixels
par.buffer = 10; %display buffer for the similarity ratings in pixels
par.griddim = 7; %size of the grid for the similarity ratings (griddim x griddim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------PART 2: Community Structure-------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%----PART 2: A-----%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%THESE ARE THE NON-SCANNED RUNS%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Community structure generation
%made only of random walks

rng(sum(100*clock)); %initialize the random number generator

par.seqLength_Part2a = 1575;
subtest_length_a = 21; %the length of a total walkthrough
numwalks_a = par.seqLength_Part2a/(subtest_length_a); %the number of walks
curr_a = Randi(21); %current position is a random number b/w 1:21
sequence_a = NaN(0,numwalks_a); %creating an empty matrix of length numwalks
index_a = 1;

comm1 = [1,2,3,18,19,20,21]; %OCEAN
comm2 = [4,5,6,7,8,9,10];   %DESERT
comm3 = [11,12,13,14,15,16,17]; %FOREST

graph_a = [2,3,18,19,20,21; 3,18,19,20,21,1; 4,19,20,21,1,2; 5,6,7,8,9,3; 6,7,8,9,10,4; 7,8,9,10,4,5;...
    8,9,10,4,5,6; 9,10,4,5,6,7; 10,4,5,6,7,8; 5,6,7,8,9,11; 12,13,14,15,16,10; 13,14,15,16,17,11;...
    14,15,16,17,11,12; 15,16,17,11,12,13; 16,17,11,12,13,14; 17,11,12,13,14,15; 12,13,14,15,16,18;...
    19,20,21,1,2,17; 20,21,1,2,3,18; 21,1,2,3,18,19; 1,2,3,18,19,20];
%the graph itself


for j = 1:numwalks_a %60
    % generating a random walk for 21 items
    for i = 1:subtest_length_a %subtest_length %21
        next_a = Randi(6); %random between 1 and 6
        curr_a = graph_a(curr_a,next_a); %changes curr index to going from one of 1/21 and then goes to the
        %column of the random of 1/6
        sequence_a(index_a) = curr_a; %sequence starts at 1 and goes down that empty sequence matrix
        index_a = index_a+1; %the index is that + 1
    end
end

%disp(sequence); %uncomment later if you want to see
%disp(reshape(sequence,21,par.total_length/21)'); %uncomment later if you want to see

nbackmarkers_a = [];
nbseq_a = sequence_a;

%Now add in repeat markers for 2-back test
for g = 1:length(nbseq_a)
    if g == 2
        comp1_a = nbseq_a(1, g-1);
        comp2_a = nbseq_a(1, g);
        
        if comp1_a == comp2_a
            nbackmarkers_a(1, 1) = 1;
            nbackmarkers_a(1, 2) = 1;
        else
            nbackmarkers_a(1, 1) = 0;
            nbackmarkers_a(1, 2) = 0;
        end
        
    elseif g > 2
        comp1_a = nbseq_a(1, g-2);
        comp2_a = nbseq_a(1, g);
        
        if comp1_a == comp2_a
            nbackmark_a = 1;
            nbackmarkers_a = [nbackmarkers_a nbackmark_a];
        else
            nbackmark_a = 0;
            nbackmarkers_a = [nbackmarkers_a nbackmark_a];
        end
        
    end
end

seqall_a = [nbseq_a; nbackmarkers_a];

%%%%%%%%%%%%%%%%
par.Seq_Part2a = seqall_a;
%a 2x1218 sequence
%row 1: the sequence (all random walks)
%row 2: the nback marks
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%now divide this into the 6 runs:
%%%%%%%%%%%%%%%%

seq_a = par.Seq_Part2a; 
seqsize_a =  length(par.Seq_Part2a);
numbreaks_a = 5;

seqbreak_a = seqsize_a/numbreaks_a;
seqarray_a = seqbreak_a*(1:numbreaks_a);
seqarray_a = [0 seqarray_a];

%now we're taking the length of the sequence and dividing it across
%columns and we will be going through all of the columns to read the seq
seqdiv_a = [];
for i = 1:numbreaks_a
    this = seqarray_a(i)+1;
    seqeach_a = seq_a(:, this:seqarray_a(i+1))';
    seqdiv_a{i} = seqeach_a;
end

%%%%%%%%%%%%%%%%
%now a 3x151 sequence, with four null trials at beg/end for 6 runs
%in each of the 6 runs: 
%row 1: Sequence | 0 = null trial 
%row 2: 1 = Random Walk 2 = True Random | 0 = null trial
%row 3: nback marker | NaN = null trial
par.Seq_Part2a_div = seqdiv_a; 
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 2B TIMING %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%THESE ARE THE NON-SCANNED TIMINGS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par.structTimeout = 1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
%----PART 2: B-----%
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%THESE ARE THE SCANNED RUNS%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%amending parsing task to fit the scanning version
%will alternate between random walks and true randoms
%will alternate between these on a set schedule

%882 total length (42 walks)

par.seqLength_Part2b = 882;
subtest_length_scan = 21;
length_scan = par.seqLength_Part2b/(subtest_length_scan);
curr_scan = Randi(21); %current position is a random number b/w 1:21
sequence_scan = []; %NaN(0,length_scan); %creating an empty matrix of length 42
sequence_scan_all = [];
index_scan = 1; %1


% community 1 = 1,2,3,18,19,20,21; %ocean
% community 2 = 4,5,6,7,8,9,10; %desert
% community 3 = 11,12,13,14,15,16,17; %forest

boundarynodes_scan = [18,3,4,10,11,17];

graph_scan = [2,3,18,19,20,21; 3,18,19,20,21,1; 4,19,20,21,1,2; 5,6,7,8,9,3; 6,7,8,9,10,4; 7,8,9,10,4,5;...
    8,9,10,4,5,6; 9,10,4,5,6,7; 10,4,5,6,7,8; 5,6,7,8,9,11; 12,13,14,15,16,10; 13,14,15,16,17,11;...
    14,15,16,17,11,12; 15,16,17,11,12,13; 16,17,11,12,13,14; 17,11,12,13,14,15; 12,13,14,15,16,18;...
    19,20,21,1,2,17; 20,21,1,2,3,18; 21,1,2,3,18,19; 1,2,3,18,19,20];
%the graph itself

truerand_runs = [1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 15, 16, 18, 19, 21, 22, 24, 25, 27, 28, 29, 30, 32, 33, 35, 36, 38, 39, 41, 42];
randonly_runs = [3, 6, 9, 12, 17, 20, 23, 26, 31, 34, 37, 40];

for j = 1:length_scan %42
    %true random walks
    truerand = randperm(21, 21);
    
    % now have to divide the walks depending on truerand or randonly
    if ismember(j, randonly_runs) == 0
        % generating a random walk for 21 items
        for i = 1:subtest_length_scan
            next_scan = Randi(6); %random between 1 and 6
            curr_scan = graph_scan(curr_scan,next_scan); %changes curr index to going from one of 1/21 and then goes to the
            %column of the random of 1/4
            sequence_scan(1, index_scan) = curr_scan; %sequence starts at 1 and goes down that empty sequence matrix: 1st ROW
            sequence_scan(2, index_scan) = 1; %says that this is the normal sequence path via [1]; 2nd ROW
            index_scan = index_scan+1; %the index is that + 1
        end
    else
        for i = 1:subtest_length_scan
            truerandmemb = truerand(i); 
            %else add in a true random sequence
            sequence_scan(1, index_scan) = truerandmemb;
            curr_scan = truerand(max(truerand)); %last index of the truerand is where the randonly will begin again
            sequence_scan(2, index_scan) = 2;
            index_scan = index_scan+1;
        end
    end
    
    % fitting the two walks together
    sequence_scan_all = sequence_scan; %#ok<*AGROW>
    
    %begin with the true random walk truly randomly (mostly)
end


nbackmarkers_scan = [];
nbseq_scan = sequence_scan_all;

%Now add in repeat markers for 2-back test
for g = 1:length(nbseq_scan)
    if g == 2
        comp1_scan = nbseq_scan(1, g-1);
        comp2_scan = nbseq_scan(1, g);
        
        if comp1_scan == comp2_scan
            nbackmarkers_scan(1, 1) = 1;
            nbackmarkers_scan(1, 2) = 1;
        else
            nbackmarkers_scan(1, 1) = 0;
            nbackmarkers_scan(1, 2) = 0;
        end
        
    elseif g > 2
        comp1_scan = nbseq_scan(1, g-2);
        comp2_scan = nbseq_scan(1, g);
        
        if comp1_scan == comp2_scan
            nbackmark_scan = 1;
            nbackmarkers_scan = [nbackmarkers_scan nbackmark_scan];
        else
            nbackmark_scan = 0;
            nbackmarkers_scan = [nbackmarkers_scan nbackmark_scan];
        end
        
    end
end

seqall_withnback_scan = [sequence_scan_all; nbackmarkers_scan];

%%%%%%%%%%%%%%%%
%now a 3x892 sequence
%this will be 147 trials in each of the 6 runs 
%row 1: Sequence 
%row 2: 1 = Random Walk 2 = True Random 
%row 3: nback marker 
par.Seq_Part2b = seqall_withnback_scan;
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%now divide this into the 6 runs:
%%%%%%%%%%%%%%%%
seq_scan = par.Seq_Part2b; 
seqsize_scan =  length(par.Seq_Part2b);
numbreaks_scan = 6;

seqbreak_scan = seqsize_scan/numbreaks_scan;
seqarray_scan = seqbreak_scan*(1:numbreaks_scan);
seqarray_scan = [0 seqarray_scan];

%now we're taking the length of the sequence and dividing it across
%columns and we will be going through all of the columns to read the seq
seqdiv_scan = [];
for i = 1:numbreaks_scan
    this = seqarray_scan(i)+1;
    seqeach_scan = seq_scan(:, this:seqarray_scan(i+1))';
    %Now, add two null trials at the beginning and end of the sequence
    %these will be for the buffer TRs that are seen at the beg and end
    null_list = [0 0 NaN; 0 0 NaN]; 
    seqeach_scan = [null_list; seqeach_scan; null_list]; 
    seqdiv_scan{i} = seqeach_scan;
end

%%%%%%%%%%%%%%%%
%now a 3x151 sequence, with four null trials at beg/end for 6 runs
%in each of the 6 runs: 
%row 1: Sequence | 0 = null trial 
%row 2: 1 = Random Walk 2 = True Random | 0 = null trial
%row 3: nback marker | NaN = null trial
par.Seq_Part2b_div = seqdiv_scan; 
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 2B TIMING %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%THESE ARE THE SCANNED TIMINGS%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Community structure onset time + ISI time
%figure out one run, then shuffle for the rest of them: 

par.pictime_Part2b = 1; %max time for an item showing
trialtime_all = repmat(par.pictime_Part2b, par.seqLength_Part2b/6, 1); %for one run 
%add in 0 trials for beg and end

trialtime_all = [0; 0; trialtime_all; 0; 0];
par.pictime_all_Part2b = trialtime_all; 

%jittered ITI between: 
%49 1s ITIs
sec1 = repmat(1, 49, 1);

%49 3s ITIs
sec3 = repmat(3, 49, 1);

%49 5s ITIs
sec5 = repmat(5, 49, 1);

allsecs = vertcat(sec1, sec3, sec5); 
toshuffle = allsecs(2:end);

alltogjit = [];
for m = 1:6
    thisshuff = Shuffle(toshuffle);
    finalshuff = [1; thisshuff]; 
    nulltrials_beg = [0; 2]; 
    nulltrials_end = [2; 2]; 
    thisshufftog = [nulltrials_beg; finalshuff; nulltrials_end]; 
    alltogjit = [alltogjit thisshufftog];
end

par.jittertime_all_Part2b = alltogjit;
%for all 6 runs together in each column, with null trials at beg + end

%add together trial and the jitter
onsetrun = [];
for p = 1:6
    runjit = par.jittertime_all_Part2b(:, p);
    onsettogall = [];
    for t = 1:length(runjit)        
        onsettog = par.pictime_all_Part2b(t) + runjit(t);
        onsettogall = [onsettogall; onsettog];
    end
    
    onsetrun = [onsetrun onsettogall]; 
end

%onsetrun gives the trial length in pure seconds with 1s stim time + 1,3,5 ITI
%each run is 147 items long

%now added two 2s null trials at the beginning and end of each of the 6 runs 
%will now be 151 onset times long

par.trialjittertime_all_Part2b  = onsetrun; 

%now adding the onsets together 
scantimetog = []; 
for p = 1:6
    addjit = par.trialjittertime_all_Part2b(:, p); 
    scantime = [];
    initial = 0;
    
    for o = 1:length(addjit) 
        newtime = initial + addjit(o);
        initial = newtime;   
        scantime = [scantime; initial];
    end
    
    scantimetog = [scantimetog scantime];
end

par.onsettime_all_Part2b = scantimetog; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------PART 3: Inductive Generalization--------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This is how the community structure is, well... structured
%Each of the communities as a whole, their primary nodes, and their boundary nodes
par.comm1 = [1,2,3,18,19,20,21]; %OCEAN
par.comm1pri = [1,2,19,20,21];
par.comm1bou = [18,3];

par.comm2 = [4,5,6,7,8,9,10]; %DESERT
par.comm2pri = [5,6,7,8,9];
par.comm2bou = [4,10];

par.comm3 = [11,12,13,14,15,16,17]; %FOREST
par.comm3pri = [12,13,14,15,16];
par.comm3bou = [11,17];

%For the induction generalization questions, we will be testing the primary
%nodes, the 1-away boundary nodes, and the 2-away boundary nodes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRIMARY NODE CUES + PRIMARY NODE CHOICES [The 1st choice is correct here!!]
%CUE, OPTION 1, OPTION 2, ANSWER SWAP

%% Question Type, Option 1 %From Excel file
prim = [1 2; 2 19; 19 20; 20 21; 21 1; 21 2; 19 21; 20 1; 2 20; 1 19;...
        5 6; 6 7; 7 8; 8 9; 9 5; 9 6; 7 9; 8 5; 6 8; 5 7;...
        12 13; 13 14; 14 15; 15 16; 16 12; 13 16; 16 14; 12 15; 15 13; 14 12;];
    
%add together the second choice for each by shuffling communities
%% Option 2 %Shuffling and putting together Option 2 (in accordance with Excel sheet)
primcom1shuff = transpose(par.comm1pri((randperm(length(par.comm1pri)))));
primcom2shuff = transpose(par.comm2pri((randperm(length(par.comm2pri)))));
primcom3shuff = transpose(par.comm3pri((randperm(length(par.comm3pri)))));

%these correspond to the red/green/purple communities that are the incorrect answers on the Excel sheet
primchoicetogether = [primcom2shuff; primcom3shuff; primcom1shuff; primcom3shuff; primcom1shuff; primcom2shuff]; 

%% Question Type, Option 1, Option 2
prim_all = [prim primchoicetogether];

%now need to add a column for the placement of the correct choice (always
%column 1) in either the first position or the second position

%so this would also correspond to the correct response they're supposed to be pressing (1 or 2)
pr1 = repmat(1, length(prim_all)/2, 1); %#ok<RPMT1> %half of the choices are '1'
pr2 = repmat(2, length(prim_all)/2, 1); %half of the choices are '2'
ptog = [pr1; pr2];
ptogshuff = randperm(length(ptog(:, 1)));
ptog = ptog(ptogshuff, :); %shuffled together the choices

primtogether = [prim_all ptog]; %adding the cue/opt1/opt2 with answer
primtogshuff = randperm(length(primtogether(:, 1)))'; 
primfinal =  primtogether(primtogshuff, :); %shuffling together all of the cue things

primname = repmat({'Prim'}, length(primfinal), 1); %adding name of 'Prim' to ID the type
primfinalcell = num2cell(primfinal); 

primfinalized = [primname, primfinalcell];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1-APART BOUNDARY NODE CUES + PRIMARY AND BOUNDARY CHOICES [The 1st choice is correct here!!]
%CUE, OPTION 1, OPTION 2, ANSWER SWAP

%% Question Type %From Excel file
bound1cue = [3; 18; 4; 10; 11; 17]; %From Excel file

%% Option 1 %From Excel file, corresponding to cue 
pricom1rand = randperm(length(par.comm1pri(1, :)));
pricom1 = par.comm1pri(:, pricom1rand); %randomly shuffled PriComm1

pricom2rand = randperm(length(par.comm2pri(1, :)));
pricom2 = par.comm2pri(:, pricom2rand); %randomly shuffled PriComm2

pricom3rand = randperm(length(par.comm3pri(1, :)));
pricom3 = par.comm3pri(:, pricom3rand); %randomly shuffled PriComm3

%taking first and second index of these randomly shuff primary community members 
bound1opt1 = [pricom1(1); pricom1(2); pricom2(1); pricom2(2); pricom3(1); pricom3(2)];

%% Option 2 %From Excel file, corresponding to cue 
bound1opt2 = [4; 17; 3; 11; 10; 18];

%% Question Type, Option 1, Option 2
bound1_all = [bound1cue bound1opt1 bound1opt2];

%now need to add a column for the placement of the correct choice (always column 1) 
%in either the first position or the second position

%so this would also correspond to the correct response they're supposed to be pressing (1 or 2)
b11 = repmat(1, length(bound1_all)/2, 1); %#ok<RPMT1> %half the choices are '1'
b12 = repmat(2, length(bound1_all)/2, 1); %half the choices are '2'
b1tog = [b11; b12];
b1shuff = randperm(length(b1tog(:, 1)));
b1tog = b1tog(b1shuff, :); %shuffled list of answer choices together

bound1i = [bound1_all, b1tog]; %putting together cue/opt1/opt2 with answers
bound1ishuff = randperm(length(bound1i(:, 1)))';
bound1final = bound1i(bound1ishuff, :); %shuffling together all of the cue things 

bound1name = repmat({'Bound1'}, length(bound1final), 1); %adding name of 'Bound1' to ID the type
bound1finalcell = num2cell(bound1final);

bound1finalized = [bound1name, bound1finalcell];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2-APART BOUNDARY NODE CUES + PRIMARY AND BOUNDARY CHOICES [The 1st choice is correct here!!]
%CUE, OPTION 1, OPTION 2, ANSWER SWAP

%% Question Type, Option 1 %From Excel file
bound2cue = [3 18; 4 10; 11 17; 18 3; 10 4; 17 11]; %From Excel file

%% Option 2 %Shuffling and putting together Option 2 (in accordance with Excel sheet)
pricom1rand = randperm(length(par.comm1pri(1, :)));
pricom1 = par.comm1pri(:, pricom1rand); %randomly shuffled PriComm1

pricom2rand = randperm(length(par.comm2pri(1, :)));
pricom2 = par.comm2pri(:, pricom2rand); %randomly shuffled PriComm2

pricom3rand = randperm(length(par.comm3pri(1, :)));
pricom3 = par.comm3pri(:, pricom3rand); %randomly shuffled PriComm3

%taking first and second index of these randomly shuff primary community members 
bound2opt2 = [pricom2(1); pricom1(1); pricom2(2); pricom3(1); pricom3(2); pricom1(2)];

%% Question Type, Option 1, Option 2
bound2_all = [bound2cue bound2opt2];

%now need to add a column for the placement of the correct choice (always
%column 1) in either the first position or the second position

%so this would also correspond to the correct response they're supposed to be pressing (1 or 2)
b21 = repmat(1, length(bound2_all)/2, 1); %#ok<RPMT1>
b22 = repmat(2, length(bound2_all)/2, 1);
b2tog = [b21; b22];
b2shuff = randperm(length(b2tog(:, 1)));
b2tog = b2tog(b2shuff, :);

bound2i = [bound2_all, b2tog]; %putting together cue/opt1/opt2 with answers
bound2ishuff = randperm(length(bound2i(:, 1)))';
bound2final = bound2i(bound2ishuff, :); %shuffling together all of the cue things 

bound2name = repmat({'Bound2'}, length(bound2final), 1); %adding name of 'Bound2' to ID the type
bound2finalcell = num2cell(bound2final);

bound2finalized = [bound2name, bound2finalcell];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%All of the induction questions together
induc = [primfinalized; bound1finalized; bound2finalized];

%Now we need to shuffle around all of these induction questions around and
%not in the chunks of primary/boundary etc.

inducshuff = randperm(length(induc(:, 1)));
inducfinal = induc(inducshuff, :); 

%now that you have a shuffled list of all of the inductive questions
%make sure that the questions have no more than 2 cues and questions
%backtoback from the same community
testcomm1 = {'1', '2', '3', '18', '19', '20', '21'}; %ocean
testcomm2 = {'4', '5', '6', '7', '8', '9', '10'}; %desert
testcomm3 = {'11', '12', '13', '14', '15', '16', '17'}; %forest

cuecomm = [];
cuecommy = [];

opt1comm = [];
opt1commy = [];

opt2comm = [];
opt2commy = [];

nodetype = [];
nodetypey = [];

%assigning the correct community for the Cues Numerical Type
for a = 1:length(inducfinal)
    if ismember(cellfun(@num2str, inducfinal(a, 2), 'UniformOutput', false), testcomm1) == 1
        cuecommy = 1;
    elseif ismember(cellfun(@num2str, inducfinal(a, 2), 'UniformOutput', false), testcomm2) == 1
        cuecommy = 2;
    elseif ismember(cellfun(@num2str, inducfinal(a, 2), 'UniformOutput', false), testcomm3) == 1
        cuecommy = 3;
    end
    cuecomm = vertcat(cuecomm, cuecommy); %#ok<AGROW>
end

%assigning the correct community for Option 1 Numerical Type
for b = 1:length(inducfinal)
    if ismember(cellfun(@num2str, inducfinal(b, 3), 'UniformOutput', false), testcomm1) == 1
        opt1commy = 1;
    elseif ismember(cellfun(@num2str, inducfinal(b, 3), 'UniformOutput', false), testcomm2) == 1
        opt1commy = 2;
    elseif ismember(cellfun(@num2str, inducfinal(b, 3), 'UniformOutput', false), testcomm3) == 1
        opt1commy = 3;
    end
    opt1comm = vertcat(opt1comm, opt1commy); %#ok<AGROW>
end

%assigning the correct community for Option 2 Numerical Type
for c = 1:length(inducfinal)
    if ismember(cellfun(@num2str, inducfinal(c, 4), 'UniformOutput', false), testcomm1) == 1
        opt2commy = 1;
    elseif ismember(cellfun(@num2str, inducfinal(c, 4), 'UniformOutput', false), testcomm2) == 1
        opt2commy = 2;
    elseif ismember(cellfun(@num2str, inducfinal(c, 4), 'UniformOutput', false), testcomm3) == 1
        opt2commy = 3;
    end
    opt2comm = vertcat(opt2comm, opt2commy); %#ok<AGROW>
end

%assigning the correct numerical question type in accordance with the labels
for d = 1:length(inducfinal)
    if strcmpi(inducfinal(d, 1), cellstr('Prim')) == 1
        nodetypey = 1;
    elseif strcmpi(inducfinal(d, 1), cellstr('Bound1')) == 1
        nodetypey = 2;
    elseif strcmpi(inducfinal(d, 1), cellstr('Bound2')) == 1
        nodetypey = 3;
    end
    nodetype = vertcat(nodetype, nodetypey); %#ok<AGROW>
end

%now putting together Cue Numerical Type, Option1 Numerical Type, Option 2 Numerical Type, and Numerical Question Type
orderingtogether = horzcat(cuecomm, opt1comm, opt2comm, nodetype); %cuecomm %opt1comm %opt2comm %nodetype

%now merging those numerical identifiers to the right of the induction questions
induct = horzcat(inducfinal, num2cell(orderingtogether));

%%%%%%%%%%%%%%%%%%
%%At this point:%%
%%%%%%%%%%%%%%%%%%

%induct: 
%bordertypeinwords (col1), 
%cueobj (col2), 
%opt1obj (col3), 
%opt2obj (col4), 
%obj2switchpos(col5),
%cuecomm# (col6), 
%opt2comm# (col7), 
%opt2comm# (col8), 
%bordertypenumerial# (col9)


%%%%%%%%%%%%%%%%%%
%NOW making sure that there's no more than 3q which are cued from same community
%%%%%%%%%%%%%%%%%%

rpchck = 1;

while rpchck > 0
    inductrand = randperm(length(induct));
    induct = induct(inductrand, :);
    bck2bck = 0;
    
    for trl = 3:length(induct)
        %first making sure that no cues of the same community follow up >3x
        %[col6]
        if isequal(induct{trl, 6}, induct{trl-1, 6}) == 1 && isequal(induct{trl-1, 6}, induct{trl-2, 6}) == 1
            bck2bck = bck2bck + 1;
            
        %now making sure that no more than two of the same border questions follow up >3x 
        %[col9]
        %elseif induct{trl, 9} && induct{trl-1, 9} == 4 || induct{trl-1, 9} && induct{trl-2, 9} == 4
        %    bck2bck = bck2bck + 1;
        end
        rpchck = bck2bck;
    end
end

%%%%%%%%%%%%%%%%
%Col1: QuestTypeInWords
%Col2: Cue
%Col3: Opt1 (Correct)
%Col4: Opt2 
%Col5: ChoiceSwitchOfOpt1
%Col6: Cue'sCommunity
%Col7: Opt1'sCommunity
%Col8: Opt2'sCommunity
%Col9: QuestType#
par.induct = induct;
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 3 TIMING %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Max time for each induction question
%Community structure onset time + ISI time
par.inductTimeout = 8; %max time for an item showing
par.inductTimefix = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------PART 4: Parser Task---------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

total_length_parse= 756; 

subtest_length_parse = 21; 
length_parse = total_length_parse/(subtest_length_parse*2); %mult by 2 here b/c this includes Ham paths too
curr_parse = Randi(21); %current position is a random number b/w 1:21
sequence_parse = NaN(0,length_parse); %creating an empty matrix of length 35
index_parse = 1; %1


% community 1 = 1,2,3,18,19,20,21; %ocean
% community 2 = 4,5,6,7,8,9,10; %desert
% community 3 = 11,12,13,14,15,16,17; %forest

boundarynodes_parse = [18,3,4,10,11,17];

graph_parse = [2,3,18,19,20,21; 3,18,19,20,21,1; 4,19,20,21,1,2; 5,6,7,8,9,3; 6,7,8,9,10,4; 7,8,9,10,4,5;...
        8,9,10,4,5,6; 9,10,4,5,6,7; 10,4,5,6,7,8; 5,6,7,8,9,11; 12,13,14,15,16,10; 13,14,15,16,17,11;...
        14,15,16,17,11,12; 15,16,17,11,12,13; 16,17,11,12,13,14; 17,11,12,13,14,15; 12,13,14,15,16,18;...
        19,20,21,1,2,17; 20,21,1,2,3,18; 21,1,2,3,18,19; 1,2,3,18,19,20];     
%the graph itself

%clockwise transitions
clockwise_parse = 1:21;

%counterclockwise transitions
counterclockwise_parse = 21:-1:1;

for j = 1:length_parse %18
    % generating a random walk for 21 items
    for i = 1:subtest_length_parse 
        next_parse = Randi(6); %random between 1 and 6
        curr_parse = graph_parse(curr_parse,next_parse); %changes curr index to going from one of 1/21 and then goes to the
        %column of the random of 1/4
        sequence_parse(1, index_parse) = curr_parse; %sequence starts at 1 and goes down that empty sequence matrix: 1st ROW
        sequence_parse(2, index_parse) = 1; %says that this is the normal sequence path via [1]; 2nd ROW
        if ismember(curr_parse, boundarynodes_parse) == 1
            sequence_parse(3, index_parse) = 1;
        else
            sequence_parse(3, index_parse) = 0;
        end
        index_parse = index_parse+1; %the index is that + 1
    end
    
    % choose a hamiltonian path for next 21 items
    direction_parse = Randi(2); %choose clockwise vs. counterclockwise, b/w 1 and 2 here
    if direction_parse == 1
        hamiltonian_parse=clockwise_parse; %the direction of the sequence is clockwise
    else
        hamiltonian_parse=counterclockwise_parse; %else the sequence goes counter
    end
    
    % fitting the two walks together
    currindex_parse = find(hamiltonian_parse == curr_parse);
    for i=1:subtest_length_parse
        currindex_parse = currindex_parse+1;
        if currindex_parse == 22; currindex_parse = 1; end
        sequence_parse(1, index_parse) = hamiltonian_parse(currindex_parse);
        if direction_parse == 1
            sequence_parse(2, index_parse) = 2;
        else
            sequence_parse(2, index_parse) = 3;
        end
        if ismember(sequence_parse(1, index_parse), boundarynodes_parse) == 1
            sequence_parse(3, index_parse) = 1;
        else
            sequence_parse(3, index_parse) = 0;
        end
        index_parse = index_parse+1;
    end
    curr_parse = hamiltonian_parse(currindex_parse);
end

%%%%%%%%%%%%%%%%
%Col1: Sequence
%Col2: Path Type: 1 = RandomWalk 2 = HamFor 3 = HamBck
%Col3: Boundary markers
par.sequence_parse = sequence_parse;
%%%%%%%%%%%%%%%%
%disp(sequence_parse) %uncomment later if you want to see

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% PART 3 TIMING %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Max time for each parse item
par.parseTimeout = 1.5; %max time for an item showing

%% SAVE IT ALL
header.parameters = par;
clearvars -except header;
save(sprintf('%s/%s_header', header.path.subjinfo, header.subjinfo))

%fin