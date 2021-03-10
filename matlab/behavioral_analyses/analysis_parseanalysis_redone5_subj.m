clc
clear all
clear all

%importing datafiles 
datanames

parsedata_all = []; 
relev_walk_data_all = []; 
final_parse_data_full = []; 

compiled_analysispath = '/Users/athula/Dropbox/Experiments/TesserScan/Behavioral Results';

cd(compiled_analysispath)

% output to text file
relev_outfname = sprintf('Parse_WalkType_Results_Overall'); 
relev_fid=fopen([relev_outfname '.txt'], 'w'); % open the file
fprintf(relev_fid, 'SUBJECT\tRand_Corr\tRand_Incorr\tHamFor_Corr\tHamFor_Incorr\tHamBack_Corr\tHamBack_Incorr\n'); %Create header

final_outfname = sprintf('Parse_OverWalks_Results_Overall'); 
final_fid=fopen([final_outfname '.txt'], 'w'); % open the file
fprintf(final_fid, 'SUBJECT\tRand_Corr\tRand_Incorr\tHam_Corr\tHam_Incorr\tBound_Prop\tOther_Prop\tParse_Diff\n'); %Create header

for s = 1:length(parsematload_all) %for each subject
    
%change directory
thisdatapath = datapaths_all{s}; 
indiv_analysispath = strcat(thisdatapath, '/analysis'); 
cd(thisdatapath);

parse_data = load(parsematload_all{s}); 
parse = parse_data.data.rundata; 

%Col1: Subj#
%Col2: Run#
%Col3: Trial#
%Col4: Obj#
%Col5: ObjSeq:  1 (random), 2 (hamfor), 3 (hamback)
%Col6: Boundary appearance (1 = boundary; 0 = non-boundary)
%Col7: Parse press (1 = pressed; 0 = no press)
%Col8: acc
%Col9: rt
%Col10: Parse Type:
    %parsetype = 'Bound_Ran';
        %parsetypenum = 11;

    %parsetype = 'Non_Ran';
        %parsetypenum = 21;

    %parsetype = 'Bound_ForH';
        %parsetypenum = 12;

    %parsetype = 'Non_ForH';
        %parsetypenum = 22;

    %parsetype = 'Bound_BckH';
        %parsetypenum = 13;

    %parsetype = 'Non_BckH';
        %parsetypenum = 23;

%the community structure used in the expriment
comm1 = [1,2,3,18,19,20,21]; %OCEAN
comm1pri = [1,2,19,20,21];
comm1bou = [18,3];

comm2 = [4,5,6,7,8,9,10]; %DESERT
comm2pri = [5,6,7,8,9];
comm2bou = [4,10];

comm3 = [11,12,13,14,15,16,17]; %FOREST
comm3pri = [12,13,14,15,16];
comm3bou = [11,17];

%relevant boundary nodes of Hamiltonian walks when entering into new one
hamforcomm_bound = [18, 4, 11];
hambckcomm_bound = [3, 17, 10];


%Now begin the compiling + analysis of walk types + parses 
randtotal = 0; %total number of items in random walk
randboundparsetotal = 0; %total numberof parses made at boundaries
rand_prim_parse_total = 0;
randincorr = 0;
randacc = 0;

hamfortotal = 0; %total number of items in hamilton forward walk
hamforboundparsetotal = 0; %total numberof parses made at boundaries
hamfor_prim_parse_total = 0;
hamforacc = 0;

hambcktotal = 0; %total number of items in hamilton backward walk
hambckboundparsetotal = 0; %total number of parses made at boundaries
hambck_prim_parse_total = 0;
hambckacc = 0;

randnon = 0;
randbou = 0;

hamfnon = 0;
hamfbou = 0;

hambnon = 0;
hambbou = 0;

%^%^%total number of relevant boundary items in random walks
rel_randbou = 0;

%^%^%total number of relevant non-boundary items in random walks
rel_randpri= 0;

%^%^%total number of relevant boundary items in ham forward walks
rel_hamforbou = 0;

%^%^%total number of relevant non-boundary items in ham forward walks
rel_hamforpri = 0;

%^%^%total number of relevant boundary items in ham backward walks
rel_hambckbou = 0;

%^%^%total number of relevant non-boundary items in ham backward walks
rel_hambckpri = 0;

for t = 1:length(parse)
    
    %getting variable names from parse to make it easier for interpretation
%Col1: Subj#
%Col2: Run#
%Col3: Trial#
%Col4: Obj#
%Col5: ObjSeq:  1 (random), 2 (hamfor), 3 (hamback)
%Col6: Boundary appearance (1 = boundary; 0 = non-boundary)
%Col7: Parse press (1 = pressed; 0 = no press)
%Col8: acc
%Col9: rt
%Col10: Parse Type:
    %parsetype = 'Bound_Ran';
        %parsetypenum = 11;

    %parsetype = 'Non_Ran';
        %parsetypenum = 21;

    %parsetype = 'Bound_ForH';
        %parsetypenum = 12;

    %parsetype = 'Non_ForH';
        %parsetypenum = 22;

    %parsetype = 'Bound_BckH';
        %parsetypenum = 13;

    %parsetype = 'Non_BckH';
        %parsetypenum = 23;
    
    subj = parse(t, 1); 
    run_num = parse(t, 2); 
    trial_num = parse(t, 3); 
    obj_num = parse(t, 4);
    walk_type = parse(t, 5);
    acc = parse(t, 8); 
    rt = parse(t, 9); 
    parse_type = parse(t, 10); 
    
    %----% getting together the last 4 trials to calculate inclusion within a community
    if t > 4
        lastfour = [parse(t-1, 4), parse(t-2, 4), parse(t-3, 4), parse(t-4, 4)];
    else
        lastfour = 0;
    end
    
    %----% calculating the overall total number of rand/hamfor/hambck trials
    if walk_type == 1
        randtotal = randtotal + 1;
    elseif walk_type == 2;
        hamfortotal = hamfortotal + 1;
    elseif walk_type == 3;
        hambcktotal = hambcktotal + 1;
    end
    
    %--------------------------------------------------------------------------------------------------------
    %----% calculating the OVERALL INCORRECT(b/c acc = 0) parses made at rand/hamfor/hambck trials
    %randwalks: aka incorrect parsing at the primary nodes
    %hamfor/hamback: aka incorrect parsing at primary and opposite bound nodes
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 1 && acc == 0 %incorrect primary parses made at all random trails at boundary items
        rand_prim_parse_total = rand_prim_parse_total + 1;
    elseif walk_type == 2 && acc == 0; %incorrect primary parses made at all hamfor trials at boundary items
        hamfor_prim_parse_total = hamfor_prim_parse_total + 1;
    elseif walk_type == 3 && acc == 0; %incorrect primary parses made at all hambck trials at boundary items
        hambck_prim_parse_total = hambck_prim_parse_total + 1;
    end
    
    %--------------------------------------------------------------------------------------------------------    
    %----% calculating the LEGITIMATE INCORRECT (b/c acc = 0) parses made
    %within Rand Walks, at boundary nodes of a community after four trials in
        %aka not parsing at the boundary nodes four trials in
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 1 %for random walks
        if ismember(lastfour(1), comm1) == 1 && ismember(lastfour(2), comm1) == 1 && ismember(lastfour(3), comm1) == 1 && ismember(lastfour(4), comm1) == 1 %if last four items member of a community
            if ismember(obj_num, comm1bou) == 1 && acc == 0 %if the current object is a boundary item of the same community, after the previous 4 members are also from the same community 
                randincorr = randincorr + 1; %add to the rand inacc count
            end
        elseif ismember(lastfour(1), comm2) == 1 && ismember(lastfour(2), comm2) == 1 && ismember(lastfour(3), comm2) == 1 && ismember(lastfour(4), comm2) == 1 %if last four items member of a community
            if ismember(obj_num, comm2bou) == 1 && acc == 0 %if the current object is a boundary item of the same community, after the previous 4 members are also from the same community 
                randincorr = randincorr + 1; %add to the rand inacc count
            end
        elseif ismember(lastfour(1), comm3) == 1 && ismember(lastfour(2), comm3) == 1 && ismember(lastfour(3), comm3) == 1 && ismember(lastfour(4), comm3) == 1 %if last four items member of a community
            if ismember(obj_num, comm3bou) == 1 && acc == 0 %if the current object is a boundary item of the same community, after the previous 4 members are also from the same community 
                randincorr = randincorr + 1; %add to the rand inacc count
            end
        end
    end

    %--------------------------------------------------------------------------------------------------------
    %----%This is calculating the ACCURACY with respect parsing in Rand
    %for random walks: aka parsing at boundary nodes 4 walks in
        %checking to see if the last 4 items in the series was within
        %one of the three communities (since only then is it counted)
        %if that's true, then add to the random parse accuracy
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 1 && acc == 1 %random walk when parsed at boundaries
        randboundparsetotal = randboundparsetotal+1; %ANY boundary parses made in Rand Walks
        if ismember(lastfour(1), comm1) == 1 && ismember(lastfour(2), comm1) == 1 && ismember(lastfour(3), comm1) == 1 && ismember(lastfour(4), comm1) == 1 %if last four items member of a community
            randacc = randacc + 1; %CORRECT boundary parses made in Rand Walks, b/c this is the true tally (at least 4 items into the same community) 
        elseif ismember(lastfour(1), comm2) == 1 && ismember(lastfour(2), comm2) == 1 && ismember(lastfour(3), comm2) == 1 && ismember(lastfour(4), comm2) == 1 %if last four items member of a community
            randacc = randacc + 1; %CORRECT boundary parses made in Rand Walks, b/c this is the true tally (at least 4 items into the same community) 
        elseif ismember(lastfour(1), comm3) == 1 && ismember(lastfour(2), comm3) == 1 && ismember(lastfour(3), comm3) == 1 && ismember(lastfour(4), comm3) == 1 %if last four items member of a community
            randacc = randacc + 1; %CORRECT boundary parses made in Rand Walks, b/c this is the true tally (at least 4 items into the same community) 
        end
    end
    
    %--------------------------------------------------------------------------------------------------------
    %----%This is calculating the ACCURACY with respect parsing in HamFor & HamBck
    %for hamfor and hamback walks: aka parsing at the specified boundary nodes
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 2 && acc == 1 %hamforward when parsed at boundaries
        hamforboundparsetotal = hamforboundparsetotal + 1; %ANY boundary parses made at hamfor trails
        if ismember(obj_num, hamforcomm_bound) == 1 %now checking if the obj# is a part of the select hamfor correct memeber
            hamforacc = hamforacc + 1; %CORRRECT boundary parses made in HamFor Walk
            %if that's true, then add to hamfor parse accuracy
        end
    elseif walk_type == 3 && acc == 1 %hambackward when parsed at boundaries
        hambckboundparsetotal = hambckboundparsetotal + 1; %ANY boundary parses made at hambck trails
        if ismember(obj_num, hambckcomm_bound) == 1 %now checking if the obj# is a part of the select hambck correct memeber
            hambckacc = hambckacc + 1; %CORRRECT boundary parses made in HamBck Walks 
            %if that's true, then add to hambck parse accuracy
        end
    end
    
    %--------------------------------------------------------------------------------------------------------    
    %----%This is calculating the TOTAL RELEVANT number of bound parses
    %(without respect to where participants parse etc. just in the raw data)
    %for random walks, for hamfor walks, for hambck walks
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 1 %for random walks
        if ismember(lastfour(1), comm1) == 1 && ismember(lastfour(2), comm1) == 1 && ismember(lastfour(3), comm1) == 1 && ismember(lastfour(4), comm1) == 1 %if last four items member of a community
            if ismember(obj_num, comm1bou) == 1 %if the current object is a boundary item of the same community
                rel_randbou = rel_randbou + 1; %add to the relevant boundary count
            elseif ismember(obj_num, comm1pri) == 1
                rel_randpri = rel_randpri + 1; %if current object not a boundary item memmber, add to primary item count
            end
        elseif ismember(lastfour(1), comm2) == 1 && ismember(lastfour(2), comm2) == 1 && ismember(lastfour(3), comm2) == 1 && ismember(lastfour(4), comm2) == 1 %if last four items member of a community
            if ismember(obj_num, comm2bou) == 1 %if the current object is a boundary item of the same community
                rel_randbou = rel_randbou + 1; %add to the relevant boundary count
            elseif ismember(obj_num, comm2pri) == 1
                rel_randpri = rel_randpri + 1; %if current object not a boundary item memmber, add to primary item count
            end
        elseif ismember(lastfour(1), comm3) == 1 && ismember(lastfour(2), comm3) == 1 && ismember(lastfour(3), comm3) == 1 && ismember(lastfour(4), comm3) == 1 %if last four items member of a community
            if ismember(obj_num, comm3bou) == 1 %if the current object is a boundary item of the same community
                rel_randbou = rel_randbou + 1; %add to the relevant boundary count
            elseif ismember(obj_num, comm3pri) == 1
                rel_randpri = rel_randpri + 1; %if current object not a boundary item memmber, add to primary item count
            end
        end
            %checking to see if the last four members are a part of comm1,2,3
                %if they are, check if the obj# is boundary of comm1,2,3
                    %if so, then add to the relative_boundary count
                    %else, then add to the relative_primary count
    elseif walk_type == 2 %for hamfor walks
        if ismember(obj_num, hamforcomm_bound) == 1
            rel_hamforbou = rel_hamforbou + 1;
        else
            rel_hamforpri = rel_hamforpri + 1;
        end
        %checking to see if the member is part of a hamfor walk
            %if they are, check if the obj# is a part of the hamfor parse set
            %if so, then add to the rel_hamforbou count
            %else, then add to the rel_hamforpri count
                %(which would contain both primary and the incorrect bound)
    elseif walk_type == 3 %for hambck walks
        if ismember(obj_num, hambckcomm_bound) == 1
            rel_hambckbou = rel_hambckbou + 1;
        else
            rel_hambckpri = rel_hambckpri + 1;
        end
        %checking to see if the member is part of a hambck walk
            %if they are, check if the obj# is a part of the hambck parse set
            %if so, then add to the rel_hambckbou count
            %else, then add to the rel_hambckpri count
                %(which would contain both primary and the incorrect bound)
    end
    
    %--------------------------------------------------------------------------------------------------------
    %----%This is calculating the TOTAL NUMBER of primary or boundary items in: Rand, HamF, and HamB
    %--------------------------------------------------------------------------------------------------------
    if walk_type == 1 %random walk
        if ismember(obj_num, comm1pri) == 1
            randnon = randnon + 1;
        elseif ismember(obj_num, comm1bou) == 1
            randbou = randbou + 1;
        elseif ismember(obj_num, comm2pri) == 1
            randnon = randnon + 1;
        elseif ismember(obj_num, comm2bou) == 1
            randbou = randbou + 1;
        elseif ismember(obj_num, comm3pri) == 1
            randnon = randnon + 1;
        elseif ismember(obj_num, comm3bou) == 1
            randbou = randbou + 1;
        end
        %if it's a randomwalk
        %then check if it's a part of the primary of the 3 comm
        %or if it's part of boundary of the 3 comm
    elseif walk_type == 2 %hamfor walk
        if ismember(obj_num, comm1pri) == 1
            hamfnon = hamfnon + 1;
        elseif ismember(obj_num, comm1bou) == 1
            hamfbou = hamfbou + 1;
        elseif ismember(obj_num, comm2pri) == 1
            hamfnon = hamfnon + 1;
        elseif ismember(obj_num, comm2bou) == 1
            hamfbou = hamfbou + 1;
        elseif ismember(obj_num, comm3pri) == 1
            hamfnon = hamfnon + 1;
        elseif ismember(obj_num, comm3bou) == 1
            hamfbou = hamfbou + 1;
        end
        %if it's a hamfor
        %then check if it's a part of the primary of the 3 comm
        %or if it's part of boundary of the 3 comm
    elseif walk_type == 3 %hambck walk
        if ismember(obj_num, comm1pri) == 1
            hambnon = hambnon + 1;
        elseif ismember(obj_num, comm1bou) == 1
            hambbou = hambbou + 1;
        elseif ismember(obj_num, comm2pri) == 1
            hambnon = hambnon + 1;
        elseif ismember(obj_num, comm2bou) == 1
            hambbou = hambbou + 1;
        elseif ismember(obj_num, comm3pri) == 1
            hambnon = hambnon + 1;
        elseif ismember(obj_num, comm3bou) == 1
            hambbou = hambbou + 1;
        end
        %if it's a hambck
        %then check if it's a part of the primary of the 3 comm
        %or if it's part of boundary of the 3 comm
    end
    
end

%--------------------------------------------------------------------------------------------------------
% consilidating relevant terms
%--------------------------------------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  RANDOM WALK TRIALS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% total number of random walks
randtotal;

%total number of boundary items in random walks
randbou;

%total number of non-boundary items in random walks
randnon;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant boundary items in random walks 
%(i.e. total real possibilities of boundary items parsed)
rel_randbou;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant non-boundary items in random walks
%(i.e. total real possibilities of non-boundary items parsed)
rel_randpri;

% parses made at ANY boundary nodes
randboundparsetotal;

% parses made at primary nodes
rand_prim_parse_total;

% total number of parses made
totalparserand = randboundparsetotal + rand_prim_parse_total;

% parses made at ONLY CORRECT boundary nodes 
%[i.e. made at boundary nodes, 4 items into a community]
randacc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% FORWARD HAMILTON WALK TRIALS %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total number of forward ham walks
hamfortotal;

%total number of boundary items in ham forward walks
hamfbou;

%total number of non-boundary items in ham forward walks
hamfnon;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant boundary items in ham forward walks
%(i.e. total real possibilities of boundary items parsed)
rel_hamforbou;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant non-boundary (including primary other boundary items) and  items in ham forward walks
%(i.e. total real possibilities of non-relevant boundary items parsed)
rel_hamforpri;

% parses made at ANY boundary nodes
hamforboundparsetotal;

% total incorrect parses made at primary nodes
hamfor_prim_parse_total;

% total number of parses made at ANY boundary nodes or primary nodes
totalparsehamfor = hamforboundparsetotal + hamfor_prim_parse_total;

% accurate parses [made boundary 9, 14 and 4 nodes]
% parses made at ONLY CORRECT boundary nodes 
hamforacc;

%total hambck incorrect
hamforincorrtotal = hamfor_prim_parse_total + (hamforboundparsetotal - hamforacc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% BACKWARD HAMILTON WALK TRIALS %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% total number of backward ham walks
hambcktotal;

%total number of boundary items in ham backward walks
hambbou;

%total number of non-boundary items in ham backward walks
hambnon;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant boundary items in ham backward walks
%(i.e. total real possibilities of boundary items parsed)
rel_hambckbou;

%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%^%%^%
%^%^%total number of relevant non-boundary items (including primary other boundary items) in ham backward walks
%(i.e. total real possibilities of non-relevant boundary items parsed)
rel_hambckpri;

% parses made at ANY boundary nodes
hambckboundparsetotal;

% total incorrect parses made at primary nodes
hambck_prim_parse_total;

% total number of parses made at ANY boundary nodes or primary nodes
totalparsehambck = hambckboundparsetotal + hambck_prim_parse_total;

%total hambck incorrect
hambckincorrtotal = hambck_prim_parse_total + (hambckboundparsetotal - hambckacc);

% accurate parses [made boundary 3, 13 and 8 nodes]
hambckacc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LISTING RELEVANT INFORMATION TOGETHER INTO ONE MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parsedata = [subj, randtotal, randbou, randnon, totalparserand, randboundparsetotal, rand_prim_parse_total, rel_randbou, rel_randpri, randacc, randincorr, hamfortotal, hamfbou, hamfnon, totalparsehamfor, hamforboundparsetotal, hamfor_prim_parse_total, rel_hamforbou, rel_hamforpri, hamforacc, hamforincorrtotal, hambcktotal, hambbou, hambnon, totalparsehambck, hambckboundparsetotal, hambck_prim_parse_total, rel_hambckbou, rel_hambckpri, hambckacc, hambckincorrtotal];
parsedata_all = [parsedata_all; parsedata]; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Utimately Relevant Comparsions %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RANDOM WALKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%rand correct
%rand correct responses at relative bound spots/relative boundaries
randcorrprop = randacc/rel_randbou;

%rand incorrect
%rand incorrect responses at relative prim spots/relative primary 
randincorrprop = rand_prim_parse_total/rel_randpri;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HAMILTONIAN FORWARD WALKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hamfor correct
hamforcorrprop = hamforacc/rel_hamforbou;

%hamfor incorrect
hamforincorrprop = hamforincorrtotal/rel_hamforpri;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HAMILTONIAN BACKWARD WALKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%hamback correct
hambackcorrprop = hambckacc/rel_hambckbou;

%hamback incorrect
hambackincorrprop = hambckincorrtotal/rel_hambckpri;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OVERALL HAMILTONIAN WALKS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hamcorr = [hamforcorrprop, hambackcorrprop];
hamcorrprop = mean(hamcorr);
hamincorr = [hamforincorrprop, hambackincorrprop];
hamincorrprop = mean(hamincorr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OVERALL CORRECT vs. INCORRECT ACROSS WALKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
commcorr = [randcorrprop, hamcorrprop];
commprop = mean(commcorr);
otherincorr = [randincorrprop, hamincorrprop];
otherprop =  mean(otherincorr);
parse_diff = commprop - otherprop; 

%%%%
relev_walk_data = [subj, randcorrprop, randincorrprop, hamforcorrprop, hamforincorrprop, hambackcorrprop, hambackincorrprop];
relev_walk_data_all  = [relev_walk_data_all; relev_walk_data]; 

final_parse_data = [subj, randcorrprop, randincorrprop, hamcorrprop, hamincorrprop, commprop, otherprop, parse_diff];
final_parse_data_full = [final_parse_data_full; final_parse_data]; 

%change directory to analysis directory
%save out overall proportions into .txt files
cd(compiled_analysispath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relevant walk data output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output to .txt file
fprintf(relev_fid, '%d\t %d\t %d\t %d\t %d\t %d\t %d\n',...
    subj, randcorrprop, randincorrprop, hamforcorrprop, hamforincorrprop, hambackcorrprop, hambackincorrprop); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall parse data output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output to .txt file
fprintf(final_fid, '%d\t %d\t %d\t %d\t %d\t %d\t %d\t %d\n',...
    subj, randcorrprop, randincorrprop, hamcorrprop, hamincorrprop, commprop, otherprop, parse_diff); 


end
 
close all

%change directory to analysis directory
%save out overall proportions into .mat files

cd(compiled_analysispath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Relevant walk data output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output to .mat file
save(sprintf('Parse_WalkType_Results_Overall'),'relev_walk_data_all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Overall parse data output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save output to .mat file
save(sprintf('Parse_OverWalks_Results_Overall'),'final_parse_data_full');
