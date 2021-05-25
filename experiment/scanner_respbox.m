function scanner_respbox
%% Check for button box
% This is just a script...

intake=PsychHID('Devices');
extKeys=0;
for n=1:length(intake)
    if (intake(n).productID == 8 && strcmp(intake(n).usageName,'Keyboard')) % UT
    %if intake(n).productID==16385 && strcmp(intake(n).usageName,'Keyboard')
        extKeys=n;
    elseif strcmp(intake(n).usageName,'Keyboard')
        intKeys=n;
    end
end
if extKeys==0, disp('No Buttonbox Detected.'); end
dev = intKeys; % change this line to intKeys for running in the lab