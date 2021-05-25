function scanner_trigger
% Waits for the subject/scanner to press 5 or %
TRIGGERED = 0;
while ~TRIGGERED
    [keyIsDown, t, keyCode] = KbCheck(-1);
    if strcmp(KbName(keyCode), '5%')
        TRIGGERED = 1;
    end
    
end