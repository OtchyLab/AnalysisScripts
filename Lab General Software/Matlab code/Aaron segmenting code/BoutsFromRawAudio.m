function num = BoutsFromRawAudio(syllStartTimes, syllEndTimes)

boutInterval = .1; %sec
boutMinNumSyllables = 2; %sec
boutMinDuration = .3; %sec

%a bout must be a minimum x milliseconds and not to syllables can be y
%apart.
boutStartSyll = []; 
boutEndSyll = [];

if(length(syllStartTimes) < 2)
    num=0;
    return;
end

%determine bouts by syllable seperated by greater than boutInterval.
intervals = syllStartTimes(2:end) - syllEndTimes(1:end-1);
boutEndSyll = (find(intervals > boutInterval));
boutStartSyll = boutEndSyll + 1;
boutStartSyll = [1,boutStartSyll];
boutEndSyll = [boutEndSyll, length(syllStartTimes)];
    
%eliminate bouts that have a duration less than minDuration.
durations = syllEndTimes(boutEndSyll) - syllStartTimes(boutStartSyll);
realBout = find(durations > boutMinDuration);
boutStartSyll = boutStartSyll(realBout);
boutEndSyll = boutEndSyll(realBout);

%eliminate bouts that have fewer than minSyllables.
numSylls = boutEndSyll - boutStartSyll;
realBout = find(numSylls >= boutMinNumSyllables);
boutStartSyll = boutStartSyll(realBout);
boutEndSyll = boutEndSyll(realBout);

num=length(boutStartSyll);

    





