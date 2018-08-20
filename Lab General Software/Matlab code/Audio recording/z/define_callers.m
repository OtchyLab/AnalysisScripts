% sets up the cell structure for callers
callers={};
callers.name={'Feedback','Playback'};
callers.prefix={{'f', 'F','S'},{'F','R'}};
callers.prefix_text={{'Motor','Feedback','Spontaneous'},{'Forward','Reverse'}};
callers.songtype={{'D','U'},{'B'}};
callers.songtype_text={{'Directed','Undirected'},{'BOS'}};

%cFeedback.name='Feedback';
%cFeedback.prefix={'f','F'};
%cFeedback.prefix_txt={'Motor','Feedback'};
%cFeedback.sontype={'D','U'};
%cFeedback.songtype_txt={'Directed','Undirected'};

%cPlayback.name='Playback';
%cPlayback.prefix={'F','R'};
%cPlayback.prefix_txt={{'Forward','Reverse'}};
%cPlayback.sontype='B';
%cPlayback.songtype_txt='BOS';
%callers=[cFeedback cPlayback];
