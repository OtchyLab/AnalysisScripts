if exist('ai');
    if isvalid(ai)
        if strcmp(get(ai, 'Running'), 'On')
            stop(ai);
        end
        delete(ai);
    end
end
if exist('ao')==1
    if isvalid(ao)
        if strcmp(get(ao, 'Running'), 'On')
            stop(ao);
        end
        delete(ao);
    end 
end
% Close the figure window.
if exist('hFig_ctr')
    if ishandle(hFig_ctr)
        delete(hFig_ctr);
    end
end
if exist('hFig')
    if ishandle(hFig)
        old_plot=do_plot; delete(hFig);   do_plot=old_plot;
    end
end
