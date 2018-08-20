function ndx = time2ndx(t, sampRate);

ndx = (t.*sampRate) + 1;