
function [times]=timeSynch
global Sexp

t0=GetSecs-Sexp.tEXP;
tedf=Eyelink('TrackerTime');
times=[t0 tedf];