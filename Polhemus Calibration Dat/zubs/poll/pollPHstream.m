% pollPHstream.m
%
% usage: [newSampe]=pollPHstream(inds);
%  inds: columns of raw datastream to return
% 
% teh wrote it. 3.14.15

function [newSample]=pollPHstream(inds) %produces an error if no inds are passed
newSample(1,:)=[10*double(cell2mat(poldata('127.0.0.1',7234))) GetSecs]; %extract latest sample (1xN output)
newSample=newSample(inds);
