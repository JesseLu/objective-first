function [template] = field_template()
global DIMS_
dims = DIMS_;
% Setup constraints
tp = ones(dims);
tp([1,dims(1)],:) = 0;
tp(:,[1,dims(2)]) = 0;
template = [tp(:); tp(:); tp(:)];
