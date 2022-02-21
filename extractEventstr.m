function eventName = extractEventstr(fullstr)
inds = union(strfind(fullstr, ':'),strfind(fullstr, ' '));
leagalInds = setdiff(1:length(fullstr), inds);
eventName=lower(fullstr(leagalInds));
eventName(eventName=='-') = '_';
eventName(eventName=='+') = 'p';

end