function rangeAdd_cllBck(fH)

% simple callback function to jump to next ion

eventName = get(fH,'SelectionType');

if strcmp(eventName,'alt')
   return; 
end
