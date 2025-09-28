function [Field] = CSR2Field(CSR, CSRtype)
% Input: CSR (number), eg 1 (which would correspond to a 50-50 CSR)
% CSRtype (string), eg "Fe++/Fe+". Not doing any checks if input makes sense


% In case user provides a char array, convert to string
CSRtype = string(CSRtype);


% process CSRtypoe string. At the ent of this, we will Element (tyoe_ion1 and 2 
% and chstate (chstat_ion1 and chstat_ion2).
% This will go wrong if CSRtype does not have a valid format.
ions = strsplit(CSRtype, "/");
chstat_ion1 = count(ions(1), "+");
chstat_ion2 = count(ions(2), "+");
if (abs(chstat_ion1 - chstat_ion2) ~= 1)
    error("Difference in charge must be +/- 1");
end
type_ion1 = extractBefore(ions(1), "+");
type_ion2 = extractBefore(ions(1), "+");
if (type_ion1 ~= type_ion2)
    error("This CSR doesn't make any sense");
end

%Not get Kingham curves for the element
KinghamCurve = kinghamcurve4Element(type_ion1);

% We only need the curves for the chstates that the user requested.
col_1 = KinghamCurve.(string(chstat_ion1) + "+");
col_2 = KinghamCurve.(string(chstat_ion2) + "+");

%the ratio of the columns...
ch_ratio = col_1./col_2;
ch_ratio(isinf(ch_ratio)) = NaN; % replace inf by nan
ch_ratio(ch_ratio == 0) = NaN; %also don;t allow for CSR of 0


% now just to a linear interpolation. Not very accurate but maybe good
% enough.
Fields = KinghamCurve.Field;
Fields(isnan(ch_ratio)) = [];
ch_ratio(isnan(ch_ratio)) = [];

if CSR>max(ch_ratio)
    Field = max(ch_ratio);
elseif CSR<min(ch_ratio)
    Field = min(ch_ratio);
end

[ch_ratio, order] = sort(ch_ratio, 'ascend');% interp1 want the input to be sorted
Fields = Fields(order);
    
Field = interp1(ch_ratio, Fields, CSR, 'linear');

end

