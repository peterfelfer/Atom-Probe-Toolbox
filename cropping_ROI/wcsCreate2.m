function wcsh = wcsCreate2(size,loc,displayName,ax)
% version of WCSs that uses lines and a grouping object. To supersede older
% versions

if not(exist('ax','var'))
    ax = gca;
end

if ~exist('displayName','var')
    displayName = 'WCS';
end

WCSzaxis = [0,0,0 ; 0,0,size];
WCSzaxis = WCSzaxis + repmat(loc,2,1);
WCSyaxis = [0,0,0 ; 0,size,0];
WCSyaxis = WCSyaxis + repmat(loc,2,1);
WCSxaxis = [0,0,0 ; size,0,0];
WCSxaxis = WCSxaxis + repmat(loc,2,1);

lh(1) = line(WCSxaxis(:,1),WCSxaxis(:,2),WCSxaxis(:,3),'LineWidth',4,'Color','r');
lh(2) = line(WCSyaxis(:,1),WCSyaxis(:,2),WCSyaxis(:,3),'LineWidth',4,'Color','g');
lh(3) = line(WCSzaxis(:,1),WCSzaxis(:,2),WCSzaxis(:,3),'LineWidth',4,'Color','b');

wcsh = hgtransform('Parent',ax);
wcsh.DisplayName = displayName;
set(lh,'Parent',wcsh);
