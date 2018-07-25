function save_fig_hgexport(filename)

%choose options
if nargin < 1
    filename = 'matlab_plot';
end
figformat = 'pdf';
fntsze = 14;

target_path = 'D:\Research\Latex';
rootdir = pwd;

%apply settings
fig = gcf;
set(gca,'FontSize',fntsze)
myStyle = hgexport('factorystyle');
myStyle.Format = figformat;
myStyle.Resolution = 300;
myStyle.Units = 'inch';
myStyle.FixedFontSize = fntsze;
myStyle.FontSizeMin = fntsze;
myStyle.FontName = 'Helvetica';
myStyle.Bounds = 'tight';
myStyle.LineWidthMin = 1.5;
myStyle.LockAxesTicks = 'on';
set(fig,'PaperUnits','inches','PaperSize',[6,3.75]/2)
pos = get(fig,'PaperPosition');
set(fig,'PaperSize',[pos(3)*1.05,pos(4)*1.05]);

%output figure
cd(target_path)
hgexport(fig,[filename,'.',figformat],myStyle,'Format',figformat)
cd(rootdir)
