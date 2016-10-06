%   PLOTTRISTREAM Plot streamlines computed by TRISTREAM
%
%   PlotTriStream(FlowP) plots streamlines FlowP produced by TRISTREAM
%
%   PlotTriStream(FlowP,style) plots streamlines using the line style 
%   specified in style. The default is style='r-'.
%
%   h=PlotTriStream(FlowP,style) additionally returns a vector h of handles
%   to each streamline which can be used to change its display properties.
function h=PlotTriStream(FlowP,style)

if(nargin<2)
    style='r-';
end

Npath=length(FlowP);

h=zeros(1,Npath);

hold on;
for p=1:Npath    
    h(p)=plot(FlowP(p).x,FlowP(p).y,style);
end
hold off

