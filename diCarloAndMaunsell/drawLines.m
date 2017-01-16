function h = drawLines(x,y, varargin)

color = 'k';
linestyle = '--';
assignopts(who, varargin);

if exist('x','var');
    if ~isempty(x)
        
        if length(x) > 1
            h = line([x; x], [repmat(get(gca,'ylim'),[length(x) 1])]','linestyle',linestyle);
            
        else
            h = line([x x],get(gca,'ylim'),'color',color,'linestyle',linestyle);
        end
    end
    
end

if exist('y','var');
    if ~isempty(y)
        h = line(get(gca,'xlim'),[y y],'color',color,'linestyle',linestyle);
    end
end


