function yyecrossplot( ...
    plothandle, maphandle, Yuse, Yeuse, proxytype, paleolat, paleolon, siteidx)


pt = unique(proxytype(siteidx));
for ii = 1:numel(pt)
    if pt(ii) == "d18a"
        marker = '^';
        color = hex2rgb('#94D2BD',1);
    elseif pt(ii) ==  "d18cforam"
        marker = 'o';
        color = hex2rgb('#0A9396',1);
    elseif pt(ii) == "d18cmacro"
        marker = 's';
        color = hex2rgb('#005F73',1);
    elseif pt(ii) == "d18p"
            marker = '^';
            color = hex2rgb('#0A9396',1);
    elseif pt(ii) ==  "mg"
        marker = 'o';
        color = hex2rgb('#B4BE65',1);
    elseif pt(ii) == "tex"
        marker = 's';
        color = hex2rgb('#FFB703',1);
    elseif pt(ii) == "uk"
        marker = 'v';
        color = hex2rgb('#9B2226',1);
    end
    pltidx = find(proxytype == pt(ii));
    yplot = repmat(Yuse(pltidx),1,size(Yeuse,2));
    yeplot = Yeuse(pltidx,:);
    s(ii) = scatter(plothandle, yplot(:), yeplot(:), 25, 'filled', ...
        'Marker', marker, 'MarkerFaceColor', color', 'MarkerEdgeColor', 'k');
    sm = scatterm(maphandle,paleolat(pltidx),paleolon(pltidx),20,'filled');
    sm.Children.Marker = marker;
    sm.Children.MarkerFaceColor = color;
    sm.Children.MarkerEdgeColor = 'none';

end 
oldlimits = axis(plothandle);
newlimits = equalaxes(oldlimits);
axis(plothandle, newlimits)
if any(contains(pt,'d18'))
    set(plothandle,'Ydir','Reverse','Xdir','Reverse')
end
xlabel(plothandle,'Y','FontName','Arial','FontSize',12,'FontWeight','bold')
ylabel(plothandle,'Ye','FontName','Arial','FontSize',12,'FontWeight','bold')
plot(plothandle,newlimits([1,2]),newlimits([3,4]),'k--')
if numel(pt)>0
    legpt = pt;
    legpt = strrep(legpt,'d18a','\delta^{18}O_{a}');
    legpt = strrep(legpt,'d18cmacro','\delta^{18}O_{cm}');
    legpt = strrep(legpt,'d18cforam','\delta^{18}O_{cf}');
    legpt = strrep(legpt,'d18p','\delta^{18}O_{p}');
    legpt = strrep(legpt,'mg','Mg/Ca');
    legpt = strrep(legpt,'tex','TEX_{86}');
    legpt = strrep(legpt,'uk',"U^{k'}_{37}");
    pause(0.5)
%     legend(plothandle,s,legpt,'Location','northwest','Orientation','vertical', ...
%         'FontName', 'Aril', 'FontSize', 8);
    pause(0.5)
end



end
