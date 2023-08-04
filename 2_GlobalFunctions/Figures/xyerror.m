function xyerror(X,Y,p)

for ii = 1:size(X,1)
    plot([X(ii,p==50),X(ii,p==50)],[Y(ii,p==95),Y(ii,p==5)],...
        '-','color',[.75 .75 .75],'LineWidth',1.5)
    plot([X(ii,p==95),X(ii,p==5)],[Y(ii,p==50),Y(ii,p==50)],...
        '-','color',[.75 .75 .75],'LineWidth',1.5)
    plot([X(ii,p==50),X(ii,p==50)],[Y(ii,p==84),Y(ii,p==16)],...
        '-','color',[.5 .5 .5],'LineWidth',3)
    plot([X(ii,p==84),X(ii,p==16)],[Y(ii,p==50),Y(ii,p==50)],...
        '-','color',[.5 .5 .5],'LineWidth',3)
end

end