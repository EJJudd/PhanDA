function makeshapefile(coords,filename)
for i = 1:length(coords)
    Shape(i).Geometry = 'Point'; 
    Shape(i).Lat=coords(i,1) ;  % latitude 
    Shape(i).Lon =coords(i,2);  % longitude 
    Shape(i).Name = i;          % some random attribute/name 
end
shapewrite(Shape, filename);
end