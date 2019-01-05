function edgexy = Crop(edgexy, a, b, c, d, ns)
%% Crop Images to remove noise
for x=1:a
    for y=1:ns
        edgexy(x,y)=0;
    end
end
for x=b:ns
    for y=1:ns
        edgexy(x,y)=0;
    end
end
for y=1:c
    for x=1:ns
        edgexy(x,y)=0;
    end
end
for y=d:ns
    for x=1:ns
        edgexy(x,y)=0;
    end
end
end
