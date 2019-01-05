function DegreeRadius = eightradii(DegreeRadius, newcenters, edgexy, ns, newradii)
%% Radius 1
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
while edgexy(Y,X) == 0
    if edgexy(Y, X) == 1
        Point1 = [X, Y];
        DegreeRadius(1) = radiusdegree(Point1(1), Point1(2), X0, Y0);
    else
        Y = Y - 1;
        if Y < 1
            X = X0 + 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
        if edgexy(Y, X) == 1
            Point1 = [X, Y];
            DegreeRadius(1) = radiusdegree(Point1(1), Point1(2), X0, Y0);
        end

    end
end
 
%% Radius 2
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
ny = 1;
while edgexy(Y,X)== 0
    if edgexy(Y, X) == 1
        Point2 = [X, Y];
        DegreeRadius(2) = radiusdegree(Point2(1), Point2(2), X0, Y0);
    else
        Y = Y - 1;
        if Y < 1
            Y = Y0 - 1.*ny;
            X = X0;
            ny = ny + 1;
        end
        if edgexy(Y,X) == 1
            Point2 = [X, Y];
            DegreeRadius(2) = radiusdegree(Point2(1), Point2(2), X0, Y0);
        else
            X = X + 1;
        end
        if X > ns
            X = X0 + 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
    end
end
 
%% Radius 3
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
ny = 1;
while edgexy(Y,X) == 0
    if edgexy(Y, X) == 1
        Point3 = [X, Y];
        DegreeRadius(3) = radiusdegree(Point3(1), Point3(2), X0, Y0);
    else
        X = X + 1;
        if X > ns
            Y = Y0 + 1.*ny;
            X = X0;
            ny = ny + 1;
        end
        if edgexy(Y, X) == 1
            Point3 = [X, Y];
            DegreeRadius(3) = radiusdegree(Point3(1), Point3(2), X0, Y0);
        end
    end
end
%% Radius 4
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
ny = 1;
while edgexy(Y,X) == 0
    if edgexy(Y, X) == 1
        Point4 = [X, Y];
        DegreeRadius(4) = radiusdegree(Point4(1), Point4(2), X0, Y0);
    else
        Y = Y + 1;
        if Y > 208
            Y = Y0 + 1.*ny;
            X = X0;
            ny = ny + 1;
        end
        if edgexy(Y,X) == 1
            Point4 = [X, Y];
            DegreeRadius(4) = radiusdegree(Point4(1), Point4(2), X0, Y0);
        else
            X = X + 1;
        end
        if X > 208
            X = X0 + 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
    end
end
 
%% Radius 5
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
while edgexy(Y,X) == 0
    if edgexy(Y, X) == 1
        Point5 = [X, Y];
        DegreeRadius(5) = radiusdegree(Point5(1), Point5(2), X0, Y0);
    else
        Y = Y + 1;
        if Y > 208
            X = X0 + 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
        if edgexy(Y, X) == 1
            Point5 = [X, Y];
            DegreeRadius(5) = radiusdegree(Point5(1), Point5(2), X0, Y0);
        end
    end
end
%% Radius 6
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
ny = 1;
while edgexy(Y,X)== 0
    if edgexy(Y, X) == 1
        Point6 = [X, Y];
        DegreeRadius(6) = radiusdegree(Point6(1), Point6(2), X0, Y0);
    else
        Y = Y + 1;
        if Y > 208
            Y = Y0 + 1.*ny;
            X = X0;
            ny = ny + 1;
        end
        if edgexy(Y,X) == 1
            Point6 = [X, Y];
            DegreeRadius(6) = radiusdegree(Point6(1), Point6(2), X0, Y0);
        else
            X = X - 1;
        end
        if X < 1
            X = X0 - 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
        
    end
end
 
%% Radius 7
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
ny = 1;
while edgexy(Y,X) == 0
    if edgexy(Y, X) == 1
        Point7 = [X, Y];
        DegreeRadius(7) = radiusdegree(Point7(1), Point7(2), X0, Y0);
    else
        X = X - 1;
        if X < 1
            Y = Y0 - 1.*ny;
            X = X0;
            ny = ny + 1;
        if edgexy(Y, X) == 1
            Point7 = [X, Y];
            DegreeRadius(7) = radiusdegree(Point7(1), Point7(2), X0, Y0);
        end
        end
    end
end
%% Radius 8
X0 = round(newcenters(1) + (1/2).*newradii);
Y0 = round(newcenters(2) + (1/2).*newradii);
X = round(newcenters(1));
Y = round(newcenters(2));
nx = 1;
ny = 1;
while edgexy(Y,X)== 0
    if edgexy(Y, X) == 1
        Point8 = [X, Y];
        DegreeRadius(8) = radiusdegree(Point8(1), Point8(2), X0, Y0);
    else
        Y = Y - 1;
        if Y < 1
            Y = Y0 - 1.*ny;
            X = X0;
            ny = ny + 1;
        end
        if edgexy(Y,X) == 1
            Point8 = [X, Y];
            DegreeRadius(8) = radiusdegree(Point8(1), Point8(2), X0, Y0);
        else
            X = X - 1;
        end
        if X < 1
            X = X0 - 1.*nx;
            Y = Y0;
            nx = nx + 1;
        end
    end
end
end