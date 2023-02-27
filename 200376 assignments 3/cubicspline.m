clear;
    method=input(' 1 -> Natural Cubic Spline\n 2 -> Not-a-Knot Cubic Spline\n 3 -> Periodic Cubic Spline\n 4 -> Clamped Cubic Spline\n','s');
    if method=='1'
        x = [];
        y = [];
        X = [];
        fileid=fopen('input.txt','r');
        data = textscan(fileid, '%*s %*s %*s %*s\n');
        data = textscan(fileid, '%f %f');
        x(1:numel(data{1}), 1) = data{1};
        y(1:numel(data{2}), 1) = data{2};
        data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s\n');
        data1 = textscan(fileid, '%f ');
        X(1:numel(data1{1}), 1) = data1{1};
        fclose(fileid);
        Naturalspline(x,y,X);
    end
    if method=='2'
        x = [];
        y = [];
        X = [];
        fileid=fopen('input.txt','r');
        data = textscan(fileid, '%*s %*s %*s %*s\n');
        data = textscan(fileid, '%f %f');
        x(1:numel(data{1}), 1) = data{1};
        y(1:numel(data{2}), 1) = data{2};
        data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s\n');
        data1 = textscan(fileid, '%f ');
        X(1:numel(data1{1}), 1) = data1{1};
        fclose(fileid);
        NotaKnotspline(x,y,X);
    end
    if method=='3'
        x = [];
        y = [];
        X = [];
        fileid=fopen('input.txt','r');
        data = textscan(fileid, '%*s %*s %*s %*s\n');
        data = textscan(fileid, '%f %f');
        x(1:numel(data{1}), 1) = data{1};
        y(1:numel(data{2}), 1) = data{2};
        data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s\n');
        data1 = textscan(fileid, '%f ');
        X(1:numel(data1{1}), 1) = data1{1};
        fclose(fileid);
        Periodicspline(x,y,X);
    end
    if method=='4'
        x = [];
        y = [];
        X = [];
        s=zeros(2,1);
        fileid=fopen('input.txt','r');
        data = textscan(fileid, '%*s %*s %*s %*s\n');
        data = textscan(fileid, '%f %f');
        x(1:numel(data{1}), 1) = data{1};
        y(1:numel(data{2}), 1) = data{2};
        data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s\n');
        data1 = textscan(fileid, '%f ');
        X(1:numel(data1{1}), 1) = data1{1};
        data = textscan(fileid, '%*s %*s %*s %*s %*s %*s %*s %*s %*s %*s\n');
        s(1,1)=fscanf(fileid,'%f',1);
        s(2,1)=fscanf(fileid,'%f',1);
        fclose(fileid);
        Clampedspline(x,y,X,s);
    end