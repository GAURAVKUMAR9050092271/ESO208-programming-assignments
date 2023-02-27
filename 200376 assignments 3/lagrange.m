clear;
fprintf("Give your no. of data and data points along with the the qwery abscissa value in a file named 'input1.txt'\n");
fprintf("Press ENTER to continue\n\n");
fid = fopen('input_lagrange.txt');
line =fgetl(fid);
n = sscanf(line, '%f ');
A = zeros(n,2);
for i=1:1:n
  line = fgetl(fid);
  B(i,1:1:2) = sscanf(line, '%f ');
end
line =fgetl(fid);
m = sscanf(line, '%f ');
Q = zeros(m,2);
for i=1:1:m
  line = fgetl(fid);
  Q(i,1) = sscanf(line, '%f ');
end

  for i=1:m
    Q(i,2) = polynomial(Q(i,1),B,n,m);
  end
  figure;
  for i=1:n
    plot(B(i,1),B(i,2),'*r',"markersize",20);
    hold on;
  end
  for i=B(1,1):0.01:B(n,1)
    plot(i,polynomial(i,B,n,m),'.y',"markersize",5);
    hold on;
  end
  hold off;
  Q;
  filename = 'output_lagrange.txt';
  fileID = fopen (filename, 'w');
  fprintf(fileID , "Lagrange Polynomial Method\n(X ,Y) is\n");
  fprintf(fileID , '%f\n',Q);
  fclose(fileID);

