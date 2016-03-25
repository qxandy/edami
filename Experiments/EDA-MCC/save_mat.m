fid=fopen('D:\mat2.txt','w');
[m,n]=size(a);
for i=1:m 
  for j=1:n 
     fprintf(fid,'%lf ',a(i,j));
  end
end
fprintf(fid,'\n');
[m,n]=size(b);
for i=1:m 
  for j=1:n 
     fprintf(fid,'%lf ',b(i,j));
  end
end
fclose(fid);