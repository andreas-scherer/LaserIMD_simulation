clear all;
clc;
symmetrie='Ci';
for j=2:200
    [phi, theta, w]=sphgrid(symmetrie,j,'c');
    file=[symmetrie,sprintf('_%i.txt',j)];
    fid=fopen(file,'w');
    fprintf(fid,'%20s %20s %20s','phi','theta', 'weights');
    fprintf(fid,'\n%20f; %20f; %20f',[phi',theta', w']');
    fclose(fid);
end
