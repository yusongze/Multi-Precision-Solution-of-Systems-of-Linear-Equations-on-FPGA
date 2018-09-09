while true
A = gallery('randsvd',1000,1,1);
% A = randn([1,9]);
temp = gallery('randsvd',1000,1,2);
b = temp(1:1000,1);
% b = randn([8,1]);
    if cond(A,2)<10
        if cond(b,2)<10
            break;
        end
    end
end
fileID = fopen('data1000by1000.h','w');
% fileID = fopen('/home/songze/data.h','w');
fprintf(fileID, "#include ""LU.h""\n");
fprintf(fileID, "float A_random[N][N]={");
fprintf(fileID, "%.60f,\n", A');
fprintf(fileID, "};\n");
fprintf(fileID, "float b_random[N]={");
fprintf(fileID, "%.60f,\n", b);
fprintf(fileID, "};\n");
cond(A,2)
