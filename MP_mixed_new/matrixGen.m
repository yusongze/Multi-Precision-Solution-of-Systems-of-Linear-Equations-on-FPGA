while true
A = gallery('randsvd',9,1,1);
A = randn([9,9]);
temp = gallery('randsvd',9,1,2);
b = temp(1:9,1);
% b = randn([8,1]);
    if cond(A,2)<10
        if cond(b,2)<10
            break;
        end
    end
end
fileID = fopen('data9by9.h','w');
% fileID = fopen('/home/songze/data.h','w');
fprintf(fileID, "#include ""LU.h""\n");
fprintf(fileID, "float A_random[N][N]={");
fprintf(fileID, "%.52f,\n", A');
fprintf(fileID, "};\n");
fprintf(fileID, "float b_random[N]={");
fprintf(fileID, "%.52f,\n", b);
fprintf(fileID, "};\n");
cond(A,2)
