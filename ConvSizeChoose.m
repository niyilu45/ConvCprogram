function ConvSizeChoose
    clear all;
    close all;
    rng(0);
    drawComplexityEn = false;

a = rand(10,1);
b = rand(3,1);
c = conv(a,b);
d = rand(25,1);
fid=fopen('SeqConvInA.txt','w');
fprintf(fid, '%f\n',a);
fclose(fid);
fid=fopen('SeqConvInB.txt','w');
fprintf(fid, '%f\n',b);
fclose(fid);
fid=fopen('SeqConvOutC.txt','w');
fprintf(fid, '%f\n',c);
fclose(fid);
fid=fopen('SeqConvInD.txt','w');
fprintf(fid, '%f\n',d);
fclose(fid);
% save('SeqConvInA.txt', 'a', '-ascii');
% save('SeqConvInB.txt', 'b', '-ascii');
% save('SeqConvOutC.txt', 'c', '-ascii');
error('This is for auto stop!');
% a = [0.1419; 0.4218; 0.9157; 0.7922; 0.9595; 0.6557; 0.0357; 0.8491; 0.9340; 0.6787; 0.7577];
% b = [0.7431; 0.3922; 0.6555; 0.1712];
% c = conv(a,b);
% A = fft(a, 16);
% B = fft(b, 16);
% Y = A.*B;
% y = ifft(Y);
% 
% s1 = a(1:5);
% s2 = a(6:10);
% s3 = [a(11);zeros(4,1)];
% S1 = fft(s1, 8);
% S2 = fft(s2, 8);
% S3 = fft(s3, 8);
% B = fft(b, 8);
% Y1 = S1 .* B;
% Y2 = S2 .* B;
% Y3 = S3 .* B;
% y1 = ifft(Y1);
% y2 = ifft(Y2);
% y3 = ifft(Y3);


    if drawComplexityEn == true
        %border = 10000;
        border = 512;

%         c1 = DirectConvComplexity(360000, 2048);
%         c2 = MaxFFTLenComplexity(360000, 2048);
%         c3 = SegFFTLenComplexity(360000, 2048);
        table = zeros(border, border);
        for n1 = 1:border
            for n2 = 1:n1
                c1 = DirectConvComplexity(n1, n2);
                c2 = MaxFFTLenComplexity(n1, n2);
                c3 = SegFFTLenComplexity(n1, n2);
                table(n1, n2) = MinTagOfThree(c1,c2,c3);
                %Copy to the another half of table
                table(n2, n1) = table(n1, n2);
            end
        end

        gca = pcolor(table);
        set(gca, 'LineStyle', 'none')
        colormap(gray(3))
    end
error('This is for auto stop!');
end

function c = DirectConvComplexity(n1, n2)
    c = n1 * n2;
end

function c = MaxFFTLenComplexity(n1, n2)
    fftOrder = ceil(log(n1+n2-1)/log(2));
    fftSize = 2^fftOrder;
    c = 4 * fftSize * (3*fftOrder + 1);
end

function c = SegFFTLenComplexity(n1, n2)
    if n1 < n2
        n  = n1;
        n1 = n2;
        n2 = n;
    end
    fftOrder = ceil(log(2*n2-1)/log(2));
    fftSize = 2^fftOrder;
    %s = ceil(n1 / fftSize); % This is error
    s = ceil(n1 / n2); % This is correct
    %s = ceil(n1 / (fftSize+1-n2)); % This is correct
    c = 4 * fftSize * ((2*s+1)*fftOrder + s);
end

function tag = MinTagOfThree(c1, c2, c3)
    tag = -1;
    if c1 <= c2 && c1 <= c3
        tag = 1;
    elseif c2 <= c1 && c2 <= c3
        tag = 2;
    elseif c3 <= c1 && c3 <= c2
        tag = 3;
    end
end
