function ConvSizeChoose
    clear all;
    close all;
    %border = 10000;
    border = 512;

    c1 = DirectConvComplexity(360000, 2048)
    c2 = MaxFFTLenComplexity(360000, 2048)
    c3 = SegFFTLenComplexity(360000, 2048)
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
