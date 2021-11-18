function ConvSizeChoose
    clear all;
    close all;
    rng(0);
    dataCheckEn = true;
    figSurfEn   = true;
    figMapEn    = true;
    border = 512;


    if dataCheckEn
        c1            = DirectConvComplexity(512, 100);
        [c2, fftSize] = SegFFTConvComplexity(512, 100);
    end

    cTable = zeros(border, border); % complexity ratio.
    rTable = zeros(border, border); % best choice result.
    sTable = zeros(border, border); % fftSize table.
    for n1 = 1:border
        for n2 = 1:border
            c1 = DirectConvComplexity(n1, n2);
            [c2, fftSize] = SegFFTConvComplexity(n1, n2);
            cTable(n1, n2) = c1 / c2;
            rTable(n1, n2) = c1 > c2;
            sTable(n1, n2) = fftSize;
            cTable(n2, n1) = cTable(n1, n2);
            rTable(n2, n1) = rTable(n1, n2);
            sTable(n2, n1) = sTable(n1, n2);
        end
    end

    if figSurfEn
        figure
        [AxisX, AxisY] = meshgrid(1:border, 1:border);
        Color = AxisX .* AxisY;
        surf(AxisX, AxisY, cTable, Color);
        colorbar
        shading interp
    end

    if figMapEn
        figure
        gca = pcolor(rTable);
        set(gca, 'LineStyle', 'none');
        colormap(gray(2));
    end

    error('This is for auto stop!');

end












function c = DirectConvComplexity(n1, n2)
    c = n1 * n2;
end

function [c, FFTSize] = SegFFTConvComplexity(n1, n2)
    if n1 < n2
        n  = n1;
        n1 = n2;
        n2 = n;
    end
    minFFTOrder = ceil(log(1+n2-1)/log(2)); % min n1 + n2 - 1
    maxFFTOrder = ceil(log(n1+n2-1)/log(2)); % max n1 + n2 - 1
    fftOrder = (minFFTOrder:maxFFTOrder).';
    fftSize = 2.^fftOrder;
    block = fftSize + 1 - n2;
    segs = ceil(n1 ./ block);
    complexity = 3 * ((2*segs+1) .* fftSize .* fftOrder / 2 + segs .* fftSize); %multiPerComplexNum*times*nlogn/2+multiNum
    [c, idx] = min(complexity);
    FFTSize = fftSize(idx);
end

