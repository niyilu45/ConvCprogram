function ConvSizeChoose
    clear all;
    close all;
    rng(0);
    dataCheckEn = true;
    figSurfEn   = true;
    figMapEn    = true;
    border = 512;

    x1 = [1451;34562;25124;2652;1414;6346];
    x2 = [5151;125;61346;7247;124;62;1616;3146];
    y = conv(x1, x2);

% [convOut] = SegFFTConv(x1, x2);
% error('auto stop');
    if dataCheckEn
        c1            = DirectConvComplexity(100, 12);
        [c2, fftSize] = SegFFTConvComplexity(100, 12);
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
%         surf(AxisX, AxisY, cTable, Color);
        surf(AxisX, AxisY, cTable);
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

function [convOut] = SegFFTConv(x1, x2)
    % 1) Calc the output length and choose the fft size.
    [x1SizeM, x1SizeN] = size(x1);
    [x2SizeM, x2SizeN] = size(x2);
    if x2SizeM > x1SizeM
        tmp = x1;
        x1  = x2;
        x2  = tmp;
        [x1SizeM, x1SizeN] = size(x1);
        [x2SizeM, x2SizeN] = size(x2);
    end
    [~, fftSize] = SegFFTConvComplexity(x1SizeM, x2SizeM);
    convOutLen = x1SizeM + x2SizeM - 1;

    % 2) segment fft conv.
    idx = 1;
    block = fftSize + 1 - x2SizeM;
    segs = ceil(x1SizeM / block);
    convOutLenTmp = convOutLen + segs*block-x1SizeM;
    convOut = zeros(convOutLenTmp, 1);
    x1 = [x1; zeros(segs*block-x1SizeM, 1)];
    X2 = fft(x2, fftSize);
    for i = 1:segs
        % maybe the reason is that x1 need more zeros, then this part segment conv seq exceed the convOutlen
        % so we should guarantee that the conv out result correct firstly.
        sX1 = fft(x1(idx:idx+block-1, 1), fftSize);
        sY = sX1 .* X2;
        sy = ifft(sY, fftSize);
        convOut(idx:idx+fftSize-1, 1) = convOut(idx:idx+fftSize-1, 1) + sy;
        idx = idx + block;
    end
    convOut = convOut(1:convOutLen);
end

