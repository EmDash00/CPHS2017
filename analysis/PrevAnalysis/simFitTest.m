function [sysr sysh] = simFitTest()
% defaultArgs = {1 2};
% defaultArgs{1:nargin} = varargin{:};

sysr.o_num = 1;
sysr.o_den = 2;

% sys.num = rand(1, sys.o_num+1);
% sys.den = rand(1, sys.o_den);

sysr.num = randn(1, sysr.o_num+1);
sysr.den = [randn(1, sysr.o_den) 1];

sysr = stabilize(sysr);

sysr.z = roots(fliplr(sysr.num));
sysr.p = roots(fliplr(sysr.den));


% z = roots(sys.num);
% p = roots(sys.den);

F = primes(30).'*0.1;
A = 1./(2*pi*F);
P = rand(size(F))*2*pi;

t = 0:0.05:10;

R = sumOfSines(F,A,P,t);

%R = double(t>3);
Y = zsim(sysr, R, 0.01);
%Y = Y+0.5*randn(size(Y));


sysh = fitModel(R, Y, 1, 4);

Yh = zsim(sysh, R, 0);

figure(1); clf;
subplot(1,2,1)
plot(t,[R;Y], '.-'); hold on;
plot(t,Yh, '.-');

legend('input', 'output', 'model predicted');

cmap = lines;

subplot(1,2,2)

plot(real(sysh.z), imag(sysh.z), '.', 'Color', cmap(3, :), 'MarkerSize', 15); hold on;
plot(real(sysh.p), imag(sysh.p), '+', 'Color', cmap(3,:), 'MarkerSize', 5);
plot(real(sysr.z), imag(sysr.z), 'o', 'Color', cmap(2,:), 'MarkerSize', 5);
plot(real(sysr.p), imag(sysr.p), 'x', 'Color', cmap(2,:), 'MarkerSize', 5);

phi = [0:0.02:1]*2*pi;
plot(cos(phi), sin(phi), 'k:');

axis equal

end

function r = sumOfSines(f,a,p,t)
    r = sum((a*ones(size(t))) .* sin(2*pi*(f*t+p*ones(size(t)))));
end

function yh = zsim(sys, r, n)
   
    rh = [zeros(1, sys.o_den) r(:).'];
    yh = zeros(size(rh));
    
    for k = (sys.o_den+1):length(rh)
        yh(k) = sum(sys.num .* rh(k-sys.o_num:k)) - ...
            sum(sys.den(1:end-1) .* yh(k-sys.o_den:k-1)) + sqrt(n)*randn(1);
    end
    
    yh = yh(sys.o_den+1: end);
end

function sys2 = stabilize(sys)
    sys2 = sys;
    
    sys2.z = roots(fliplr(sys.num));
    sys2.z = sys2.z/(max(abs(sys2.z))+1);
    sys2.num = fliplr(poly(sys2.z));
    
    sys2.p = roots(fliplr(sys.den));
    sys2.p = sys2.p/(max(abs(sys2.p))+1);
    sys2.den = fliplr(poly(sys2.p));
    
end

function sysh = fitModel(inp, out, onum, oden)
    sidx = max([onum,oden])+1;
    
    inpMat = hankel(inp(sidx-onum:sidx), inp(sidx:end));
    outMat = -1*hankel(out(sidx-oden:sidx-1), out(sidx-1:end-1));
    
    %params = out(sidx:end)*pinv([inpMat; outMat]);
    params = out(sidx:end)/[inpMat; outMat];

    sysh.o_num = onum;
    sysh.o_den = oden;
    sysh.num = params(1:onum+1);
    sysh.den = [params(onum+2:end) 1];
    
    sysh.z = roots(fliplr(sysh.num));
    sysh.p = roots(fliplr(sysh.den));
end