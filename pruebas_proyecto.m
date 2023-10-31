clear all; close all;
%Seleccionar Imagen
[Nom_I,PathName1] = uigetfile('C:\Users\Joan\Documents\Proyecto_Vision\*.*','Seleccione el archivo imagen');
%Leer Imagen
I=imread(Nom_I);
%información de la imagen
InfoI=imfinfo(Nom_I); 
W_I=InfoI.Width ;
H_I=InfoI.Height ;
Form=InfoI.Format;
%conversión de RGB a GrayScale.
I= rgb2gray(I);

%ploteo de la información
figure(1);
subplot(2,2,1);imshow(I);title(['Image, W=' num2str(W_I) 'H=' num2str(H_I)]);
subplot(2,2,2);imhist(I);title('General Histogram');
SH = sum(I);
SV = sum(I');
subplot(2,2,3);plot(0:W_I-1,SH);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
subplot(2,2,4);plot(0:H_I-1,SV);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;

%filtro paso bajo

Mask = 131; %131
order = 2;
diffV = zeros(1,length(SV));
aux = SV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 21; %131
order = 2;
aux = diffV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 131; %131
order = 2;
diffH = zeros(1,length(SH));
aux = SH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end

Mask = 21; %131
order = 2;
aux = diffH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end
    
%derivo
order = 2;
aux = diffV;
for i=1:order
    init = aux(1);
    for j=1:length(diffV)
        if j-1==0
            diffV(j) = 0;
        else
            diffV(j) = aux(j) - aux(j-1);
        end
    end
    aux = diffV;
end

order = 2;
aux = diffH;
for i=1:order
    init = aux(1);
    for j=1:length(diffH)
        if j-1==0
            diffH(j) = 0;
        else
            diffH(j) = aux(j) - aux(j-1);
        end
    end
    aux = diffH;
end

[pcks,locs] = findpeaks((-1)*diffV);

porcentaje = 0.1;

umbral = porcentaje*max(pcks);

j0 = 1;
j1 = length(pcks);
aux = pcks(1);

while (aux<umbral)
    j0 = j0 + 1;
    aux = pcks(j0);
end
aux = pcks(length(pcks));
while (aux<umbral)
    j1 = j1-1;
    aux = pcks(j1);
end

j0 = locs(j0);
j1 = locs(j1);

[pcks,locs] = findpeaks((-1)*diffH);

umbral = porcentaje*max(pcks);

i0 = 1;
i1 = length(pcks);
aux = pcks(1);

while (aux<umbral)
    i0 = i0 + 1;
    aux = pcks(i0);
end
aux = pcks(length(pcks));
while (aux<umbral)
    i1 = i1-1;
    aux = pcks(i1);
end

i0 = locs(i0);
i1 = locs(i1);

IR = I(j0:j1,i0:i1);

I(j0:j1,i0:i1) = 255;

[f,c] = size(IR);

hist = zeros(256,1);

for i=1:H_I
    for j=1:W_I
        hist(uint16(I(i,j))+1) = hist(uint16(I(i,j))+1) + 1;
    end
end

hist(256) = hist(256) - f*c;

for i=1:length(hist)
    if i>1
        hist(i) = hist(i) + hist(i-1);
    end
end
hist = hist/hist(256);

umbral = 0.2;
i=1;
while hist(i) <=umbral
    i = i + 1;
end

i = i - 1;

hist = zeros(256,1);

for k=1:f
    for j=1:c
        hist(uint16(IR(k,j))+1) = hist(uint16(IR(k,j))+1) + 1;
    end
end

for k=1:length(hist)
    if k>1
        hist(k) = hist(k) + hist(k-1);
    end
end
hist = hist/hist(256);

umbral = 0.005;
j=1;

while hist(j)<=umbral
    j = j + 1;
end
j = j - 2;

m = 255/(i-j);

T = zeros(256,1);
for k=1:256
    T(k) = uint8((k-1-j)*m);
end

for i=1:f
    for j=1:c
        IR(i,j) = T(IR(i,j)+1);
    end
end


figure(2);
subplot(2,2,1);imshow(IR);title(['Image, W=' num2str(c) 'H=' num2str(f)]);
subplot(2,2,2);imhist(IR);title('General Histogram');
SH = sum(IR);
SV = sum(IR');
subplot(2,2,3);plot(0:c-1,SH);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
subplot(2,2,4);plot(0:f-1,SV);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;


Mask = 11; %131
order = 2;
diffV = zeros(1,length(SV));
aux = SV;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffV)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffV(j) = diffV(j) + init;
            else
                if j+k>length(aux)
                    diffV(j) = diffV(j) + finit;
                else
                    diffV(j) = diffV(j) + aux(j+k);
                end
            end
        end
        diffV(j) = diffV(j)/Mask;
    end
    aux = diffV;
end

Mask = 11; %131
order = 2;
diffH = zeros(1,length(SH));
aux = SH;
for i=1:order
    init = aux(1);
    finit = aux(length(aux));
    for j=1:length(diffH)
        for k=(0:Mask-1)-(Mask-1)/2
            if j + k<=0
                diffH(j) = diffH(j) + init;
            else
                if j+k>length(aux)
                    diffH(j) = diffH(j) + finit;
                else
                    diffH(j) = diffH(j) + aux(j+k);
                end
            end
        end
        diffH(j) = diffH(j)/Mask;
    end
    aux = diffH;
end

umbral = 0.99*max([diffH(1), diffH(length(diffH))]);
HP = zeros(10,2);

[pcks,locs] = findpeaks((-1)*diffH);

k=0;

for i=1:length(locs)
    auxm = locs(i);
    auxM = locs(i);
    while diffH(auxm)<=umbral
        auxm = auxm - 1;
    end
    while diffH(auxM)<=umbral
        auxM = auxM + 1;
    end
    if auxm ~= auxM
        if k>=1
            if HP(k)~=auxm
                k = k + 1;
                HP(k,1) = auxm;
                HP(k,2) = auxM;
            end
        else
            k = k + 1;
            HP(k,1) = auxm;
            HP(k,2) = auxM;
        end
    end
end

umbral = 0.99*max([diffV(1), diffV(length(diffV))]);
VP = zeros(10,2);

[pcks,locs] = findpeaks((-1)*diffV);

k=0;

for i=1:length(locs)
    auxm = locs(i);
    auxM = locs(i);
    while diffV(auxm)<=umbral
        auxm = auxm - 1;
    end
    while diffV(auxM)<=umbral
        auxM = auxM + 1;
    end
    if auxm ~= auxM
        if k>=1
            if VP(k)~=auxm
                k = k + 1;
                VP(k,1) = auxm;
                VP(k,2) = auxM;
            end
        else
            k = k + 1;
            VP(k,1) = auxm;
            VP(k,2) = auxM;
        end
    end
end

figure(3);
plot(0:f-1,diffV);
