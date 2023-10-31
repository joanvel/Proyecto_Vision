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

%defino un delta para la diferencia
delta = 1;

%Filtro paso bajo (promedio)
Mask0 = 1;
SHf = zeros(length(SH)+1-Mask0,1);

for i=1:length(SHf)
    for j=1:Mask0
        SHf(i,1) = SHf(i,1) + SH(i+j-1);
    end
    SHf(i,1) = SHf(i,1)/Mask0;
end
Mask1 = 500;

diffH = zeros(length(SHf)-delta,1);
diffHf = zeros(length(diffH)+1-Mask1,1);
diff2H = zeros(length(diffHf)-delta,1);
diffV = zeros(length(SV)-delta,1);
diff2V = zeros(length(diffV)-delta,1);

umbral0 = -0.5*10^4;
aux0 = 0;
aux1 = 0;
aux2 = 0;
aux3 = 0;
j0 = 0;
j1 = 1;
k0 = length(diff2H);
k1 = length(diff2V);
for i=1:(length(diffH))
    diffH(i,1) = SHf(i+delta)-SHf(i);
end
%filtro pasobajo a la primera derivada (promedio)

for i=1:length(diffHf)
    for j=1:Mask1
        diffHf(i,1) = diffHf(i,1) + diffH(i+j-1,1);
    end
    diffHf(i,1) = diffHf(i,1)/Mask1;
end


for i=1:(length(diff2H))
    diff2H(i,1) = diffHf(i+delta)-diffHf(i);
end
for i=1:(length(diffV))
    diffV(i,1) = SV(i+delta)-SV(i);
end
for i=1:(length(diff2V))
    diff2V(i,1) = diffV(i+delta)-diffV(i);
end
% while(aux0>umbral0)
%     j0 = j0 + 1;
%     aux0 = diff2H(j0);
% end
% while(aux1>umbral0)
%     aux1 = diff2H(k0);
%     k0 = k0 - 1;
% end
% while(aux2>umbral0)
%     j1 = j1 + 1;
%     aux2 = diff2V(j1);
% end
% while(aux3>umbral0)
%     aux3 = diff2V(k1);
%     k1 = k1 - 1;
% end
% j0 = j0 + delta;
% j1 = j1 + delta;
% k0 = k0 + delta;
% k1 = k1 + delta;
% 
% IR = I(j1:k1,j0:k0);
% IE = I;
% IE(j1:k1,j0:k0) = 255;

%ploteo de la información
% W_IR = k0-j0+1;
% H_IR = k1-j1+1;
% figure(1);
% subplot(2,2,1);imshow(IR);title(['Imagen recortada, W=' num2str(W_IR) 'H=' num2str(H_IR)]);
% subplot(2,2,2);imhist(IR);title('General Histogram');
% SHR = sum(IR);
% SVR = sum(IR');
% subplot(2,2,3);plot(0:W_IR-1,SHR);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
% subplot(2,2,4);plot(0:H_IR-1,SVR);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;
% 
% 
% figure(2);
% subplot(2,2,1);imshow(IE);title(['Image, W=' num2str(W_I) 'H=' num2str(H_I)]);
% subplot(2,2,2);imhist(IE);title('General Histogram');
% SHE = sum(IE);
% SVE = sum(IE');
% subplot(2,2,3);plot(0:W_I-1,SHE);xlabel('Horizontal Axe');ylabel('Cumulative sum of grays');title('Horizomtal profile');grid on;
% subplot(2,2,4);plot(0:H_I-1,SVE);xlabel('Vertical Axe');ylabel('Cumulative sum of grays');title('Vertical profile');grid on;
% 
% figure(3);
% imshow(IR);
figure(3);
plot(0:W_I-Mask0-delta,diffH);
figure(4);
plot(0:W_I-Mask0-Mask1+1-2*delta,diff2H);