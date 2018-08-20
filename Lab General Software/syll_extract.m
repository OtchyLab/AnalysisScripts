%sound=A23_382775429_10_17_13_1_54.data;figure; specgram1(sound,512,44100,400,350);ylim([0 10000]);
syllA_index=[1.55 1.72];
syllB_index=[1.7 2.05];
syllC_index=[2.05 2.2 3.74 3.89 4.4 4.55 5.3 5.5 7.17 7.35 7.89 8.02];
syllD_index=[0.01 0.011];
syllE_index=[0.01 0.011];



syllA_index=syllA_index*22050;

for k = 1:(size(syllA_index,2)/2)
    a=syllA_index(2*k-1);
    b=syllA_index(2*k);
    syllA=cat(1,syllA,sound(a:b));
end;
figure; specgram1(syllA,512,22050,400,350);


syllB_index=syllB_index*22050;

for k = 1:(size(syllB_index,2)/2)
    a=syllB_index(2*k-1);
    b=syllB_index(2*k);
    syllB=cat(1,syllB,sound(a:b));
end;
figure; specgram1(syllB,512,22050,400,350);



syllC_index=syllC_index*22050;


for k = 1:(size(syllC_index,2)/2)
    a=syllC_index(2*k-1);
    b=syllC_index(2*k);
    syllC=cat(1,syllC,sound(a:b));
end;
figure; specgram1(syllC,512,22050,400,350);


syllD_index=syllD_index*22050;


for k = 1:(size(syllD_index,2)/2)
    a=syllD_index(2*k-1);
    b=syllD_index(2*k);
    syllD=cat(1,syllD,sound(a:b));
end;
figure; specgram1(syllD,512,22050,400,350);




syllE_index=syllE_index*22050;


for k = 1:(size(syllE_index,2)/2)
    a=syllE_index(2*k-1);
    b=syllE_index(2*k);
    syllE=cat(1,syllE,sound(a:b));
end;
figure; specgram1(syllE,512,22050,400,350);

