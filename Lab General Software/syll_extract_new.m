%sound=A24_382765148_10_16_12_21_24.data;figure; specgram1(sound,512,44100,400,350);ylim([0 10000]);
syllA_index=[1.4 2];




syllA_index=syllA_index*44100;

for k = 1:(size(syllA_index,2)/2)
    a=syllA_index(2*k-1);
    b=syllA_index(2*k);
    syllA=cat(1,syllA,sound(a:b));
end;
figure; specgram1(syllA,512,44100,400,350);


