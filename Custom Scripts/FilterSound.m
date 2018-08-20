function [sound]=FilterSound(sound)


%Initial matrices
A=zeros(1,3);
B=zeros(1,3);

%Filtering loop
  for i=1:length(sound)
        A(1) = A(2);
        A(2) = A(3);
        A(3) = sound(i) /1.020353514; % GAIN value;
        B(1) = B(2);
        B(2) = B(3);
        B(3) =   (A(1) + A(3)) - 2 * A(2) + ( -0.9605029194 * B(1)) + (  1.9597070338 * B(2)); %Filter coeff
        sound(i) = B(3);
  end
