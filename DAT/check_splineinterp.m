load before.dat;
load between.dat;
load after.dat;

figure(1)
surf(before);

figure(2)
surf(between);

figure(3)
surf(after);

figure(4)
surf(after-between);
