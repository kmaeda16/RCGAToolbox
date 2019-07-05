diary('diary.txt');

cd('Goldbeter1995_CircClock');
simulateModel;
cd ..;

cd('Tyson1991_CellCycle6var');
simulateModel;
cd ..;

cd('Chassagnole2002_CarbonMetabolism');
simulateModel;
cd ..;

cd('Maeda2019_AmmoniumTransportAssimilation');
simulateModel;
cd ..;

diary off;
