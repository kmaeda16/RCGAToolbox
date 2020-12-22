% This script runs simulateModel for all the folders
clear all;

diary('diary.txt');
cwd = pwd;
addpath(cwd);

%% Simulation
cd('Goldbeter1995_CircClock');
simulateModel('BIOMD0000000016_url.xml');
cd ..;

cd('Tyson1991_CellCycle6var');
simulateModel('BIOMD0000000005_url.xml');
cd ..;

cd('Chassagnole2002_CarbonMetabolism');
simulateModel('BIOMD0000000051_url.xml');
cd ..;

cd('Maeda2019_AmmoniumTransportAssimilation');
simulateModel('Maeda2019_RefinedActive_Kim_20190413_1.xml');
cd ..;

cd('Threestep');
simulateModel('threestep_SBML.xml');
cd ..;

cd('HIV');
simulateModel('hiv_SBML.xml');
cd ..;

%% 
rmpath(cwd);
diary off;
