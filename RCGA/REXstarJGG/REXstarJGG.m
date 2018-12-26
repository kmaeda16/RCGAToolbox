function [ best, Population ] = REXstarJGG(Param)

Param = checkInputs(Param,'REXstarJGG');

[ best, Population ] = RCGA_Main(Param,@JGG);
