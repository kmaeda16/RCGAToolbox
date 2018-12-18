function [ best, Population ] = REXstarJGG(Param)

[ best, Population ] = RCGA_Main(Param,@JGG);
