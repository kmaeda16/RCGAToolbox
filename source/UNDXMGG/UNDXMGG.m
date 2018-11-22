function [ best, Population ] = UNDXMGG(Param)

[ best, Population ] = RCGA_Main(Param,@MGG);
