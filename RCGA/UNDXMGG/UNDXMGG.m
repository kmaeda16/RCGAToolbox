function [ best, Population ] = UNDXMGG(Param)

Param = checkInputs(Param,'UNDXMGG');

[ best, Population ] = RCGA_Main(Param,@MGG);
