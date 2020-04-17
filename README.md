# Covid19_EERAModel
Simple COVID-19 simulation model with ABC-smc inference


version 0.2.2: inference COVID model over:
    * number of new daily hospitalised cases (NHS data) in Scotland/health board.
    * compute number of death as a sum of hospital deaths and community deaths
    * output only number of deaths from Hospitalised cases.

Bug reports to be solved:
    * issues in initiating simulations (can require multiple restarts)
    * multithreading not working on macOS
    *  Because of memory allocation limits for struct., the model is not able to output the total number of deaths (community + hospital).
    

    
Compilation:
g++ -Igsl -O3 -Wall -std=c++11 -std=gnu++11 -stdlib=libc++ IniFile.cpp main_version_0.2.2.cpp -lgsl -lgslcblas  -o covid-mod


