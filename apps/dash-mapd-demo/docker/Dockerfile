FROM omnisci/core-os-cpu
ENV INSERT_SAMPLE_DATA=true
# startomnisci requires some inputs to install sample data, so we provide 
# it using "printf".
CMD printf "Y\n1\n" | /omnisci/startomnisci 

