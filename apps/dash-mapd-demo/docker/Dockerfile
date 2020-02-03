FROM omnisci/core-os-cpu:v4.5.0
ENV INSERT_SAMPLE_DATA=true
# startomnisci requires some inputs to install sample data, so we provide
# it using "printf".
CMD printf "Y\n2\n" | /omnisci/startomnisci

