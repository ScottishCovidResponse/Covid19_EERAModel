FROM ubuntu:latest
RUN apt update -y
RUN apt upgrade -y
RUN apt install build-essential python3 python3-pip git vim -y
RUN python3 -m pip install cmake
RUN git clone https://github.com/ScottishCovidResponse/Covid19_EERAModel.git /home/Covid19_EERAModel
RUN git clone https://github.com/ScottishCovidResponse/data_pipeline_api.git /home/Covid19_EERAModel/api
RUN python3 -m pip install -r /home/Covid19_EERAModel/api/bindings/cpp/requirements.txt
WORKDIR /home/Covid19_EERAModel
ENV DEBIAN_FRONTEND=noninteractive
RUN apt install -y libgsl-dev pkg-config
RUN cmake -H. -Bbuild -DDATA_PIPELINE=/home/Covid19_EERAModel/api -DPython3_EXECUTABLE=$(which python3)
RUN cmake --build build
ENV EERAMODEL_ROOT=/home/Covid19_EERAModel
ENTRYPOINT [ "/usr/bin/bash" ]