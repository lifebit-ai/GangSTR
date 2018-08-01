FROM ubuntu:xenial
MAINTAINER Pablo Prieto <pablo@lifebit.ai>

LABEL \
    description="GangSTR container"

RUN apt-get update -y && apt-get install -y \
    wget pkg-config zlib1g-dev

RUN wget -O /opt/GangSTR-1.4.tar.gz https://github.com/gymreklab/GangSTR/releases/download/v1.4/GangSTR-1.4.tar.gz && \
	cd /opt && tar -xzvf GangSTR-1.4.tar.gz && cd GangSTR-1.4 && ./install-gangstr.sh $PWD && mv /opt/GangSTR-1.4/bin/GangSTR /usr/local/bin/

ENTRYPOINT ["GangSTR"]
