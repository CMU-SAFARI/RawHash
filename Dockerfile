FROM gcc:latest

RUN apt-get update && apt-get install -y \
    cmake \
    make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /rawhash2
COPY . /rawhash2

RUN mkdir build && cd build \
    && cmake .. \
    && make -j

ENTRYPOINT ["./build/bin/rawhash2"]
CMD ["--help"]

LABEL Name=rawhash2 Version=0.0.1
