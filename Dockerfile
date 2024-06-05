FROM gcc:latest

RUN apt-get update && apt-get install -y \
    cmake make mold ccache \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /rawhash2
COPY . /rawhash2

RUN mkdir -p build && cd build \
    && cmake .. \
    && make -j 1

ENTRYPOINT ["./build/bin/rawhash2"]

LABEL Name=rawhash2 Version=0.0.1
