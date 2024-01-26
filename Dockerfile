FROM debian:bullseye-20240110-slim as build

RUN apt-get update && \
    apt-get install -y \
    pkg-config \
    cmake \
    make \
    libboost-program-options-dev \
    libfftw3-dev \
    libgsl-dev \
    libgomp1 \
    clang \
    git

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /deps

RUN git clone https://github.com/jkriege2/TinyTIFF && \
    cd TinyTIFF && \
    mkdir build && \
    cd build && \
    cmake .. -DTinyTIFF_BUILD_SHARED_LIBS=ON && \
    cmake --build . --target install 

WORKDIR /build

COPY . .

RUN mkdir build  && \
    make x86_64build &&  \
    mv main /bin/ddm 

FROM busybox:latest as run

WORKDIR /files

COPY --from=build /bin/ddm /bin/ddm
COPY --from=build /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.74.0 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libomp.so.5 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libgomp.so.1 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libstdc++.so.6 /usr/lib/x86_64-linux-gnu/
COPY --from=build /lib/x86_64-linux-gnu/libgcc_s.so.1 /lib/x86_64-linux-gnu/
COPY --from=build /lib/x86_64-linux-gnu/libdl.so.2 /lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libgsl.so.25 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libgslcblas.so.0 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libfftw3f_omp.so.3 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libfftw3.so.3 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/lib/x86_64-linux-gnu/libfftw3f.so.3 /usr/lib/x86_64-linux-gnu/
COPY --from=build /usr/local/lib/libTinyTIFFShared_Release.so.3.0.0.0 /usr/local/lib/

ENV LD_LIBRARY_PATH=/usr/local/lib

CMD ["ddm"]