FROM osgeo/gdal:ubuntu-small-3.6.3 as build
RUN \
    apt-get update && \
    apt-get install -y \
    libboost-program-options1.74-dev \
    libboost-test1.74-dev \
    libdlib-dev \
    libblas-dev \
    g++ \
    make \
    cmake


RUN mkdir /app
COPY ./ /app
WORKDIR /app
RUN cmake -DCMAKE_BUILD_TYPE=Release -S . -B build && cd build && make -j VERBOSE=1
RUN cd build && make test VERBOSE=1

FROM osgeo/gdal:ubuntu-small-3.6.3 as result
ENV PATH="/app:${PATH}"
RUN \
    apt-get update && \
    apt-get install -y \
    libboost-program-options1.74.0 \
    libdlib19 \
    libblas3 \
    libgomp1

RUN mkdir /app
WORKDIR /app
COPY --from=build /app/build/geosheightcorrection /app
COPY --from=build /app/build/geostablegenerator /app
COPY --from=build /app/perform_geos_height_correction.sh /app
RUN chmod +x /app/*
CMD geosheightcorrection
