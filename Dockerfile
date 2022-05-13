FROM osgeo/gdal:alpine-small-3.4.3 as build
RUN \
    apk update && \
    apk add \
    boost-program_options \
    boost-dev \
    g++ \
    make \
    cmake

RUN mkdir /app
COPY * /app
WORKDIR /app
RUN cmake -S . -B build && cd build && make

FROM osgeo/gdal:alpine-normal-3.4.3 as result
ENV PATH="/app:${PATH}"
RUN \
    apk update && \
    apk add \
    bash \
    boost-program_options

RUN mkdir /app
WORKDIR /app
COPY --from=build /app/build/geosheightcorrection /app
COPY --from=build /app/perform_geos_height_correction.sh /app
RUN chmod +x /app/*
CMD geosheightcorrection
