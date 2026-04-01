# Build with docker buildx build -f Dockerfile -t spec:TAG --platform linux/amd64 . --load
FROM containers.qarnot.com/qrntadgaussf/roce.aocl:latest AS lib

# OCI compatible labels
LABEL org.opencontainers.image.created='2025-06-06T08:30:00.000Z'
LABEL org.opencontainers.image.authors='Samuel Lazerson (samuel.lazerson@gauss-fusion.com) Gauss Fusion GmbH'
LABEL org.opencontainers.image.url='https://containers.qarnot.com/samuelgaussf/spec'
LABEL org.opencontainers.image.documentation='Image of SPEC code for running on Qarnot machines.'
LABEL org.opencontainers.image.source='https://gitlab.gaussfusion.internal/system-codes/SPEC'
LABEL org.opencontainers.image.version='1.0'
LABEL org.opencontainers.image.revision='1'
LABEL org.opencontainers.image.vendor='Gauss Fusion GmbH'
LABEL org.opencontainers.image.licenses='MIT'
LABEL org.opencontainers.image.ref.name=''
LABEL org.opencontainers.image.title='SPEC code for Qarnot'
LABEL org.opencontainers.image.description='This image contains the SPEC software for MRxMHD equilibria.'

# Copy in hdf5
#COPY --from=containers.qarnot.com/qrntadgaussf/roce.hdf5:latest /usr/local/hdf5 /usr/local/hdf5

# Install m4
RUN apt-get update \
  && apt-get install -y \
  libhdf5-dev \
  m4 \
  python3 \
  python-is-python3 \
  python3-pip \
  python3-setuptools \
  python3-scikit-build-core

# Now start build
FROM lib AS build

# We add the HDF libraries to the path
#ENV LD_LIBRARY_PATH="/usr/local/hdf5/lib/:$LD_LIBRARY_PATH"

# Copy the app directory to the build
COPY . /app
WORKDIR /app

ENV CMAKE_ARGS_BLA_VENDOR="AOCL"
ENV FFTW3_ROOT_DIR="${AOCL_ROOT}"

RUN pip install --break-system-packages -v . 2>&1 | tee compile.log

# Now Check
RUN python -c "import spec; print(dir(spec))"

# Copy all built executables onto global path
RUN cp -RP /app/xspec /usr/local/bin/

