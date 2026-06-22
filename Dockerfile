# =============================================================================
# SURFR Pipeline — Dockerfile (linux/amd64)
#
# Built for x86-64 (Dardel / Intel HPC). Build on Apple Silicon with:
#   docker buildx build --platform linux/amd64 --tag surfr_pipeline:latest --load .
#
# Tool versions:
#   samtools        1.23.1
#   KMC             3.2.4
#   miRTrace        1.0.1  (Java 21 runtime)
#   R               4.4    + tidyverse, paletteer, arrow, ggvenn, MASS
#   dekupl-run      1.3.5  (mergeTags binary only)
#   pigz            system (Ubuntu 22.04 apt)
# =============================================================================

FROM --platform=linux/amd64 ubuntu:22.04

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
LABEL maintainer="SURFR-pipeline"
LABEL version="1.0.0"
LABEL samtools="1.23.1"
LABEL KMC="3.2.4"
LABEL miRTrace="1.0.1"
LABEL dekupl-run="1.3.5"
LABEL R="4.4"

# ---------------------------------------------------------------------------
# Environment — set once, available in every subsequent RUN and at runtime
# ---------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" \
    R_LIBS_SITE=/opt/R/library \
    PATH=/opt/kmc/bin:/opt/mirtrace:/opt/dekupl/bin:/usr/local/bin:$PATH

# ---------------------------------------------------------------------------
# 1. Base system packages
# ---------------------------------------------------------------------------
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        autoconf \
        automake \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        dirmngr \
        git \
        gnupg2 \
        libbz2-dev \
        libgomp1 \
        libcurl4-openssl-dev \
        libgsl-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        libzstd-dev \
        make \
        perl \
        pigz \
        python3 \
        python3-pip \
        python3-psutil \
        software-properties-common \
        unzip \
        wget \
        zlib1g-dev \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff5-dev \
        libuv1-dev \
        libwebp-dev \
        libxml2-dev \
        pkg-config && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 2. Java 21 (Eclipse Temurin) — required by miRTrace
# ---------------------------------------------------------------------------
RUN wget -qO - https://packages.adoptium.net/artifactory/api/gpg/key/public \
        | gpg --dearmor -o /etc/apt/trusted.gpg.d/adoptium.gpg && \
    echo "deb https://packages.adoptium.net/artifactory/deb jammy main" \
        > /etc/apt/sources.list.d/adoptium.list && \
    apt-get update -qq && \
    apt-get install -y --no-install-recommends temurin-21-jdk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 3. samtools 1.23.1 — built from source
# ---------------------------------------------------------------------------
RUN SAMTOOLS_VERSION=1.23.1 && \
    wget -q "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
        -O /tmp/samtools.tar.bz2 && \
    tar -xjf /tmp/samtools.tar.bz2 -C /tmp && \
    cd /tmp/samtools-${SAMTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --without-curses && \
    make -j$(nproc) && \
    make install && \
    rm -rf /tmp/samtools*

# ---------------------------------------------------------------------------
# 4. KMC 3.2.4 — built from source via git clone (submodules required)
#    The GitHub source tarball omits the cloudflare/zlib submodule that
#    KMC's Makefile needs. A full clone with --recurse-submodules is required.
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/kmc/bin && \
    git clone --depth 1 --branch v3.2.4 --recurse-submodules \
        https://github.com/refresh-bio/KMC.git /tmp/KMC && \
    cd /tmp/KMC && \
    make kmc kmc_tools && \
    mv bin/kmc       /opt/kmc/bin/kmc && \
    mv bin/kmc_tools /opt/kmc/bin/kmc_tools && \
    chmod +x /opt/kmc/bin/kmc /opt/kmc/bin/kmc_tools && \
    rm -rf /tmp/KMC

# ---------------------------------------------------------------------------
# 5. miRTrace 1.0.1 — JAR + Python wrapper
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/mirtrace && \
    wget -q "https://github.com/friedlanderlab/mirtrace/releases/download/v1.0.1/mirtrace-v1.0.1.zip" \
        -O /tmp/mirtrace.zip && \
    unzip -q /tmp/mirtrace.zip -d /tmp/mirtrace_extract && \
    MIRTRACE_DIR=$(find /tmp/mirtrace_extract -maxdepth 1 -mindepth 1 -type d | head -1) && \
    cp "${MIRTRACE_DIR}/mirtrace.jar" /opt/mirtrace/mirtrace.jar && \
    cp "${MIRTRACE_DIR}/mirtrace"     /opt/mirtrace/mirtrace && \
    chmod +x /opt/mirtrace/mirtrace && \
    sed -i '1s|#!/usr/bin/env python$|#!/usr/bin/env python3|' /opt/mirtrace/mirtrace && \
    rm -rf /tmp/mirtrace*

# ---------------------------------------------------------------------------
# 6. R 4.4 — installed from CRAN apt repo (Posit/CRAN official)
# ---------------------------------------------------------------------------
RUN apt-get update -qq && \
    apt-get install -y --no-install-recommends \
        software-properties-common \
        dirmngr \
        apt-transport-https && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc \
        | gpg --dearmor -o /usr/share/keyrings/r-project.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/r-project.gpg] https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/" \
        > /etc/apt/sources.list.d/r-project.list && \
    apt-get update -qq && \
    apt-get install -y --no-install-recommends r-base r-base-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# ---------------------------------------------------------------------------
# 7. R packages — write install script to file to avoid shell-escaping issues
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/R/library && \
    echo 'R_LIBS_SITE=/opt/R/library' >> /etc/R/Renviron
COPY install_r_packages.R /tmp/install_r_packages.R
RUN Rscript --vanilla /tmp/install_r_packages.R && \
    rm /tmp/install_r_packages.R

# ---------------------------------------------------------------------------
# 8. dekupl-mergeTags 1.3.5 — build from source
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/dekupl/bin && \
    git clone --depth 1 https://github.com/Transipedia/dekupl-mergeTags.git /tmp/dekupl-mergeTags && \
    cd /tmp/dekupl-mergeTags && \
    make && \
    mv mergeTags /opt/dekupl/bin/mergeTags && \
    chmod +x /opt/dekupl/bin/mergeTags && \
    rm -rf /tmp/dekupl-mergeTags

# ---------------------------------------------------------------------------
# 9. Smoke tests — build fails here if any tool is broken
# ---------------------------------------------------------------------------
RUN echo "=== Smoke tests ===" && \
    samtools --version | head -1 && \
    pigz --version 2>&1 | head -1 && \
    /opt/kmc/bin/kmc     2>&1 | head -1 || true && \
    /opt/kmc/bin/kmc_tools 2>&1 | head -1 || true && \
    java -version 2>&1 && \
    java -jar /opt/mirtrace/mirtrace.jar --version 2>&1 | head -1 && \
    Rscript --version && \
    Rscript --vanilla -e "for(p in c('tidyverse','paletteer','arrow','ggvenn','MASS')){library(p,character.only=TRUE,lib.loc='/opt/R/library');cat(p,'OK\n')}" && \
    /opt/dekupl/bin/mergeTags --help 2>&1 | head -3 || true && \
    echo "=== All smoke tests passed ==="

# ---------------------------------------------------------------------------
# 10. Default entrypoint — prints a usage summary (mirrors %runscript)
# ---------------------------------------------------------------------------
RUN printf '#!/bin/bash\n\
echo ""\n\
echo "SURFR Pipeline Container"\n\
echo "========================"\n\
echo ""\n\
echo "Tool paths:"\n\
echo "  samtools         /usr/local/bin/samtools"\n\
echo "  pigz             /usr/bin/pigz"\n\
echo "  kmc              /opt/kmc/bin/kmc"\n\
echo "  kmc_tools        /opt/kmc/bin/kmc_tools"\n\
echo "  mirtrace         /opt/mirtrace/mirtrace"\n\
echo "  Rscript          /usr/bin/Rscript"\n\
echo "  dekupl-mergeTags /opt/dekupl/bin/mergeTags"\n\
echo ""\n\
echo "Convert to Singularity sandbox for Dardel:"\n\
echo "  docker save surfr_pipeline:latest -o surfr_pipeline_docker.tar"\n\
echo "  docker run --rm --privileged --platform linux/amd64 \\\\"\n\
echo "    -v \$(pwd):/work quay.io/singularity/singularity:v4.2.0 \\\\"\n\
echo "    build --sandbox /work/surfr_pipeline/ docker-archive:///work/surfr_pipeline_docker.tar"\n\
' > /entrypoint.sh && chmod +x /entrypoint.sh

CMD ["/bin/bash", "/entrypoint.sh"]
