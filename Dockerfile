# =============================================================================
# SURFR Pipeline — Dockerfile (linux/amd64)
#
# Built for x86-64 (Dardel / Intel HPC). Build on Apple Silicon with:
#   docker buildx build --platform linux/amd64 --tag surfr_pipeline:latest --load .
#
# Tool versions:
#   samtools        1.23.1 (preprint used 1.20; merge/fastq output unaffected,
#                           see README "Notes on the preprint text")
#   KMC             3.2.4
#   miRTrace        1.0.1  (Java 21 runtime)
#   R               4.4.1  (rocker/r-ver:4.4.1, matches the preprint)
#   R packages      tidyverse, paletteer, arrow, ggvenn, ggrastr, MASS, installed from
#                   the Posit Package Manager CRAN snapshot of 2024-10-30
#                   (last day R 4.4.1 was the current release). Exact installed
#                   versions are written to /opt/surfr/R_package_versions.tsv.
#   dekupl-mergeTags commit 4cdad2c (2017-06-27, last commit of the repository)
#   pigz            system (Ubuntu 22.04 apt)
# =============================================================================

# rocker/r-ver:4.4.1 is Ubuntu 22.04 (jammy) with R 4.4.1 built from source.
# For full immutability, replace the tag with its digest, obtained with:
#   docker buildx imagetools inspect rocker/r-ver:4.4.1
FROM --platform=linux/amd64 rocker/r-ver:4.4.1

# ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------
LABEL maintainer="SURFR-pipeline"
LABEL version="1.0.0"
LABEL samtools="1.23.1"
LABEL KMC="3.2.4"
LABEL miRTrace="1.0.1"
LABEL dekupl-mergeTags="4cdad2c5ce45c3a30458aa73ce970e31c7646699"
LABEL R="4.4.1"
LABEL cran_snapshot="2024-10-30"

# ---------------------------------------------------------------------------
# Environment — set once, available in every subsequent RUN and at runtime
# CRAN: pinned snapshot used by install_r_packages.R. Set explicitly here so
#       the pin does not depend on what the base image happens to configure.
# ---------------------------------------------------------------------------
ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    JAVA_TOOL_OPTIONS="-Djava.awt.headless=true" \
    CRAN="https://p3m.dev/cran/__linux__/jammy/2024-10-30" \
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
        libcairo2-dev \
        libxt-dev \
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
# 6. R 4.4.1 — provided by the rocker/r-ver:4.4.1 base image
#    (R_HOME=/usr/local/lib/R, Rscript at /usr/local/bin/Rscript).
#    Do NOT apt-install r-base here: it would pull the latest CRAN R.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 7. R packages — installed from the pinned snapshot into rocker's default
#    site library. Rscript is run WITHOUT --vanilla on purpose: --vanilla
#    skips Rprofile.site/Renviron.site, where rocker configures the library
#    path and the user agent needed for Package Manager binaries.
# ---------------------------------------------------------------------------
RUN mkdir -p /opt/surfr
COPY install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && \
    rm /tmp/install_r_packages.R

# Enforced R check, kept in its own RUN so the '|| true' chain below
# cannot swallow a failure.
RUN Rscript -e "stopifnot(getRversion() == '4.4.1'); for (p in c('tidyverse','paletteer','arrow','ggvenn','ggrastr','MASS')) { suppressPackageStartupMessages(library(p, character.only = TRUE)); cat(p, as.character(packageVersion(p)), 'OK\n') }"

# ---------------------------------------------------------------------------
# 8. dekupl-mergeTags — built from source, pinned to a fixed commit
# ---------------------------------------------------------------------------
RUN MERGETAGS_COMMIT=4cdad2c5ce45c3a30458aa73ce970e31c7646699 && \
    mkdir -p /opt/dekupl/bin && \
    git clone https://github.com/Transipedia/dekupl-mergeTags.git /tmp/dekupl-mergeTags && \
    cd /tmp/dekupl-mergeTags && \
    git checkout --quiet "${MERGETAGS_COMMIT}" && \
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
    /opt/dekupl/bin/mergeTags 2>&1 | head -3 || true && \
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
echo "  Rscript          /usr/local/bin/Rscript (R 4.4.1)"\n\
echo "  dekupl-mergeTags /opt/dekupl/bin/mergeTags"\n\
echo "  R package list   /opt/surfr/R_package_versions.tsv"\n\
echo ""\n\
echo "Convert to a Singularity image (SIF) for Dardel:"\n\
echo "  docker save surfr_pipeline:latest -o surfr_pipeline_docker.tar   (on your machine)"\n\
echo "  copy the .tar to the cluster, then on the cluster:"\n\
echo "  singularity build surfr_pipeline.sif docker-archive://\$PWD/surfr_pipeline_docker.tar"\n\
echo "  run with: SINGULARITY_TMPDIR=/tmp singularity exec --no-mount bind-paths -B /cfs/klemming surfr_pipeline.sif ..."\n\
echo "  (--no-mount bind-paths is required on Dardel; see README, Container)"\n\
' > /entrypoint.sh && chmod +x /entrypoint.sh

CMD ["/bin/bash", "/entrypoint.sh"]
