FROM mcr.microsoft.com/devcontainers/python:1-3.12-bookworm

USER root

ARG IGBLAST_VERSION=1.21.0
# The development base image includes a Yarn APT source.  This workshop does
# not use Yarn, and its rotating signing key can make unrelated APT installs
# fail in Codespaces, so remove only that source before updating package lists.
RUN grep -rl 'dl.yarnpkg.com' /etc/apt/sources.list /etc/apt/sources.list.d 2>/dev/null | xargs -r rm -f \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        default-jre-headless \
        fastqc \
        vsearch \
    && rm -rf /var/lib/apt/lists/*

# Pin IgBLAST so every Codespace starts from the same executable and baseline
# germline data.  The workshop's ferret files are overlaid below.
RUN curl --fail --location --silent --show-error \
        "https://ftp.ncbi.nlm.nih.gov/blast/executables/igblast/release/${IGBLAST_VERSION}/ncbi-igblast-${IGBLAST_VERSION}-x64-linux.tar.gz" \
    | tar --extract --gzip --directory /opt

ENV IGBLAST_HOME=/opt/ncbi-igblast-1.21.0
ENV PATH=${IGBLAST_HOME}/bin:${PATH}

COPY requirements-workshop.txt /tmp/requirements-workshop.txt
RUN python -m pip install --no-cache-dir -r /tmp/requirements-workshop.txt

# Keep only the custom ferret overlay in Git; it replaces/adds the reference
# database files distributed by the pinned upstream IgBLAST archive.
COPY igblast-overlay/internal_data/ferret/ ${IGBLAST_HOME}/internal_data/ferret/
COPY igblast-overlay/optional_file/ferret_gl.aux ${IGBLAST_HOME}/optional_file/ferret_gl.aux

USER vscode
