FROM mambaorg/micromamba:1.4.9
LABEL maintainer="Remi Denise (rdenise@ucc.ie)"
LABEL version="0.2.7"

RUN micromamba install -y -n base -c conda-forge -c bioconda taxmyphage==0.2.7 && \
    micromamba clean --all --yes
WORKDIR /app
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "taxmyphage"]