# docker/Dockerfile
FROM condaforge/mambaforge:24.3.0-0

# where we’ll work inside the container
WORKDIR /work

# (optional) tiny OS deps; keep lean
RUN apt-get update && apt-get install -y --no-install-recommends \
    git && rm -rf /var/lib/apt/lists/*

# Ensure strict priority, add channels
SHELL ["/bin/bash", "-lc"]

# copy env file first (so Docker caches it if unchanged)
COPY envs/avml.yml /tmp/env.yml

# create env from file
RUN mamba env create -f /tmp/env.yml && mamba clean -afy

# use that env by default
SHELL ["/bin/bash", "-lc"]
ENV PATH=/opt/conda/envs/avml/bin:$PATH

# copy repo last (so code changes don’t force env rebuild)
COPY . /work

# default command if you just run the image
ENTRYPOINT ["snakemake"]
CMD ["--cores", "4"]
