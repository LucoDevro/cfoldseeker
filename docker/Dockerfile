FROM vanallenlab/miniconda:3.12

LABEL description="cfoldseeker - Find gene clusters using protein structure similarity"
LABEL org.opencontainers.image.authors="lucas.devrieze@kuleuven.be"

# Use bash as shell
SHELL ["/bin/bash", "-c"]

# Install CAGEcleaner and add environment to PATH
RUN conda install -c bioconda -c conda-forge cfoldseeker -y
ENV PATH=/usr/miniconda3/envs/cfoldseeker/bin/:$PATH
