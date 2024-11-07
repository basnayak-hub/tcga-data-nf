# bring in the micromamba image so we can copy files from it
FROM mambaorg/micromamba:1.5.1 as micromamba

# This is the image we are going add micromaba to:
FROM violafanfani/tcga-data-base-r:0.0.4

USER root

# adding opencontainer labels to link registry to github repository
LABEL org.opencontainers.image.title="tcga-data-nf"
LABEL org.opencontainers.image.description="Workflow to download and prepare TCGA data"
LABEL org.opencontainers.image.url="https://hub.docker.com/r/violafanfani/tcga-data-nf"
LABEL org.opencontainers.image.documentation="https://github.com/QuackenbushLab/tcga-data-nf"
LABEL org.opencontainers.image.source="https://github.com/QuackenbushLab/tcga-data-nf"
LABEL org.opencontainers.image.vendor="https://github.com/QuackenbushLab"
LABEL org.opencontainers.image.authors="Viola Fanfani"
LABEL org.opencontainers.image.revision="v0.0.14"

# if your image defaults to a non-root user, then you may want to make the
# next 3 ARG commands match the values in your image. You can get the values
# by running: docker run --rm -it my/image id -a
ARG MAMBA_USER=root
ARG MAMBA_USER_ID=0
ARG MAMBA_USER_GID=0
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_ROOT_PREFIX="/opt/conda"
ENV MAMBA_EXE="/bin/micromamba"

COPY --from=micromamba "$MAMBA_EXE" "$MAMBA_EXE"
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY --from=micromamba /usr/local/bin/_entrypoint.sh /usr/local/bin/_entrypoint.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh || (echo "Error in _dockerfile_initialize_user_accounts.sh" && false)
RUN /usr/local/bin/_dockerfile_setup_root_prefix.sh || (echo "Error in _dockerfile_setup_root_prefix.sh" && false)

USER $MAMBA_USER

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# Optional: if you want to customize the ENTRYPOINT and have a conda
# environment activated, then do this:
# ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "my_entrypoint_program"]

# You can modify the CMD statement as needed....
CMD ["/bin/bash"]

#COPY --chown=micromamba:micromamba ${CONDA_FILE} /tmp
ARG CONDA_FILE=containers/env.base.python.yml
COPY --chown=micromamba:micromamba ${CONDA_FILE} /tmp
# Optional: you can now populate a conda environment:
RUN micromamba install -y -n base -f /tmp/`basename ${CONDA_FILE}` && \
    micromamba clean --all --yes

USER root

RUN R -e "devtools::install_github('pmandros/TCGAPurityFiltering')"
RUN R -e "devtools::install_github('immunogenomics/presto')"
RUN R -e "devtools::install_github('aet21/EpiSCORE')"

#ADD  NetSciDataCompanion/ /opt/NetSciDataCompanion/
#RUN R -e "devtools::install_local('/opt/NetSciDataCompanion/', force=TRUE, dependencies=T)"
RUN R -e "devtools::install_github('https://github.com/QuackenbushLab/NetworkDataCompanion')"