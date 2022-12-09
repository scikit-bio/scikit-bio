FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ENV MPLBACKEND=Agg
ENV USE_CYTHON=TRUE
ARG PYTHON_VERSION
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda create -n test -c conda-forge --yes python=$PYTHON_VERSION"
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && conda install -c conda-forge --yes --file ci/aarch64.conda_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && conda install --yes -c conda-forge gxx_linux-aarch64"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && pip install -r ci/aarch64.requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && pip install . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && make test"
