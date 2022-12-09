FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ENV MPLBACKEND=Agg
ENV USE_CYTHON=TRUE
ARG PYTHON_VERSION
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda create -n testing -c conda-forge --yes python=$PYTHON_VERSION gxx_linux-aarch64"
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda env update -q -f ci/conda_host_env.yml"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda install -q --yes --file ci/aarch64.conda_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && pip install -r ci/aarch64.requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && pip install . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && conda list"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate testing && make test"
