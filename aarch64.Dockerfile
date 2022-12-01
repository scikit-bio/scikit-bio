FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ENV MPLBACKEND=Agg
ENV USE_CYTHON=TRUE
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && conda install -c conda-forge --yes --file ci/conda_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && pip install -r ci/pip_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && pip install . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate test && make test"
