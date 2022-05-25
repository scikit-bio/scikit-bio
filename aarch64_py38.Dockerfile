FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ENV MPLBACKEND=Agg
ENV USE_CYTHON=TRUE
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install -c conda-forge --yes python=3.8"
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install -c conda-forge --yes --file ci/conda_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install --yes -c conda-forge gxx_linux-aarch64"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && pip install -r ci/pip_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && pip install . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && make test"
