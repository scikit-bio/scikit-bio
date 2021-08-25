FROM --platform=linux/arm64 condaforge/linux-anvil-aarch64
RUN sudo yum update -y && \
	sudo yum install -y make git && \
	sudo yum clean all
ENV MPLBACKEND=Agg
ENV USE_CYTHON=TRUE
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install --yes python=3.6"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install --yes -c conda-forge gxx_linux-aarch64=7.5.0"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && \
	conda install --yes -c conda-forge pandas=1.1.5 scipy=1.5.3 \
	ipython=7.16.1 matplotlib=3.3.3"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install --yes cython pip"
COPY . /work
WORKDIR /work
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && conda install --yes --file ci/conda_requirements.txt"
# ^^ Can't do this directly as many of the latest versions are missing a python3.6 build for aarch64
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && pip install -r ci/pip_requirements.txt"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && pip install . --no-deps"
RUN bash -c ". /opt/conda/etc/profile.d/conda.sh && conda activate base && make test"
