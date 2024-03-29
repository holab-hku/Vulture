FROM r-base:4.1.0
#  $ docker build . -t continuumio/anaconda3:latest -t continuumio/anaconda3:5.3.0
#  $ docker run --rm -it continuumio/anaconda3:latest /bin/bash
#  $ docker push continuumio/anaconda3:latest
#  $ docker push continuumio/anaconda3:5.3.0
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
ENV DEBIAN_FRONTEND=noninteractive 
ENV PATH=/opt/samtools/bin:$PATH
ENV PATH=/STAR-2.7.9a/bin/Linux_x86_64_static:$PATH
ENV PATH="/opt/sratoolkit.2.11.0-ubuntu64/bin/:${PATH}"

COPY . /code

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates \
    libglib2.0-0 libxext6 libsm6 libxrender1 curl grep sed dpkg libcurl4-openssl-dev libssl-dev libhdf5-dev \
    git mercurial subversion procps build-essential libc6-dev \
    libxml-libxml-perl pigz awscli uuid-runtime time tini

RUN Rscript /code/r/scvh_dependencies.r

RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz -O /tmp/sratoolkit.tar.gz \
	&& tar zxvf /tmp/sratoolkit.tar.gz -C /opt/ && rm /tmp/sratoolkit.tar.gz
RUN mkdir -p /root/.ncbi
RUN printf '/LIBS/GUID = "%s"\n' `uuidgen` > /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/file/volumes/flatAd = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/refseq/volumes/refseqAd = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/sra/volumes/sraAd = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/sraPileup/volumes/ad = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/sraRealign/volumes/ad = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/apps/wgs/volumes/wgsAd = "."\n' >> /root/.ncbi/user-settings.mkfg
RUN printf '/repository/user/ad/public/root = "."\n' >> /root/.ncbi/user-settings.mkfg

RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh -O ~/anaconda.sh && \
    /bin/bash ~/anaconda.sh -b -p /opt/conda && \
    rm ~/anaconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc && \
    echo "umask 000" >> ~/.bashrc && \
    echo "ulimit -n 4096" >> ~/.bashrc && \
    # apt-get install -y curl grep sed dpkg && \
    # TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    # curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    # dpkg -i tini.deb && \
    # rm tini.deb && \
    apt-get clean && \
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -vxjf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make && \
    cd ../ \
    pip install kb-python==0.27.3 pysam==0.19.1 umi_tools==1.1.2 velocyto==0.17.17 scvelo==0.2.4 boto3==1.24.91 awscli==1.25.92 leidenalg==0.9.0 bbknn==1.5.1 && \
    wget https://github.com/alexdobin/STAR/archive/2.7.9a.tar.gz && \
    tar -xzf 2.7.9a.tar.gz && \
    echo "export PATH=/STAR-2.7.9a/bin/Linux_x86_64_static:\$PATH" >> ~/.bashrc && \
    wget https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xf samtools-1.13.tar.bz2 && \
    cd samtools-1.13 && \
    ./configure --prefix=/opt/samtools && \
    make && \
    make install && \
    echo "export PATH=/opt/samtools/bin:\$PATH" >> ~/.bashrc && \
    . ~/.bashrc && \
    chmod 755 /code
ENTRYPOINT [ "/usr/bin/tini","--"]
CMD [ "/bin/bash" ]
