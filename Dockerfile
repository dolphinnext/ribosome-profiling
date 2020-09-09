# Don't upgrade nfcore/base, it creates "Kernel too old" error for singularity (because of the debian image)
FROM nfcore/base:1.7
LABEL author="onur.yukselen@umassmed.edu" description="Docker image containing all requirements for the dolphinnext/riboseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN apt-get update && apt-get install -y sed gcc zlib1g-dev g++ tree
# Install dolphin-tools
RUN git clone https://github.com/dolphinnext/tools /usr/local/bin/dolphin-tools
RUN git clone https://github.com/jrw24/G418_readthrough.git /usr/local/bin/G418_readthrough
RUN sed -i '#!/usr/bin/env python' /usr/local/bin/G418_readthrough/riboseq/densebuilder_main.py
RUN sed -i -e 's/if not chrom in validChrs/#if not chrom in validChrs/g' /usr/local/bin/G418_readthrough/riboseq/densebuilder_main.py
ENV PATH /opt/conda/envs/dolphinnext-riboseq-1.0/bin:/usr/local/bin/dolphin-tools/:/usr/local/bin/G418_readthrough/riboseq:/usr/local/bin/G418_readthrough/utils:/usr/local/bin/G418_readthrough/RNAseq:$PATH

RUN mv /usr/local/bin/G418_readthrough/GFF /opt/conda/envs/dolphinnext-riboseq-1.0/lib/python2.7/site-packages/.

RUN mkdir -p /project /nl /mnt /share

RUN pip install 'biopython==1.58'
RUN python -c "import Bio.Seq"
RUN python -c "import GFF"

### SRA-toolkit
RUN cd /usr/bin && wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.7/sratoolkit.2.10.7-ubuntu64.tar.gz && \
    tar -xvzf sratoolkit.2.10.7-ubuntu64.tar.gz
ENV PATH $PATH:/usr/bin/sratoolkit.2.10.7-ubuntu64/bin

### Kraken tally
RUN cd /usr/bin && wget http://wwwdev.ebi.ac.uk/enright-dev/kraken/reaper/binaries/reaper-13-100/linux/tally && chmod 755 tally

### Seqtk
RUN cd /usr/bin && git clone https://github.com/lh3/seqtk.git && cd seqtk && make
ENV PATH $PATH:/usr/bin/seqtk

### kentUtils
RUN cd /usr/bin && wget https://github.com/ENCODE-DCC/kentUtils.git
ENV PATH $PATH:/usr/bin/kentUtils/bin/linux.x86_64

### skewer
RUN cd /usr/bin && git clone https://github.com/relipmoc/skewer.git && cd skewer && make && make install
ENV PATH $PATH:/usr/bin/skewer

### pigz
RUN cd /usr/bin && wget https://zlib.net/pigz/pigz-2.4.tar.gz && tar xvzf pigz-2.4.tar.gz && cd pigz-2.4 && make
ENV PATH $PATH:/usr/bin/pigz-2.4