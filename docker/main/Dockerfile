FROM continuumio/anaconda3:latest

LABEL maintainer="Shane Brubaker <shanebrubaker55@gmail.com>"

# Fix: Remove asterisks from ENV variable
ENV DEBIAN_FRONTEND=noninteractive

# Install system packages and bioinformatics tools in a single layer
RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    # System packages
    git \
    wget \
    curl \
    build-essential \
    unzip \
    graphviz \
    libgraphviz-dev \
    pkg-config \
    swig \
    libx11-dev \
    libgsl0-dev \
    libopenblas-dev \
    liblapacke-dev \
    # Bioinformatics tools
    samtools \
    mafft \
    muscle \
    raxml \
    tabix \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean

# Clone the repository
RUN git clone https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-fourth-Edition.git

# Update conda and configure channels
RUN conda update -n base conda && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Create bioinformatics environment
RUN conda env create -f /Bioinformatics-with-Python-Cookbook-fourth-Edition/Ch01/bioinformatics_base.yml

# Install additional packages if needed
# RUN conda run -n bioinformatics_base pip install pyarrow==8.0.0

# Initialize conda for bash and set up environment activation
RUN conda init bash && \
    echo "conda activate bioinformatics_base" >> /root/.bashrc && \
    echo "setterm -foreground magenta" >> /etc/bash.bashrc

# Create workspace directory
RUN mkdir -p /workspace

# Set working directory
WORKDIR /Bioinformatics-with-Python-Cookbook-fourth-Edition

# Create a non-root user for security (optional but recommended)
# RUN useradd -m -s /bin/bash biouser && \
#     chown -R biouser:biouser /Bioinformatics-with-Python-Cookbook-fourth-Edition /workspace
# USER biouser

# Expose port
EXPOSE 9875

# Set environment variables for Jupyter
ENV JUPYTER_ENABLE_LAB=yes
ENV JUPYTER_TOKEN=""

# Start Jupyter Lab with better formatting
CMD ["conda", "run", "--no-capture-output", "-n", "bioinformatics_base", \
     "jupyter-lab", \
     "--ip=0.0.0.0", \
     "--no-browser", \
     "--allow-root", \
     "--port=9875", \
     "--NotebookApp.token=", \
     "--NotebookApp.password="]