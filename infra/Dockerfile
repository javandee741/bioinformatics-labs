# infra/Dockerfile
FROM python:3.11-slim

# 1) системные пакеты по минимуму
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl ca-certificates build-essential && \
    rm -rf /var/lib/apt/lists/*

# 2) не-root пользователь (безопаснее)
RUN useradd -ms /bin/bash appuser
WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 3) только нужный код
COPY scripts/ scripts/
COPY workflows/ workflows/
#COPY Makefile README.md ./

# 4) права
RUN chown -R appuser:appuser /app
USER appuser

# 5) default CMD можно оставить пустым — запускать через docker run/compose
#CMD ["bash", "-lc", "python scripts/gc_content.ipynb -i data/example.fasta -o results/revcomp.fasta"]
CMD ["bash", "-c", "echo 'Hello from bioinfo:week4 container!' > /app/results/test_output.txt"]
