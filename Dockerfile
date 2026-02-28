# Stage 1: Test — install with dev deps and run pytest
FROM python:3.11-slim AS test

WORKDIR /app

COPY pyproject.toml .
COPY src/ src/
COPY tests/ tests/
COPY examples/ examples/

RUN pip install --no-cache-dir -e ".[dev]"
RUN pytest -v

# Stage 2: Production — clean image with runtime deps only
FROM python:3.11-slim AS production

RUN apt-get update && \
    apt-get install -y --no-install-recommends procps && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY pyproject.toml .
COPY src/ src/
COPY examples/ examples/

RUN pip install --no-cache-dir .

ENTRYPOINT ["/bin/bash"]
