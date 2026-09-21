FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git clang-format unzip ca-certificates nlohmann-json3-dev python3 python3-numpy python3-numba \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /workspace
CMD ["bash", "run_regression.sh"]
