FROM ghcr.io/astral-sh/uv:python3.12-bookworm-slim AS builder
WORKDIR /app
RUN useradd --create-home --shell /bin/bash appuser
COPY --chown=appuser:appuser api/requirements.txt .
RUN --mount=type=cache,target=/root/.cache/uv \
    uv venv && \
    uv pip install -r requirements.txt
COPY --chown=appuser:appuser api/ ./api/
COPY --chown=appuser:appuser src_py/ ./src_py/
COPY --chown=appuser:appuser assets/ ./assets/
FROM python:3.12-slim-bookworm AS final
WORKDIR /app
RUN useradd --create-home --shell /bin/bash appuser
COPY --from=builder --chown=appuser:appuser /app .
USER appuser
EXPOSE 8000
CMD ["/app/.venv/bin/python", "-m", "uvicorn", "api.index:app", "--host", "0.0.0.0", "--port", "8000"]
