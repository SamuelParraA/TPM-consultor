FROM python:3.12-slim
WORKDIR /app
COPY app.py config.example.json sample_metadata.csv ./
COPY web ./web
EXPOSE 8765
ENV TRANSCRIPTOMICA_ROOT=/data DB_PATH=/app/data/rnaseq_index.sqlite
CMD ["sh", "-c", "python app.py --host 0.0.0.0 --port ${PORT:-8765}"]
