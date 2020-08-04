FROM python:3.7
LABEL maintainer "Ivo Leist"
WORKDIR /code
COPY requirements.txt /
RUN pip install -r /requirements.txt
COPY ./ ./
EXPOSE 8051
CMD ["python", "./usage.py"]
