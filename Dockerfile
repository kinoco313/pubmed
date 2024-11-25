FROM python:3

WORKDIR /usr/dev

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
RUN pip install --upgrade deepl