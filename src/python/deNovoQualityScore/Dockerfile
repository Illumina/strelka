FROM python:2.7

# install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends parallel

# update python package manager
RUN pip install -U pip setuptools

WORKDIR /opt/denovo

# install runtime dependencies (avoid cache invalidation)
COPY requirements.txt  /opt/denovo/
RUN pip install -U -r requirements.txt

# copy package content
COPY . .

# install test dependencies and run the tests
RUN pip install -U -r test/test_requirements.txt
RUN python -m pytest

# define bash as entry point
CMD ["/usr/bin/env", "bash"]
