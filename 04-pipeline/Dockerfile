FROM alpine:3.15

RUN apk update && apk add R bash

ADD src /src
RUN chmod +x /src/*

WORKDIR /data

CMD ["/src/pipeline.sh"]