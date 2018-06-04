FROM sagemath/sagemath:8.2

RUN sage -pip install colorama

WORKDIR /home/sage/project
COPY ./src/ /home/sage/project/.

ENTRYPOINT [ "sage-entrypoint" ]
