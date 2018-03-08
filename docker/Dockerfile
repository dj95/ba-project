FROM sagemath/sagemath:8.0

#USER root

RUN sage -pip install colorama
#RUN apt-get install python-colorama

USER sage

WORKDIR /home/sage

COPY ./src/ /home/sage/.

ENTRYPOINT [ "sage-entrypoint" ]
