FROM sagemath/sagemath:7.6

USER root

RUN apt-get update 
RUN apt-get install python-colorama

USER sage

COPY ./src/ /home/sage/.

ENTRYPOINT [ "sage-entrypoint" ]
CMD [ "sage" ]
