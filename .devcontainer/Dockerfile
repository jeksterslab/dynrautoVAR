FROM jeksterslab/dynr-rocker

# add additional packages
ADD tools /usr/src/local/src
RUN chmod 777 /usr/src/local/src/setup.sh && \
    ./usr/src/local/src/setup.sh          && \
    rm -rf /usr/src/local/src

# add user with sudo priviledges
RUN apt-get update && apt-get -y install sudo
RUN useradd -ms /bin/bash user && echo "user:user" | chpasswd && adduser user sudo
USER user
WORKDIR /home/user

# create dotfiles
ADD tools /home/user/src
RUN chmod 777 /home/user/src/dotfiles.sh && \
    ./home/user/src/dotfiles.sh          && \
    rm -rf /home/user/src

# extra metadata
LABEL author="Ivan Jacob Agaloos Pesigan <learn.jeksterslab@gmail.com>"
LABEL description="jeksterslab/dynr-rocker container."
