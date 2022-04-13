#!/bin/bash

CID=$(docker run --cap-add=SYS_PTRACE --security-opt seccomp=unconfined -v /Users/horta/code/hmmer-pick-profile:/hmmer-pick-profile -w /hmmer-pick-profile --platform linux/amd64 -it --detach --name hmmer ubuntu:21.10)

TZ=Europe/London
ENV="-e TZ=$TZ -e DEBIAN_FRONTEND=noninteractive -e DEBCONF_TERSE=yes"
APT="apt-get -y -q --no-install-recommends"

docker exec "$ENV" "$CID" $APT update
docker exec "$ENV" "$CID" $APT install tzdata
docker exec "$ENV" "$CID" ln -snf /usr/share/zoneinfo/$TZ /etc/localtime
docker exec "$ENV" "$CID" dpkg-reconfigure --frontend noninteractive tzdata

docker exec "$ENV" "$CID" $APT upgrade
docker exec "$ENV" "$CID" $APT install apt-utils
docker exec "$ENV" "$CID" $APT install ca-certificates
docker exec "$ENV" "$CID" $APT install curl
docker exec "$ENV" "$CID" $APT install build-essential
docker exec "$ENV" -it "$CID" /bin/bash
docker stop "$CID" >/dev/null 2>&1 &
disown

echo "Enter: docker exec $ENV -it $CID /bin/bash"
