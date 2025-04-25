#!/bin/bash
mkdir -p work

echo "UID=$(id -u)" > .env
echo "GID=$(id -g)" >> .env
echo "HOST_DIR=$PWD" >> .env
