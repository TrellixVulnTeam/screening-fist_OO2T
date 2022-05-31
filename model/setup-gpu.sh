#!/bin/sh

sudo apt update && sudo apt upgrade -y && sudo apt install git neovim -y

# cuda
sudo apt install build-essential linux-headers-$(uname -r)
wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-debian11-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo dpkg -i cuda-repo-debian11-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo cp /var/cuda-repo-debian11-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo add-apt-repository contrib
sudo apt-get update
sudo apt-get -y install cuda

mkdir src
cd src/

git clone https://github.com/jamesengleback/dotfiles
cd dotfiles
./install-debian.sh
./setup-config.sh
./install-miniconda.sh
