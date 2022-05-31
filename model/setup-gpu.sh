#!/bin/bash
mkdir src
cd src/

sudo apt update && sudo apt upgrade -y && sudo apt install git neovim -y

# cuda
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.7.0/local_installers/cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-7-local_11.7.0-515.43.04-1_amd64.deb
sudo cp /var/cuda-repo-ubuntu2004-11-7-local/cuda-*-keyring.gpg /usr/share/keyrings/
sudo apt-get update
sudo apt-get -y install cuda


git clone https://github.com/jamesengleback/dotfiles
cd dotfiles
./install-debian.sh
./setup-config.sh
./install-miniconda.sh
source ~/.bashrc
source ~/miniconda3/etc/profile.d/conda.sh
conda update -n base -c defaults conda

echo 'PS1="\[\033[01;93m\] [gpu] \w \[\033[01;93m\]\$ \[\033[00m\]"' >> ~/.bashrc
cd ~/
git clone https://github.com/jamesengleback/screening-fist
cd screening-fist/model
#conda env create -f rdk-tch.yml
conda create -n tch -y
conda activate tch
conda install cuda pytorch rdkit -c nvidia -c pytorch -c rdkit -y
pip install wandb einops fair-esm tqdm
