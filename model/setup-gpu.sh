#!/bin/bash
cd ~
mkdir src
cd src/

sudo apt update && sudo apt upgrade -y && sudo apt install git neovim gcc nvidia-cuda-toolkit linux-headers-$(uname -r)

# cuda
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.3.0/local_installers/cuda-repo-ubuntu2004-11-3-local_11.3.0-465.19.01-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-3-local_11.3.0-465.19.01-1_amd64.deb
sudo apt-key add /var/cuda-repo-ubuntu2004-11-3-local/7fa2af80.pub
sudo apt-get update
sudo apt-get -y install cuda

sudo apt install nvtop

## reboot here

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
conda create -n nn -y
conda activate nn
conda install pytorch torchvision torchaudio cudatoolkit=11.3 rdkit -c pytorch -c rdkit -y
pip install wandb einops fair-esm tqdm ipython
