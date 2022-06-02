# links:
# https://developer.nvidia.com/cuda-11.3.0-download-archive
# https://pytorch.org/
# https://www.linode.com/docs/products/compute/gpu/guides/install-nvidia-cuda/
mkdir src
cd src/

cd ~
mkdir src
cd src/

sudo apt update && sudo apt upgrade -y && sudo apt install git neovim gcc nvidia-cuda-toolkit linux-headers-$(uname -r)

# cuda 11.3
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2004/x86_64/cuda-ubuntu2004.pin
sudo mv cuda-ubuntu2004.pin /etc/apt/preferences.d/cuda-repository-pin-600
wget https://developer.download.nvidia.com/compute/cuda/11.3.0/local_installers/cuda-repo-ubuntu2004-11-3-local_11.3.0-465.19.01-1_amd64.deb
sudo dpkg -i cuda-repo-ubuntu2004-11-3-local_11.3.0-465.19.01-1_amd64.deb
sudo apt-key add /var/cuda-repo-ubuntu2004-11-3-local/7fa2af80.pub
sudo apt-get update
sudo apt-get -y install cuda nvtop

sudo reboot
