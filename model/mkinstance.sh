#!/bin/sh
linode-cli linodes create --type 'g1-gpu-rtx6000-1'  --root_pass $1  --region 'us-east'  --label 'gpu'  --image 'linode/ubuntu22.04'  --authorized_keys 'ssh-rsa'
