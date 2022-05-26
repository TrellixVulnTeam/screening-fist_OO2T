#!/bin/sh
# debian
# sudo apt-get install gnupg
# sudo apt update && sudo apt install mongodb-org
# https://www.mongodb.com/docs/mongocli/stable/install/
wget -qO - https://www.mongodb.org/static/pgp/server-5.0.asc | sudo apt-key add -
echo "deb http://repo.mongodb.org/apt/debian buster/mongodb-org/5.0 main" | sudo tee /etc/apt/sources.list.d/mongodb-org-5.0.list
sudo apt update
sudo apt install -y mongocli

# https://stackoverflow.com/questions/64608581/mongodb-code-exited-status-14-failed-but-not-any-clear-errors
sudo chown mongodb:mongodb mongodb-27017.sock
sudo systemctl restart mongod.service
sudo systemctl status mongod.service
