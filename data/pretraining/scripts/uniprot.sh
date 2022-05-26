#!/bin/sh
# big one:

MONGODB_ADDR="mongodb://localhost:27017/uniprot"

#curl -s https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.dat.gz | gunzip -c -d | ./uniprot_to_json.py | mongoimport -v $MONGODB_ADDR -c trembl
curl -s https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz | gunzip -c -d | ./uniprot_to_json.py | mongoimport -v $MONGODB_ADDR -c sprot
#curl -s https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz | gunzip -c -d | ./uniprot_to_json.py 
