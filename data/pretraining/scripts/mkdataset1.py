#!/usr/bin/env python
import sys
import re
import pymongo as mongo

def main(args):
    client = mongo.MongoClient('mongodb://localhost:27017')
    db = client.get_database('uniprot')

    bdb = db.get_collection('bindingdb')  # 865000
    ssm = db.get_collection('seq-smiles') # sprot p450s - 22989

    sys.stdout.write('seq,smiles,hit\n')
    for doc in bdb.find({}):
        sys.stdout.write(f"{','.join([doc['seq'], doc['smiles'], str(doc['hit'])])}\n")
    for doc in ssm.find({}):
        if re.search('[a-z,:.=\(\)]', doc['seq']) is None:
            for smiles in doc['substrate_smiles']:
                if smiles is not None:
                    sys.stdout.write(\
                            f"{','.join([doc['seq'], smiles, str(True)])}\n")


if __name__ == '__main__':
    main(sys.argv[1:])
