#!/usr/bin/env python
import sys
import re
import json
from tqdm import tqdm
import pymongo as mongo

def cid_smiles(cid, sm):
    s = sm.find_one({'cpd':cid}, {'smiles':1})
    if s is not None:
        return s['smiles'] 

def ec_smiles(ec_q, ec, sm):
    ec_cids = ec.find_one({'ec':ec_q}, {'substrate_ids':1, 'product_ids':1, 'ec':1})
    if ec_cids is not None:
        if 'substrate_ids' in ec_cids and 'product_ids' in ec_cids:
            return {'ec':ec_q, 
                    'substrate_smiles':[cid_smiles(i, sm) for i in ec_cids['substrate_ids']],
                    'product_smiles':[cid_smiles(i, sm) for i in ec_cids['product_ids']],
                    }
        elif 'substrate_ids' in ec_cids:
            return {'ec':ec_q, 
                    'substrate_smiles':[cid_smiles(i, sm) for i in ec_cids['substrate_ids']],
                    }
        elif 'product_ids' in ec_cids:
            return {'ec':ec_q, 
                    'product_smiles':[cid_smiles(i, sm) for i in ec_cids['product_ids']],
                    }
        else:
            return {'ec':ec_q, 
                    'substrate_smiles':['unk'],
                    'product_smiles':['unk'],
                    }
    else:
        return {'ec':ec_q, 
                'substrate_smiles':['unk'],
                'product_smiles':['unk'],
                }

def main(args):
    client = mongo.MongoClient('mongodb://localhost:27017')
    db = client.get_database('uniprot')
    #p450s = db.get_collection('p450s')
    sprot = db.get_collection('sprot')

    ecdb = client.get_database('ec_nums')
    ec = ecdb.get_collection('ec')
    sm = ecdb.get_collection('smiles')


    if 'seq-smiles' not in db.list_collection_names():
        ssm = db.create_collection('seq-smiles')
    else:
        ssm = db.get_collection('seq-smiles')

    # trying to find
    todo = sprot.find({'_id':{'$nin':list(ssm.find({},{'_id':1}))}})

    lensprot = sprot.count_documents({})
    iteration = 0
    #with tqdm(sprot.find({}), mininterval=0.1) as bar:
    with tqdm(todo, mininterval=0.1) as bar:
        for doc in bar:
            iteration += 1
            if ssm.find_one({'_id':doc['_id']}) is not None:
                # '(EC=1\.14\.[0-9]+\.[0-9]+)' - p450s
                ec_nums =   re.findall('EC=([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)', 
                                json.dumps(\
                                     {i:j for i, j in zip(doc.keys(), doc.values()) if i != '_id'})
                                )
                for i in ec_nums:
                    o = {'_id':doc['_id'], 'seq':doc['seq'], **ec_smiles(i, ec, sm)}
                    ssm.update_one({'_id':o['_id']},
                                  {'$set':o},
                                  upsert=True)
                if iteration % 128 == 0:
                    a = ssm.count_documents({})
                    bar.set_postfix({"size seq-smiles":a,
                                    "of sprot": f"{round(a/lensprot * 100, 4)}%",
                                    })


if __name__ == '__main__':
    main(sys.argv[1:])
