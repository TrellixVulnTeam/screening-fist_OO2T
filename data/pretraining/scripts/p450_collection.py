#!/usr/bin/env python
import sys
import re
import pymongo as mongo

def main(args):  
    MONGODB_ADDR="mongodb://localhost:27017/"
    DBNAME = 'uniprot'
    client = mongo.MongoClient(MONGODB_ADDR)
    assert DBNAME in client.list_database_names()
    db = client.get_database(DBNAME)

    sprot = db['sprot']

    qr_p450 = sprot.find({'$or': [
                          {'CC':{'$regex':'.*[Pp]450.*'}},
                          {'KW':{'$regex':'.*[Pp]450.*'}},
                          {'KW':{'$regex':'.*[Pp]450.*'}},
                          {'DR':{'$regex':'.*[Pp]450.*'}},
                          {'DE':{'$regex':'.*[Pp]450.*'}},
                          {'FT':{'$regex':'.*[Hh]eme.*'}},
                          {'CC':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                          {'KW':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                          {'KW':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                          {'DR':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                          {'DE':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                          {'FT':{'$regex':'EC=1\.14\.[0-9]+\.[0-9]+.*'}},
                           ]})

    for doc in qr_p450:
        ec_nums =   [re.findall('(EC=1\.14\.[0-9]+\.[0-9]+)', str(i))\
                        for i in doc.values()]
        ec_nums = list(set([i[0] for i in ec_nums if len(i) != 0]))
        ac = doc['AC']
        if len(ec_nums) > 0:
            doc['ec'] = ec_nums
        if len(list(db['p450s'].find({'AC':doc['AC']}).limit(1))) == 0:
            db['p450s'].insert_one(doc)
        else:
            db['p450s'].update_one({'_id':doc['_id']}, {'$set':doc})

if __name__ == '__main__':
    main(sys.argv[1:])
