#!/bin/sh

MONGODB_ADDR="mongodb://localhost:27017/uniprot"


if !  test -f "BindingDB_All_2022m4.tsv.zip" ; then
	curl 'https://www.bindingdb.org/bind/downloads/BindingDB_All_2022m4.tsv.zip' \
	  -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.9' \
	  -H 'Accept-Language: en-GB,en-US;q=0.9,en;q=0.8' \
	  -H 'Connection: keep-alive' \
	  -H 'Cookie: JSESSIONID=90E4A6EF868539C0EB657EA50B9779DA' \
	  -H 'Referer: https://www.bindingdb.org/bind/chemsearch/marvin/SDFdownload.jsp?download_file=/bind/downloads/BindingDB_All_2022m4.tsv.zip' \
	  -H 'Sec-Fetch-Dest: document' \
	  -H 'Sec-Fetch-Mode: navigate' \
	  -H 'Sec-Fetch-Site: same-origin' \
	  -H 'Sec-Fetch-User: ?1' \
	  -H 'Upgrade-Insecure-Requests: 1' \
	  -H 'User-Agent: Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/101.0.4951.64 Safari/537.36' \
	  -H 'sec-ch-ua: " Not A;Brand";v="99", "Chromium";v="101"' \
	  -H 'sec-ch-ua-mobile: ?0' \
	  -H 'sec-ch-ua-platform: "Linux"' \
	  --compressed > BindingDB_All_2022m4.tsv.zip
fi

#unzip BindingDB_All_2022m4.tsv.zip
zcat -dcf BindingDB_All_2022m4.tsv.zip > BindingDB_All.tsv
source ../venv/bin/activate
./bindingdb.py BindingDB_All.tsv | mongoimport -v $MONGODB_ADDR -c bindingdb
