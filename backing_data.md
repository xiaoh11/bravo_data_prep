# BRAVO API Backing Data.

Data required by the manage.py commands 

## gene

#### Genenames: HUGO Gene Nomenclature Commitee (HGNC)
Can be obtained from [genenames custom downloads](https://www.genenames.org/download/custom/)
Check the Approved symbol, Approved name, Previous symbols, Alias symbols, Ensembl gene ID boxes.
In advanced filtering section, use this where clause in the text box to cut down on the number of results to illustrate this example.
```
gd_app_sym IN ('HBB', 'HBD', 'HBM', 'HBZ', 'HBA1', 'HBA2' ,'HBE1' ,'HBG1' ,'HBG2' ,'HBQ1' ,'HBAP1' ,'HBBP1' ,'HBXP1')
```
The above query can be [url encoded](https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&order_by=gd_app_sym_sort&format=text&where=gd_app_sym%20IN%20(%27HBB%27,%20%27HBD%27,%20%27HBM%27,%20%27HBZ%27,%20%27HBA1%27,%20%27HBA2%27%20,%27HBE1%27%20,%27HBG1%27%20,%27HBG2%27%20,%27HBQ1%27%20,%27HBAP1%27%20,%27HBBP1%27%20,%27HBXP1%27)&submit=submit)

The genenames set from all the chromosomes with the above columns is currently 2.8M
 and can be obtained using curl.
```sh
curl 'https://www.genenames.org/cgi-bin/download/custom?col=gd_app_sym&col=gd_app_name&col=gd_prev_sym&col=gd_aliases&col=gd_pub_ensembl_id&status=Approved&order_by=gd_app_sym_sort&format=text&submit=submit' > hgcn_custom_results.txt
```

Tab delim results
```tab
Approved symbol	Approved name	Previous symbols	Alias symbols	Ensembl gene ID
HBB	hemoglobin subunit beta		CD113t-C, beta-globin	ENSG00000244734
HBBP1	hemoglobin subunit beta pseudogene 1		HBH1, HBHP	ENSG00000229988
HBD	hemoglobin subunit delta			ENSG00000223609
HBE1	hemoglobin subunit epsilon 1		HBE	ENSG00000213931
HBG1	hemoglobin subunit gamma 1		HBG-T2	ENSG00000213934
HBG2	hemoglobin subunit gamma 2		HBG-T1	ENSG00000196565
```

Headers need to be rewritten to match expected column names and the results gzipped
```sh
echo -e "symbol\tname\talias_symbol\tprev_symbol\talias_symbol" > hgcn_genenames.txt
tail -n +2 hgcn_custom_result.txt >> hgcn_genenames.txt
gzip hgcn_genenames.txt
```

#### Canonical Transcripts: Ensembl ID to Ensembl Transcript ID mapping
- No headers
- Two columns
- Whitespace delimited
- Order: gene\_id transcript\_id
- Gzipped

Sourced from Ensembl's [Biomart tool](https://www.ensembl.org/info/data/biomart/)
The following is the query xml for the Gene stable ID and Transcript stable ID of GRCh38.p13,
and the curl command to retrieve the results from the REST endpoint.
```sh
read -r -d '' RAW_QUERY <<-HEREDOC 
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName = "default" formatter = "TSV" header = "0" uniqueRows = "0" count = "" datasetConfigVersion = "0.6" >
<Dataset name = "hsapiens_gene_ensembl" interface = "default" >
<Attribute name = "ensembl_gene_id" />
<Attribute name = "ensembl_transcript_id" />
</Dataset>
</Query>
HEREDOC

QUERY=$(echo "$RAW_QUERY" | tr -d '\n')

curl "http://www.ensembl.org/biomart/martservice" \
  --data-urlencode "query=${QUERY}" > canonical_transcripts.txt

QUERY_URL='http://useast.ensembl.org/biomart/martview/57c40cc27b63f390d259d9fd4894919f?VIRTUALSCHEMANAME=default&ATTRIBUTES=hsapiens_gene_ensembl.default.feature_page.ensembl_gene_id|hsapiens_gene_ensembl.default.feature_page.ensembl_transcript_id&FILTERS=&VISIBLEPANEL=attributepanel'
curl "${QUERY_URL}" > canonical_transcripts.txt
gzip canonical_transcripts.txt
```



