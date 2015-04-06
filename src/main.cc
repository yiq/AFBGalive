#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include "htslib/hts.h"
#include "htslib/vcf.h"
#include "RedisKVStore.h"

int compress_output = 0;

int main(int argc, char **argv){

	int c;
	while((c = getopt(argc, argv, "z")) != -1) {
		switch (c) {
			case 'z':
				compress_output = 1;
				break;
			default:
				break;
		}
	}
	

	//YiCppLib::RedisKVStore::pointer rStore(new YiCppLib::RedisKVStore("/var/run/redis.sock"));
	YiCppLib::RedisKVStore::pointer rStore(new YiCppLib::RedisKVStore("127.0.0.1", 6379));

	htsFile *in_fp = hts_open("-", "rz");
	htsFile *out_fp;
	if(compress_output) {
		out_fp = hts_open("-", "wz");
	}
	else {
		out_fp = hts_open("-", "w");
	}

	bcf_hdr_t *hdr = vcf_hdr_read(in_fp);
	bcf_hdr_t *out_hdr = bcf_hdr_dup(hdr);

	int ret;

	ret = bcf_hdr_append(out_hdr, "##INFO=<ID=BGAF_1KG,Number=A,Type=Float,Description=\"Background Allele Frequency in the 1000 Genomes Project\">");
	ret = bcf_hdr_append(out_hdr, "##INFO=<ID=BGAF_EXAC,Number=A,Type=Float,Description=\"Background Allele Frequency in the Exome Aggregation Consortium (ExAC) Project\">");
	
	bcf_hdr_write(out_fp, out_hdr);

	bcf1_t *rec = bcf_init1();

	while(vcf_read1(in_fp, hdr, rec) >= 0) {
		int32_t chrom = rec->rid;
		int32_t pos = rec->pos;

		bcf_unpack(rec, BCF_UN_ALL);

		int nafs = rec->n_allele-1;

		float *afs_1kg = (float *)malloc(sizeof(float) * nafs);
		float *afs_exac = (float *)malloc(sizeof(float) * nafs);


		for(size_t i=1; i<rec->n_allele; i++) {
			std::string key = std::to_string(chrom) + ":" + \
							  std::to_string(pos) + ":" + \
							  rec->d.allele[0] + ":" + \
							  rec->d.allele[i];
			
			std::string AFStr1KG = rStore->stringValueForKeyInNamespace(key, "1KG");
			std::string AFStrExac = rStore->stringValueForKeyInNamespace(key, "EXAC");

			afs_1kg[i-1] = AFStr1KG != "" ? std::stof(AFStr1KG) : 0;
			afs_exac[i-1] = AFStrExac != "" ? std::stof(AFStrExac) : 0;
		}

		ret = bcf_update_info_float(out_hdr, rec, "BGAF_1KG", afs_1kg, nafs);
		ret = bcf_update_info_float(out_hdr, rec, "BGAF_EXAC", afs_exac, nafs);

		vcf_write1(out_fp, out_hdr, rec);

		free(afs_1kg);
		free(afs_exac);

	}

	bcf_hdr_destroy(out_hdr);
	bcf_hdr_destroy(hdr);
	hts_close(in_fp);
	hts_close(out_fp);
}
